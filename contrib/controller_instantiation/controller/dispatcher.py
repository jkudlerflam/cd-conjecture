#!/usr/bin/env python3
"""A finite, locally reversible interpreter for a preloaded quantum gate tape.

This is an explicit compiler *interface*, not a length-to-program generator.
The control reads only neutral classical symbols.  Quantum gates are supplied
as whole-space two-site unitaries and may depend on q without changing the
number of classical states.  No third-party packages are needed.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Callable, Iterable, Sequence


DEFAULT_ORACLES = (
    "SWAP", "X", "H", "Z", "CNOT", "CV", "CV_DAG", "TYPE",
    "SHIFT", "SHIFT_DAG", "QUERY", "LOAD", "EXCHANGE", "ROT",
    "ROT_DAG", "CH", "INIT",
)
MOVES = ("MOVE_R", "MOVE_L")
STATIC = ("START", "END", "DATA")


@dataclass(frozen=True, order=True)
class Symbol:
    program: str
    program_mark: int = 0
    data_mark: int = 0


@dataclass(frozen=True)
class Control:
    name: str
    direction: int
    oracle: str | None = None
    adjoint: bool = False


@dataclass
class Configuration:
    tape: list[Symbol]
    position: int
    control: str

    def snapshot(self) -> tuple:
        return self.position, self.control, tuple(self.tape)


def name(phase: str, sign: int, command: str | None = None) -> str:
    return f"{phase}:{'+' if sign == 1 else '-'}" + (
        "" if command is None else ":" + command
    )


class Dispatcher:
    """Generate an actual finite collision permutation and oracle move table.

    Each control label has one fixed movement direction.  A machine step is:
    move the head one site, applying that label's two-site oracle if any;
    then apply the collision permutation to (control, arrived-at symbol).

    A whole-space unitary implementation uses fresh/updated head phases.
    The on-site collision exchanges fresh and updated through C and C^-1.
    The ordered-sweep G blocks are those of SWEEP_INTERFACE.md.
    """

    def __init__(self, oracles: Sequence[str] = DEFAULT_ORACLES):
        if not oracles or len(set(oracles)) != len(oracles):
            raise ValueError("Oracle names must be nonempty and distinct")
        if set(oracles) & set(MOVES + STATIC):
            raise ValueError("Reserved names used as oracle labels")
        self.oracles = tuple(oracles)
        self.commands = self.oracles + MOVES
        self.symbols = tuple(
            Symbol(p, a, b)
            for p in self.commands + STATIC
            for a in range(2)
            for b in range(2)
        )
        self.symbol_index = {s: i for i, s in enumerate(self.symbols)}
        controls: list[Control] = []
        for sign in (1, -1):
            controls.append(Control(name("ADV", sign), sign))
            for command in self.commands:
                controls.append(Control(name("GO", sign, command), 1))
                controls.append(Control(name("RET", sign, command), -1))
                if command in self.oracles:
                    controls.append(Control(
                        name("EXEC", sign, command), 1, command, sign == -1
                    ))
                    controls.append(Control(name("BACK", sign, command), -1))
                else:
                    direction = sign * (1 if command == "MOVE_R" else -1)
                    controls.append(Control(name("SHIFT", sign, command), direction))
        self.controls = tuple(controls)
        self.control_index = {c.name: i for i, c in enumerate(controls)}
        self.control_by_name = {c.name: c for c in controls}
        self.B = len(self.symbols)
        self.K = len(self.controls)
        self.local_factor = (1 + 2 * self.K) * self.B
        self.partial: dict[int, int] = {}
        self._ranges: dict[int, int] = {}
        self._build_partial()
        self.collision = self._complete_permutation()
        self.inverse = [0] * len(self.collision)
        for src, dst in enumerate(self.collision):
            self.inverse[dst] = src

    def encode(self, control: str, symbol: Symbol) -> int:
        return self.control_index[control] * self.B + self.symbol_index[symbol]

    def decode(self, value: int) -> tuple[str, Symbol]:
        c, s = divmod(value, self.B)
        return self.controls[c].name, self.symbols[s]

    def _add(self, q: str, a: Symbol, qp: str, b: Symbol) -> None:
        src, dst = self.encode(q, a), self.encode(qp, b)
        if src in self.partial:
            raise AssertionError(f"Collision domain overlap: {(q, a)}")
        if dst in self._ranges:
            raise AssertionError(
                f"Collision range overlap: {(qp, b)} from "
                f"{self.decode(self._ranges[dst])} and {(q, a)}"
            )
        self.partial[src] = dst
        self._ranges[dst] = src

    def _build_partial(self) -> None:
        for sign in (1, -1):
            advance = name("ADV", sign)
            for command in self.commands:
                go = name("GO", sign, command)
                ret = name("RET", sign, command)
                # ADV arrives at an instruction.  Copy its opcode into the
                # neutral finite control and mark this program location.
                self._add(advance, Symbol(command), go, Symbol(command, 1, 0))
                # The fetched program marker is the only range distinction
                # needed between entry into GO and its scanning self-loop.
                for p in self.commands + STATIC:
                    self._add(go, Symbol(p), go, Symbol(p))
                    self._add(ret, Symbol(p), ret, Symbol(p))
                # At the remembered instruction, erase the opcode using the
                # read-only program symbol; clear P and advance it by one.
                self._add(ret, Symbol(command, 1, 0), advance, Symbol(command))
                if command in self.oracles:
                    execute = name("EXEC", sign, command)
                    back = name("BACK", sign, command)
                    self._add(go, Symbol("DATA", 0, 1), execute, Symbol("DATA", 0, 1))
                    self._add(execute, Symbol("DATA"), back, Symbol("DATA"))
                    self._add(back, Symbol("DATA", 0, 1), ret, Symbol("DATA", 0, 1))
                else:
                    shift = name("SHIFT", sign, command)
                    self._add(go, Symbol("DATA", 0, 1), shift, Symbol("DATA"))
                    self._add(shift, Symbol("DATA"), ret, Symbol("DATA", 0, 1))
            # Endpoint bounces are entries of a partial injection.  They are
            # not one-way updates with the output otherwise fixed.  The
            # completion below supplies every remaining inverse case.
            endpoint = "END" if sign == 1 else "START"
            self._add(advance, Symbol(endpoint), name("ADV", -sign), Symbol(endpoint))

    def _complete_permutation(self) -> list[int]:
        size = self.K * self.B
        result = [-1] * size
        for src, dst in self.partial.items():
            result[src] = dst
        unused_domain = set(range(size)) - self.partial.keys()
        unused_range = set(range(size)) - self._ranges.keys()
        # Complete separately within each read-only program-symbol sector.
        # Thus malformed inputs cannot change the program either.  Each
        # prescribed entry already preserves this symbol, so the numbers
        # of unmatched domain/range elements in each sector are equal.
        for program in self.commands + STATIC:
            domain = {x for x in unused_domain if self.decode(x)[1].program == program}
            image = {x for x in unused_range if self.decode(x)[1].program == program}
            assert len(domain) == len(image)
            common = domain & image
            for x in common:
                result[x] = x
            for src, dst in zip(sorted(domain - common), sorted(image - common)):
                result[src] = dst
        assert sorted(result) == list(range(size))
        return result

    def collide(self, control: str, symbol: Symbol, inverse: bool = False) -> tuple[str, Symbol]:
        table = self.inverse if inverse else self.collision
        return self.decode(table[self.encode(control, symbol)])

    def make_configuration(self, program: Sequence[str], width: int) -> Configuration:
        if width < 2:
            raise ValueError("At least two data sites are required")
        if any(command not in self.commands for command in program):
            raise ValueError("Unknown program command")
        tape = [Symbol("START")]
        tape.extend(Symbol(command) for command in program)
        tape.append(Symbol("END"))
        tape.extend(Symbol("DATA", 0, int(i == 0)) for i in range(width))
        return Configuration(tape, 0, name("ADV", 1))

    def step(
        self,
        config: Configuration,
        callback: Callable[[str, int, bool], None] | None = None,
        require_prescribed: bool = True,
    ) -> tuple[str, int, bool] | None:
        """One move plus collision; optional callback receives physical index.

        Requiring a prescribed collision is useful for checking promised
        initialized programs.  The completed collision is defined on all
        inputs regardless; unused completion rules have no intended program
        semantics and are not a substitute for initialization penalties.
        """
        control = self.control_by_name[config.control]
        target = config.position + control.direction
        if not 0 <= target < len(config.tape):
            raise ValueError("Program attempted to leave the physical tape")
        event = None
        if control.oracle is not None:
            left = min(config.position, target)
            event = (control.oracle, left, control.adjoint)
            if callback:
                callback(*event)
        config.position = target
        inp = self.encode(config.control, config.tape[target])
        if require_prescribed and inp not in self.partial:
            raise ValueError("Program left the prescribed initialized sector")
        config.control, config.tape[target] = self.decode(self.collision[inp])
        return event

    def run_cycle(
        self, program: Sequence[str], width: int,
        callback: Callable[[str, int, bool], None] | None = None,
    ) -> dict:
        config = self.make_configuration(program, width)
        initial = config.snapshot()
        events = []
        snapshots = []
        midpoint = None
        # Every instruction costs at most 2n+4 moves and occurs twice.
        limit = 2 * (len(program) + 1) * (2 * len(config.tape) + 6)
        for tick in range(1, limit + 1):
            event = self.step(config, callback)
            if event:
                events.append(event)
            snapshots.append(config.snapshot())
            if config.control == name("ADV", -1) and config.position == len(program) + 1:
                if midpoint is not None:
                    raise AssertionError("Repeated forward endpoint")
                midpoint = tick
            if config.snapshot() == initial:
                assert midpoint is not None
                return {
                    "steps": tick,
                    "midpoint": midpoint,
                    "events": events,
                    "snapshots": snapshots,
                    "final": config,
                }
        raise AssertionError("No closed forward/inverse cycle within proved bound")

    def manifest(self) -> dict:
        return {
            "schema": "cd-neutral-dispatcher-v1",
            "scope": "preloaded program interpreter; no length-to-program generator",
            "oracles": list(self.oracles),
            "neutral_alphabet_size_B": self.B,
            "direction_fixed_controls_K": self.K,
            "head_including_single_site_factor": self.local_factor,
            "collision_size": self.K * self.B,
            "prescribed_collision_entries": len(self.partial),
            "neutral_symbols": [
                [s.program, s.program_mark, s.data_mark] for s in self.symbols
            ],
            "controls": [
                {"name": c.name, "direction": c.direction,
                 "oracle": c.oracle, "adjoint": c.adjoint}
                for c in self.controls
            ],
            "collision_permutation": self.collision,
        }


def iter_lower_gate_stream(
    stream: Iterable[tuple[str, int] | dict], width: int,
    allowed_oracles: Sequence[str] = DEFAULT_ORACLES,
) -> Iterable[str]:
    """Lower (opcode, left-address) to a fixed-alphabet cursor instruction tape.

    Integer addresses are compiled into repeated cursor moves; they never
    become controller states.  This expands by at most width per gate.
    Only a polynomial-in-p circuit should be expanded this way.  A physical
    length-n endpoint walk must use the separate compact marked-walk routine.
    """
    cursor = 0
    for instruction in stream:
        if isinstance(instruction, dict):
            opcode, left = instruction["opcode"], instruction["left"]
        else:
            opcode, left = instruction
        if opcode not in allowed_oracles:
            raise ValueError(f"Unknown oracle: {opcode}")
        if not isinstance(left, int) or not 0 <= left < width - 1:
            raise ValueError(f"Invalid left operand address: {left}")
        direction = "MOVE_R" if left >= cursor else "MOVE_L"
        for _ in range(abs(left - cursor)):
            yield direction
        yield opcode
        cursor = left


def lower_gate_stream(
    stream: Iterable[tuple[str, int] | dict], width: int,
    allowed_oracles: Sequence[str] = DEFAULT_ORACLES,
) -> list[str]:
    """Materialize the streaming lowerer when the program is reasonably small."""
    return list(iter_lower_gate_stream(stream, width, allowed_oracles))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    dispatcher = Dispatcher()
    manifest = dispatcher.manifest()
    if args.output:
        args.output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k not in (
        "neutral_symbols", "controls", "collision_permutation"
    )}, indent=2))


if __name__ == "__main__":
    main()
