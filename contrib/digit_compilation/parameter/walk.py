#!/usr/bin/env python3
"""Constant-state parameter extraction, midpoint, and autonomous cleanup.

The six gates from parameter_oracle.py form a fixed on-site collision
subword; there is no expanded instruction tape or numerical loop counter.
All local F/G gates are full-space permutations, including malformed input.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "autonomous_controller" / "audit"))
from parameter_oracle import Site as OracleSite, alphabet, local, species


MARKS = ("LEFT", "BULK", "RIGHT")
WRITE = ("PARAM", "COPY_DIGIT", "PARAM", "END", "COPY_END", "END")
# (name, direction, movement oracle). PUSH shifts at the departure then
# swaps Sigma; POP swaps Sigma then inverse-shifts at the arrival.
CONTROLS = (
    ("F:LAUNCH", -1, "I"), ("F:START", 1, "I"),
    ("F:OUT", -1, "PUSH"), ("F:BACK", 1, "POP"),
    ("R:LAUNCH", 1, "I"), ("R:START", -1, "I"),
    ("R:OUT", 1, "POP"), ("R:BACK", -1, "PUSH"),
)


@dataclass(frozen=True)
class Site:
    mark: str
    sigma: tuple
    digit: int = 0
    end: int = 0
    control: str | None = None
    phase: int = 0


class ParameterWalk:
    def __init__(self, p, q, m):
        if not 2 <= m or not 0 <= p < 2 ** m or q < 1:
            raise ValueError("Need q>=1, m>=2, and 0<=p<2**m")
        self.p, self.q, self.m = p, q, m
        self.alphabet = tuple(alphabet(q, m))
        self.symbols = tuple((mark, digit, end) for mark in MARKS
                             for digit in range(2) for end in range(2))
        self.controls = tuple(q for q, _, _ in CONTROLS)
        self.direction = {q: d for q, d, _ in CONTROLS}
        self.oracle = {q: gate for q, _, gate in CONTROLS}
        self.control_index = {q: j for j, q in enumerate(self.controls)}
        self.symbol_index = {a: j for j, a in enumerate(self.symbols)}
        self.B, self.K = len(self.symbols), len(self.controls)
        self.local_factor = (1 + 2 * self.K) * self.B
        self.dimension = self.local_factor * len(self.alphabet)
        self.partial, self.images = {}, {}
        self._build()
        self.collision = self._complete()
        self.inverse = [0] * len(self.collision)
        for a, b in enumerate(self.collision):
            self.inverse[b] = a

    def encode(self, control, symbol):
        return self.control_index[control] * self.B + self.symbol_index[symbol]

    def decode(self, index):
        q, a = divmod(index, self.B)
        return self.controls[q], self.symbols[a]

    def _add(self, q, mark, end, target):
        for digit in range(2):
            a = self.encode(q, (mark, digit, end))
            b = self.encode(target, (mark, digit, end))
            assert a not in self.partial, (q, mark, end)
            assert b not in self.images, (target, mark, end)
            self.partial[a], self.images[b] = b, a

    def _build(self):
        # These are literal entries, not rules generated from p,m or n.
        for mark in ("LEFT", "BULK"):
            for end in (0, 1):
                self._add("F:LAUNCH", mark, end, "F:START")
                self._add("R:START", mark, end, "R:LAUNCH")
            self._add("F:OUT", mark, 0, "F:OUT")
            self._add("F:OUT", mark, 1, "F:BACK")
            self._add("F:BACK", mark, 0, "F:BACK")
            self._add("R:OUT", mark, 0, "R:OUT")
            self._add("R:BACK", mark, 1, "R:OUT")
            self._add("R:BACK", mark, 0, "R:BACK")
        self._add("F:START", "RIGHT", 0, "F:OUT")
        self._add("F:START", "RIGHT", 1, "F:BACK")
        self._add("R:OUT", "RIGHT", 0, "R:START")
        self._add("R:BACK", "RIGHT", 1, "R:START")
        # Reversible endpoint bounces, with distinct unused source/range.
        self._add("F:BACK", "RIGHT", 0, "R:BACK")
        self._add("R:LAUNCH", "RIGHT", 0, "F:LAUNCH")

    def _complete(self):
        table = list(range(self.B * self.K))
        for a, b in self.partial.items():
            table[a] = b
        # Complete separately in EACH neutral-symbol sector. Thus all
        # static marks and both dynamic bits are conserved by this C table.
        for symbol in self.symbols:
            sector = {self.encode(q, symbol) for q in self.controls}
            domain = sector - self.partial.keys()
            image = sector - self.images.keys()
            common = domain & image
            for a, b in zip(sorted(domain - common), sorted(image - common)):
                table[a] = b
        assert sorted(table) == list(range(len(table)))
        return table

    def apply_local(self, site, op):
        result = local(op, OracleSite(site.sigma, site.digit, site.end), self.p, self.m)
        return replace(site, sigma=result.sigma, digit=result.digit, end=result.end)

    def write(self, site, reverse=False):
        for op in reversed(WRITE) if reverse else WRITE:
            site = self.apply_local(site, op)
        return site

    @staticmethod
    def write_forward(site):
        return ((site.control == "F:START" and site.mark == "RIGHT") or
                (site.control == "F:OUT" and site.mark != "RIGHT"))

    @staticmethod
    def write_reverse(site):
        return ((site.control == "R:START" and site.mark == "RIGHT") or
                (site.control == "R:OUT" and site.mark != "RIGHT"))

    def collide(self, site, reverse=False, prescribed=False):
        """P = controlled WRITE^-1_R C controlled WRITE_F, or its inverse.

        C is the saved 96-entry neutral permutation. Controlled WRITE only
        depends on a finite control/static mark and does not change either.
        Hence every factor is a permutation on the WHOLE local alphabet.
        """
        if (self.write_reverse(site) if reverse else self.write_forward(site)):
            site = self.write(site, reverse=False)
        source = self.encode(site.control, (site.mark, site.digit, site.end))
        if prescribed and source not in self.partial:
            raise AssertionError("Initialized walk entered a completion-only collision")
        table = self.inverse if reverse else self.collision
        q, (mark, digit, end) = self.decode(table[source])
        site = replace(site, mark=mark, digit=digit, end=end, control=q)
        if (self.write_forward(site) if reverse else self.write_reverse(site)):
            site = self.write(site, reverse=True)
        return site

    def f(self, site):
        """Whole-space on-site involution, fixing no-head sectors."""
        if site.control is None:
            return site
        return replace(self.collide(site, reverse=site.phase == 0), phase=1-site.phase)

    def _transport(self, a, b, on_left, gate, reverse):
        if gate == "I":
            return a, b
        # Forward PUSH: SHIFT departure, SWAP. Forward POP: SWAP,
        # SHIFT_DAG arrival. Reverse branches are their exact adjoints.
        departure, arrival = (0, 1) if on_left else (1, 0)
        pair = [a, b]
        if gate == "PUSH" and not reverse:
            pair[departure] = self.apply_local(pair[departure], "SHIFT")
        elif gate == "POP" and reverse:
            pair[departure] = self.apply_local(pair[departure], "SHIFT")
        pair = [replace(pair[0], sigma=pair[1].sigma),
                replace(pair[1], sigma=pair[0].sigma)]
        if gate == "POP" and not reverse:
            pair[arrival] = self.apply_local(pair[arrival], "SHIFT_DAG")
        elif gate == "PUSH" and reverse:
            pair[arrival] = self.apply_local(pair[arrival], "SHIFT_DAG")
        return tuple(pair)

    def g(self, a, b):
        """Whole-space two-site involution; zero/two heads are fixed."""
        if (a.control is None) == (b.control is None):
            return a, b
        on_left = a.control is not None
        head = a if on_left else b
        direction = self.direction[head.control]
        forward = head.phase == 0 and ((direction == 1) == on_left)
        reverse = head.phase == 1 and ((direction == -1) == on_left)
        if not (forward or reverse):
            return a, b
        a, b = self._transport(a, b, on_left, self.oracle[head.control], reverse)
        if on_left:
            return (replace(a, control=None, phase=0),
                    replace(b, control=head.control, phase=1-head.phase))
        return (replace(a, control=head.control, phase=1-head.phase),
                replace(b, control=None, phase=0))

    def initialize(self, payload):
        """Payload includes its pinned H_(0,0) at RIGHT; n>=m>=2."""
        if len(payload) < self.m or payload[-1] != ("H", 0, 0):
            raise ValueError("Need n>=m and pinned H_(0,0) at RIGHT")
        return tuple(Site("LEFT" if j == 0 else "RIGHT" if j == len(payload)-1 else "BULK",
                          sigma, control="F:LAUNCH" if j == len(payload)-1 else None)
                     for j, sigma in enumerate(payload))

    def micro_sweep(self, initial):
        tape = list(initial)
        for at in range(len(tape)-1):
            tape[at], tape[at+1] = self.g(tape[at], tape[at+1])
            yield tuple(tape), ("G", at)
        for at in reversed(range(len(tape))):
            tape[at] = self.f(tape[at])
            yield tuple(tape), ("F", at)

    def step(self, initial):
        """One prescribed sweep, independently callable at every step."""
        tape = list(initial)
        heads = [j for j, a in enumerate(tape) if a.control is not None]
        assert len(heads) == 1
        at = heads[0]
        assert tape[at].phase == 0
        destination = at + self.direction[tape[at].control]
        if not 0 <= destination < len(tape):
            raise ValueError("Head left segment (initialization contract violated)")
        bond = min(at, destination)
        tape[bond], tape[bond+1] = self.g(tape[bond], tape[bond+1])
        assert tape[destination].phase == 1
        tape[destination] = replace(self.collide(tape[destination], prescribed=True), phase=0)
        return tuple(tape)

    def manifest(self):
        return {
            "schema": "cd-parameter-walk-v1",
            "scope": "Pinned indexed helper to parameter bits at midpoint and exact cleanup; no generator/controller composition",
            "initial_contract": "n>=m>=2; LEFT BULK^(n-2) RIGHT; H_(0,0) at RIGHT; neutral tracks zero; F:LAUNCH fresh head at RIGHT",
            "neutral_alphabet_B": self.B,
            "direction_fixed_controls_K": self.K,
            "single_site_factor_times_Sigma": self.local_factor,
            "prescribed_collision_entries": len(self.partial),
            "midpoint_sweeps": "2m",
            "full_cycle_sweeps": "4m",
            "controls": [{"name": q, "direction": d, "oracle": op} for q, d, op in CONTROLS],
            "neutral_symbols": self.symbols,
            "collision_permutation_C": self.collision,
            "full_collision_P": "controlled WRITE^-1 on reverse output control, after C, after controlled WRITE on forward input control",
            "WRITE": WRITE,
            "WRITE_input_controls": [["F:START", "RIGHT"], ["F:OUT", "LEFT"], ["F:OUT", "BULK"]],
            "WRITE_inverse_output_controls": [["R:START", "RIGHT"], ["R:OUT", "LEFT"], ["R:OUT", "BULK"]],
            "movement_oracles": {"I": [], "PUSH": ["SHIFT departure", "Sigma SWAP"],
                                  "POP": ["Sigma SWAP", "SHIFT_DAG arrival"]},
            "completion": "Within each (mark,digit,end), fix unused common domain/range; match remaining sorted indices",
        }


if __name__ == "__main__":
    output = Path(__file__).with_name("transition_table.json")
    output.write_text(json.dumps(ParameterWalk(2, 2, 2).manifest(), indent=2) + "\n")
    print(output)
