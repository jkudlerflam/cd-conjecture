"""Explicit one-head spatial lowering of an asynchronous multitape table.

Each virtual position is a block of k physical cells, tagged 0,...,k-1.
Each cell carries one source symbol and one virtual-head marker bit. Only
the physical head's current cell is inspected or changed by a microstep.
This is a DETERMINISTIC, generally IRREVERSIBLE machine. The separate
multitape reversible lift must not be mistaken for its reversible completion.
"""
from dataclasses import dataclass

from tm import Configuration

LEFT = "LEFT"
RIGHT = "RIGHT"
REJECT = ("REJECT",)


@dataclass(frozen=True)
class Cell:
    tape: int
    symbol: str
    marked: bool


class Compiler:
    def __init__(self, source):
        self.source = source
        self.k = source.tape_count
        self.controls = {( "SCAN", q) for q in source.active_tape}
        self.controls |= {("RETURN", q) for q in source.states}
        self.controls |= {("HALT", q) for q in source.halting_states}
        self.controls.add(REJECT)
        for (q, _), (r, _, d) in source.transitions.items():
            if d:
                i = source.active_tape[q]
                self.controls.update(("MOVE", r, i, d, left) for left in range(self.k))

    def local_rule(self, control, cell):
        """Total finite table evaluator; output (control, symbol, unit move).

        This function has no access to the tape, a virtual-head coordinate,
        a stored integer value, or the length of the physical tape.
        """
        if control == REJECT or control[0] == "HALT":
            return control, cell, 0
        phase, q, *rest = control
        if phase == "SCAN":
            if cell == LEFT:
                return control, cell, 1
            if cell == RIGHT:
                return REJECT, cell, 0
            i = self.source.active_tape[q]
            if cell.tape != i or not cell.marked:
                return control, cell, 1
            rule = self.source.transitions.get((q, cell.symbol))
            if rule is None:
                return REJECT, cell, 0
            r, b, d = rule
            if d == 0:
                return ("RETURN", r), Cell(i, b, True), 0
            return ("MOVE", r, i, d, self.k - 1), Cell(i, b, False), d
        if phase == "MOVE":
            i, direction, remaining = rest
            if cell in (LEFT, RIGHT):
                return REJECT, cell, 0
            if remaining:
                return ("MOVE", q, i, direction, remaining - 1), cell, direction
            if cell.tape != i or cell.marked:
                return REJECT, cell, 0
            return ("RETURN", q), Cell(i, cell.symbol, True), 0
        if phase == "RETURN":
            if cell == LEFT:
                phase = "HALT" if q in self.source.halting_states else "SCAN"
                return (phase, q), cell, 0
            return control, cell, -1
        raise AssertionError(control)

    def alphabet(self):
        return (LEFT, RIGHT) + tuple(Cell(i, a, mark)
                for i in range(self.k) for a in self.source.alphabet for mark in (False, True))

    def manifest(self):
        return {
            "source_tapes": self.k,
            "source_alphabet": len(self.source.alphabet),
            "source_states_including_halts": len(self.source.states),
            "source_rules": len(self.source.transitions),
            "physical_alphabet_size": 2 + 2 * self.k * len(self.source.alphabet),
            "physical_control_states": len(self.controls),
            "full_table_entries": len(self.controls) * (2 + 2 * self.k * len(self.source.alphabet)),
            "physical_cells_for_L_virtual_positions": f"{self.k}*L+2",
            "microstep_bound_per_source_step": "2*S+k+4, S=physical cells including endpoints",
            "scope": "One physical head, deterministic table; not reversible and not yet a quantum gate",
        }


class Runner:
    def __init__(self, compiler, source_config, lo, hi):
        """Boundary adapter. lo..hi is a preallocated virtual interval.

        This adapter is NOT an initialization algorithm for the final core.
        """
        assert lo <= hi and len(source_config.tapes) == compiler.k
        assert all(lo <= h <= hi for h in source_config.heads)
        assert all(lo <= p <= hi for t in source_config.tapes for p in t)
        self.compiler = compiler
        self.lo, self.hi = lo, hi
        s = compiler.source
        self.cells = [LEFT] + [Cell(i, source_config.tapes[i].get(p, s.blank),
                                   source_config.heads[i] == p)
                for p in range(lo, hi + 1) for i in range(compiler.k)] + [RIGHT]
        self.position = 0
        phase = "HALT" if source_config.state in s.halting_states else "SCAN"
        self.control = (phase, source_config.state)
        self.steps = 0

    def step(self):
        if self.control == REJECT or self.control[0] == "HALT":
            return False
        # The only dynamic cell access is at the one physical head.
        control, symbol, move = self.compiler.local_rule(self.control, self.cells[self.position])
        self.cells[self.position] = symbol
        self.position += move
        self.control = control
        self.steps += 1
        assert 0 <= self.position < len(self.cells)
        return True

    def macrostep(self):
        assert self.position == 0 and self.control[0] == "SCAN"
        before = self.steps
        limit = 2 * len(self.cells) + self.compiler.k + 4
        while self.steps - before <= limit:
            self.step()
            if self.control == REJECT:
                raise ValueError("One-head machine rejected: malformed tape or insufficient capacity")
            if self.position == 0 and self.control[0] in ("SCAN", "HALT"):
                assert self.steps - before <= limit
                return self.steps - before
        raise AssertionError("Proved sweep bound exceeded")

    def decode(self):
        """Independent boundary comparison only; never called by step()."""
        tapes = [{} for _ in range(self.compiler.k)]
        heads = [None] * self.compiler.k
        for j, cell in enumerate(self.cells[1:-1]):
            p, i = divmod(j, self.compiler.k)
            p += self.lo
            assert cell.tape == i
            if cell.symbol != self.compiler.source.blank:
                tapes[i][p] = cell.symbol
            if cell.marked:
                assert heads[i] is None
                heads[i] = p
        assert all(h is not None for h in heads)
        assert self.position == 0 and self.control[0] in ("SCAN", "HALT")
        return Configuration(self.control[1], tapes, heads)
