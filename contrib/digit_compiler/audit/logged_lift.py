#!/usr/bin/env python3
"""An explicit, full-space reversible lift of a finite multitape table.

This module is deliberately NOT a spatial single-head implementation. A
collision may inspect the symbols at the k distinct work heads. Every such
collision has a finite table; there is no integer arithmetic in the machine
being lifted. Python dictionaries below represent finite transition tables
and tapes in the test simulator, not a physical addressed-memory instruction.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product


HISTORY_BLANK = "_"


@dataclass(frozen=True)
class Rule:
    source: str
    read: tuple[str, ...]
    target: str
    write: tuple[str, ...]
    move: tuple[int, ...]


class Tape:
    def __init__(self, blank="_", cells=None, head=0):
        self.blank = blank
        self.cells = {i: x for i, x in (cells or {}).items() if x != blank}
        self.head = head

    def read(self):
        return self.cells.get(self.head, self.blank)

    def write(self, x):
        if x == self.blank:
            self.cells.pop(self.head, None)
        else:
            self.cells[self.head] = x

    def snapshot(self):
        return self.head, tuple(sorted(self.cells.items()))


class LoggedLift:
    """Four permutation primitives, with exact inverse on all configurations.

    C: swap each legal (q, old symbols, history blank) with its uniquely
       labelled (transition-control, new symbols, transition-record).
    M: move every work head by the displacement stored in its control.
    E: using the history record, swap transition-control and target-control.
    H: shift the history head right by one cell.

    Forward macro = H E M C. Inverse macro = C M^-1 E H^-1.
    """

    def __init__(self, rules, alphabets, states=None):
        self.rules = tuple(rules)
        self.alphabets = tuple(tuple(x) for x in alphabets)
        self.k = len(self.alphabets)
        self.states = tuple(sorted(set(states or ()) | {
            r.source for r in self.rules} | {r.target for r in self.rules}))
        self.delta = {}
        for i, r in enumerate(self.rules):
            assert len(r.read) == len(r.write) == len(r.move) == self.k
            assert all(d in (-1, 0, 1) for d in r.move)
            assert all(a in g and b in g for a, b, g in zip(r.read, r.write, self.alphabets))
            assert (r.source, r.read) not in self.delta, "Source table must be deterministic"
            self.delta[r.source, r.read] = i
        self.controls = tuple(("Q", q) for q in self.states) + tuple(
            ("P", i) for i in range(len(self.rules)))
        self.history_alphabet = (HISTORY_BLANK,) + tuple(range(len(self.rules)))

    def collision(self, control, symbols, history):
        """C, including its explicit action on every malformed input."""
        if control[0] == "Q" and history == HISTORY_BLANK:
            i = self.delta.get((control[1], symbols))
            if i is not None:
                r = self.rules[i]
                return ("P", i), r.write, i
        elif control[0] == "P":
            i = control[1]
            r = self.rules[i]
            if history == i and symbols == r.write:
                return ("Q", r.source), r.read, HISTORY_BLANK
        return control, symbols, history

    def erase_control(self, control, history):
        """E: one control transposition per nonblank history symbol."""
        if history == HISTORY_BLANK:
            return control
        i = history
        q = ("Q", self.rules[i].target)
        p = ("P", i)
        if control == p:
            return q
        if control == q:
            return p
        return control

    def _collision_on_tapes(self, control, tapes, history):
        control, symbols, h = self.collision(control, tuple(t.read() for t in tapes), history.read())
        for t, a in zip(tapes, symbols):
            t.write(a)
        history.write(h)
        return control

    def _move(self, control, tapes, sign=1):
        if control[0] == "P":
            for tape, displacement in zip(tapes, self.rules[control[1]].move):
                tape.head += sign * displacement

    def forward(self, control, tapes, history):
        assert len(tapes) == self.k
        control = self._collision_on_tapes(control, tapes, history)
        self._move(control, tapes)
        control = self.erase_control(control, history.read())
        history.head += 1
        return control

    def backward(self, control, tapes, history):
        assert len(tapes) == self.k
        history.head -= 1
        control = self.erase_control(control, history.read())
        self._move(control, tapes, -1)
        return self._collision_on_tapes(control, tapes, history)

    def finite_tables(self):
        """Enumerate the complete C and E tables, not merely valid paths."""
        collision = []
        erasure = []
        for control, symbols, history in product(
                self.controls, product(*self.alphabets), self.history_alphabet):
            x = control, tuple(symbols), history
            collision.append((x, self.collision(*x)))
        for control, history in product(self.controls, self.history_alphabet):
            erasure.append(((control, history), (self.erase_control(control, history), history)))
        return collision, erasure

    def manifest(self):
        return {
            "source_states": len(self.states),
            "source_rules": len(self.rules),
            "source_work_tapes": self.k,
            "source_work_alphabet_sizes": [len(a) for a in self.alphabets],
            "lift_control_states": len(self.controls),
            "history_alphabet_size": len(self.history_alphabet),
            "permutation_primitives_per_source_step": 4,
            "history_cells_per_source_step": 1,
            "scope": "Finite multitape permutations, not a one-head nearest-neighbor controller",
        }


@dataclass(frozen=True)
class SerialRule:
    source: str
    tape: int
    read: str
    target: str
    write: str
    move: int


class SerialLoggedLift(LoggedLift):
    """The same lift without Cartesian expansion of untouched tape symbols.

    Each source state selects one work tape, and each source rule reads/writes
    just that tape. This is the interface used by the serial-memory routines.
    All other work-head symbols are literal spectator identities.
    """

    def __init__(self, rules, alphabets, states=None):
        self.rules = tuple(rules)
        self.alphabets = tuple(tuple(x) for x in alphabets)
        self.k = len(self.alphabets)
        self.states = tuple(sorted(set(states or ()) | {
            r.source for r in self.rules} | {r.target for r in self.rules}))
        self.delta = {}
        self.active = {}
        for i, r in enumerate(self.rules):
            assert 0 <= r.tape < self.k and r.move in (-1, 0, 1)
            assert r.read in self.alphabets[r.tape] and r.write in self.alphabets[r.tape]
            assert self.active.setdefault(r.source, r.tape) == r.tape
            assert (r.source, r.read) not in self.delta
            self.delta[r.source, r.read] = i
        self.controls = tuple(("Q", q) for q in self.states) + tuple(
            ("P", i) for i in range(len(self.rules)))
        self.history_alphabet = (HISTORY_BLANK,) + tuple(range(len(self.rules)))

    def collision(self, control, symbols, history):
        symbols = list(symbols)
        if control[0] == "Q" and history == HISTORY_BLANK:
            q = control[1]
            if q in self.active:
                i = self.delta.get((q, symbols[self.active[q]]))
                if i is not None:
                    r = self.rules[i]
                    symbols[r.tape] = r.write
                    return ("P", i), tuple(symbols), i
        elif control[0] == "P":
            i = control[1]
            r = self.rules[i]
            if history == i and symbols[r.tape] == r.write:
                symbols[r.tape] = r.read
                return ("Q", r.source), tuple(symbols), HISTORY_BLANK
        return control, tuple(symbols), history

    def _move(self, control, tapes, sign=1):
        if control[0] == "P":
            r = self.rules[control[1]]
            tapes[r.tape].head += sign * r.move

    def supported_collision_pairs(self):
        """Complete nonidentity C action, with spectator tapes implicit."""
        return [
            {"left": {"control": ["Q", r.source], "active_tape": r.tape,
                      "work_symbol": r.read, "history_symbol": HISTORY_BLANK},
             "right": {"control": ["P", i], "active_tape": r.tape,
                       "work_symbol": r.write, "history_symbol": i}}
            for i, r in enumerate(self.rules)
        ]

    @classmethod
    def from_specification(cls, specification):
        """Adapter for the root/ALU asynchronous finite-table specification."""
        rules = []
        for (q, a), (target, b, move) in sorted(specification["transitions"].items()):
            rules.append(SerialRule(q, specification["active_tape"][q], a,
                                    target, b, move))
        return cls(rules, [specification["alphabet"]] * specification["tape_count"],
                   specification.get("halting_states", ()))
