"""Autonomous reversible ONE-head simulation of any finite one-head table.

The history alphabet stores a fixed source transition label (q, read-symbol),
never a runtime integer or pointer. The local table is factorized to avoid an
unnecessary quadratic-size materialization. local_rule sees only one control
and its scanned symbol. inverse_rule supplies the complete partial inverse.
"""
from dataclasses import dataclass


LEFT, SEPARATOR, BLANK, FRONTIER, RIGHT = range(5)


@dataclass(frozen=True)
class Work:
    symbol: object
    marked: bool = False


@dataclass(frozen=True)
class Record:
    state: object
    symbol: object


@dataclass(frozen=True)
class History:
    record: Record


class Source:
    def __init__(self, states, alphabet, halting, rule):
        self.states = tuple(states)
        self.alphabet = tuple(alphabet)
        self.halting = frozenset(halting)
        self._rule = rule

    def rule(self, state, symbol):
        if state in self.halting:
            return None
        return self._rule(state, symbol)

    @classmethod
    def from_dictionary(cls, alphabet, transitions, halting):
        states = set(halting) | {q for q, _ in transitions} | {r for r, _, _ in transitions.values()}
        return cls(sorted(states, key=repr), alphabet, halting,
                   lambda q, a: transitions.get((q, a)))

    @classmethod
    def from_spatial_compiler(cls, compiler):
        halting = {q for q in compiler.controls if q[0] in ("HALT", "REJECT")}
        return cls(sorted(compiler.controls, key=repr), compiler.alphabet(), halting,
                   compiler.local_rule)


class Compiler:
    def __init__(self, source):
        self.source = source

    @staticmethod
    def plain(cell):
        return ((isinstance(cell, Work) and not cell.marked)
                or isinstance(cell, History) or cell in (SEPARATOR, BLANK))

    def event(self, record):
        return self.source.rule(record.state, record.symbol)

    def local_rule(self, control, cell):
        """One symbol read/write and at most one physical head move."""
        phase, label = control
        if phase == "F":
            if not isinstance(cell, Work) or not cell.marked:
                return None
            t = Record(label, cell.symbol)
            rule = self.event(t)
            if rule is None:
                return None
            _, b, d = rule
            if d:
                return ("MOVE", t), Work(b), d
            return ("DEPART", t), Work(b, True), 0
        if phase in ("MOVE", "DEPART", "SCAN"):
            t = label
            rule = self.event(t)
            if rule is None:
                return None
            r, _, d = rule
            if phase == "MOVE":
                if d and isinstance(cell, Work) and not cell.marked:
                    return ("DEPART", t), Work(cell.symbol, True), 0
            elif phase == "DEPART":
                if isinstance(cell, Work) and cell.marked:
                    return ("SCAN", t), cell, 1
            elif self.plain(cell):
                return control, cell, 1
            elif cell == FRONTIER:
                return ("APPEND", r), History(t), 1
            return None
        if phase == "APPEND":
            if cell == BLANK:
                return ("RETSTART", label), FRONTIER, 0
        elif phase == "RETSTART":
            if cell == FRONTIER:
                return ("RETURN", label), cell, -1
        elif phase == "RETURN":
            if self.plain(cell):
                return control, cell, -1
            if isinstance(cell, Work) and cell.marked:
                return ("F", label), cell, 0
        return None

    def incoming_direction(self, target):
        phase, label = target
        if phase == "MOVE":
            rule = self.event(label)
            return None if rule is None or rule[2] == 0 else rule[2]
        return {"F": 0, "DEPART": 0, "SCAN": 1, "APPEND": 1,
                "RETSTART": 0, "RETURN": -1}.get(phase)

    def inverse_rule(self, target, written):
        """Returns predecessor control, old symbol, forward move, or None.

        No history stack is popped by the host. The inverse reads the symbol
        at the unique preceding physical position determined by target control.
        """
        phase, label = target
        if phase == "F":
            if isinstance(written, Work) and written.marked:
                return ("RETURN", label), written, 0
        elif phase in ("MOVE", "DEPART", "SCAN"):
            t = label
            rule = self.event(t)
            if rule is None:
                return None
            _, b, d = rule
            if phase == "MOVE":
                if d and written == Work(b):
                    return ("F", t.state), Work(t.symbol, True), d
            elif phase == "DEPART":
                if d and isinstance(written, Work) and written.marked:
                    return ("MOVE", t), Work(written.symbol), 0
                if not d and written == Work(b, True):
                    return ("F", t.state), Work(t.symbol, True), 0
            elif isinstance(written, Work) and written.marked:
                return ("DEPART", t), written, 1
            elif self.plain(written):
                return ("SCAN", t), written, 1
        elif phase == "APPEND":
            if isinstance(written, History):
                rule = self.event(written.record)
                if rule is not None and rule[0] == label:
                    return ("SCAN", written.record), FRONTIER, 1
        elif phase == "RETSTART":
            if written == FRONTIER:
                return ("APPEND", label), BLANK, 0
        elif phase == "RETURN":
            if written == FRONTIER:
                return ("RETSTART", label), written, -1
            if self.plain(written):
                return ("RETURN", label), written, -1
        return None

    def materialize(self):
        """Exact full finite partial table; intended only for small fixtures."""
        events = [Record(q, a) for q in self.source.states for a in self.source.alphabet
                  if self.source.rule(q, a) is not None]
        controls = [("F", q) for q in self.source.states]
        controls += [(phase, t) for t in events for phase in ("DEPART", "SCAN")]
        controls += [("MOVE", t) for t in events if self.event(t)[2]]
        controls += [(phase, q) for q in self.source.states for phase in ("APPEND", "RETSTART", "RETURN")]
        alphabet = [LEFT, SEPARATOR, BLANK, FRONTIER, RIGHT]
        alphabet += [Work(a, m) for a in self.source.alphabet for m in (False, True)]
        alphabet += [History(t) for t in events]
        rules = {(q, a): out for q in controls for a in alphabet
                 if (out := self.local_rule(q, a)) is not None}
        return controls, alphabet, rules

    def manifest(self):
        q, g = len(self.source.states), len(self.source.alphabet)
        m = (q-len(self.source.halting))*g
        return {"source_states": q, "source_alphabet": g,
                "source_transition_count_upper_bound": m,
                "reversible_controls_upper_bound": 4*q+3*m,
                "reversible_alphabet_upper_bound": 2*g+m+5,
                "history_cells_per_source_step": 1,
                "microstep_bound_per_source_step": "2*S+8; S is the entire allocated work and history interval",
                "scope": "Finite single-head reversible partial transition table; initialized work/history format promised"}


class Runner:
    def __init__(self, compiler, state, work, position, history_capacity):
        assert 0 <= position < len(work) and history_capacity >= 1
        self.compiler = compiler
        self.cells = [LEFT]+[Work(a, i == position) for i, a in enumerate(work)]
        self.cells += [SEPARATOR, FRONTIER]+[BLANK]*history_capacity+[RIGHT]
        self.position = position+1
        self.control = ("F", state)
        self.steps = 0
        self.work_length = len(work)

    def snapshot(self):
        return self.control, self.position, tuple(self.cells)

    def step(self):
        out = self.compiler.local_rule(self.control, self.cells[self.position])
        if out is None:
            return False
        self.control, self.cells[self.position], move = out
        self.position += move
        assert 0 <= self.position < len(self.cells)
        self.steps += 1
        return True

    def backward(self):
        direction = self.compiler.incoming_direction(self.control)
        if direction is None:
            return False
        previous = self.position-direction
        if not 0 <= previous < len(self.cells):
            return False
        out = self.compiler.inverse_rule(self.control, self.cells[previous])
        if out is None:
            return False
        self.control, self.cells[previous], move = out
        assert move == direction
        self.position = previous
        self.steps -= 1
        return True

    def source_step(self):
        assert self.control[0] == "F"
        start = self.steps
        while self.step():
            if self.control[0] == "F":
                assert self.steps-start <= 2*len(self.cells)+8
                return self.steps-start
        raise ValueError("Source halted or work/history format/capacity insufficient")

    def decode_work(self):
        return [cell.symbol for cell in self.cells[1:self.work_length+1]]
