"""Independent interpreter for fixed asynchronous multitape transition tables.

Payload arithmetic is absent from step(): it only looks up a finite rule,
writes one symbol, and moves one head by -1, 0, or 1. Integer tape coordinates
and counters are simulator bookkeeping. They are not machine instructions.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class AsyncTable:
    alphabet: tuple
    blank: str
    tape_count: int
    active_tape: dict
    transitions: dict
    halting_states: frozenset
    initial_state: str

    @classmethod
    def from_export(cls, data):
        rules = data["transitions"]
        if not isinstance(rules, dict):
            rows = rules
            rules = {(q, a): (r, b, d) for q, a, r, b, d in rows}
            assert len(rules) == len(rows), "Duplicate transition"
        table = cls(tuple(data["alphabet"]), data["blank"], data["tape_count"],
                    dict(data["active_tape"]), dict(rules),
                    frozenset(data["halting_states"]), data["initial_state"])
        table.validate()
        return table

    @property
    def states(self):
        return set(self.active_tape) | set(self.halting_states)

    def validate(self):
        states = self.states
        assert self.blank in self.alphabet and len(set(self.alphabet)) == len(self.alphabet)
        assert self.tape_count > 0 and self.initial_state in states
        assert not set(self.active_tape) & self.halting_states
        assert all(0 <= i < self.tape_count for i in self.active_tape.values())
        for (q, a), (r, b, d) in self.transitions.items():
            assert q in self.active_tape and r in states
            assert a in self.alphabet and b in self.alphabet and d in (-1, 0, 1)


class Configuration:
    def __init__(self, state, tapes, heads):
        self.state = state
        self.tapes = [dict(t) for t in tapes]
        self.heads = list(heads)
        self.steps = 0

    def snapshot(self):
        return self.state, tuple(tuple(sorted(t.items())) for t in self.tapes), tuple(self.heads)


def step(table, config):
    if config.state in table.halting_states:
        return False
    i = table.active_tape[config.state]
    p = config.heads[i]
    a = config.tapes[i].get(p, table.blank)
    rule = table.transitions.get((config.state, a))
    if rule is None:
        raise ValueError(("Undefined source transition", config.state, a))
    r, b, d = rule
    if b == table.blank:
        config.tapes[i].pop(p, None)
    else:
        config.tapes[i][p] = b
    config.heads[i] += d
    config.state = r
    config.steps += 1
    return True
