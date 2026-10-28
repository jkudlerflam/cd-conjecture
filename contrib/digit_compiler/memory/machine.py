#!/usr/bin/env python3
"""Finite-state, serial binary-word array access; no indexed array primitive.

Every base transition touches exactly one cell of one tape and moves that tape's
head by at most one. Python dictionaries below represent physical tape cells;
the only dictionary key used by an execution step is that tape's current head.
Address arithmetic is performed by the transition table on the counter tape.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path

ALPHABET = ("_", "^", "0", "1", "#", "$")
TAPES = ("M", "A", "V", "C", "R")


@dataclass(frozen=True)
class Rule:
    tape: str
    read: str
    write: str
    move: int
    target: str


class Table:
    def __init__(self):
        self.rules = {}
        self.active = {}

    def add(self, state, tape, read, target, write=None, move=0):
        assert read in ALPHABET and (write is None or write in ALPHABET)
        assert move in (-1, 0, 1)
        assert self.active.setdefault(state, tape) == tape
        assert (state, read) not in self.rules
        self.rules[state, read] = Rule(tape, read, read if write is None else write, move, target)

    def same(self, state, tape, target, move=0):
        for s in ALPHABET:
            self.add(state, tape, s, target, move=move)

    def rewind(self, state, tape, target):
        self.add(state, tape, "^", target)
        for s in ("0", "1", "#", "$"):
            self.add(state, tape, s, state, move=-1)

    def payload(self):
        tape_names = list(dict.fromkeys(r.tape for r in self.rules.values()))
        states = set(self.active) | {r.target for r in self.rules.values()}
        return {
            "alphabet": ALPHABET, "blank": "_", "tapes": tape_names,
            "tape_count": len(tape_names), "initial_state": "GET",
            "entry_states": [q for q in ("GET", "SET", "APPEND", "LEN", "COPY") if q in states],
            "halting_states": sorted(states - set(self.active)),
            "active_tape": {q: tape_names.index(tape) for q, tape in sorted(self.active.items())},
            "transitions": [[q, s, r.target, r.write, r.move] for (q, s), r in sorted(self.rules.items())],
            "entries": [dict(state=q, **r.__dict__) for (q, _), r in sorted(self.rules.items())],
            "states": len(self.active), "rules": len(self.rules),
            "default": "missing transition rejects; execution history permits exact rollback",
        }


def build_table():
    t = Table()
    # Initial address copy: A is never modified, and C begins completely blank.
    for op in ("GET", "SET", "APPEND", "LEN"):
        t.add(op, "R", "_", op + ":rstart", "^", 1)
        t.add(op + ":rstart", "C", "_", op + ":cstart", "^", 1)
        if op in ("GET", "SET"):
            t.add(op + ":cstart", "A", "^", op + ":afirst", move=1)
            for state in ("afirst", "acopy"):
                for bit in ("0", "1"):
                    t.add(op + ":" + state, "A", bit, op + ":put" + bit, move=1)
            for bit in ("0", "1"):
                t.add(op + ":put" + bit, "C", "_", op + ":acopy", bit, 1)
            t.add(op + ":acopy", "A", "$", op + ":cend")
            t.add(op + ":cend", "C", "_", op + ":crewind", "$")
            t.rewind(op + ":crewind", "C", op + ":mstart")
            t.add(op + ":mstart", "M", "^", op + ":czero", move=1)
            t.add(op + ":czero", "C", "^", op + ":z0", move=1)
            for flag in (0, 1):
                for bit in ("0", "1"):
                    t.add(op + ":z" + str(flag), "C", bit, op + ":z" + str(int(flag or bit == "1")), move=1)
            t.add(op + ":z0", "C", "$", op + ":zero_back")
            t.add(op + ":z1", "C", "$", op + ":nonzero_back")
            t.rewind(op + ":zero_back", "C", op + ":target_first")
            t.rewind(op + ":nonzero_back", "C", op + ":skip_first")
            for state in ("skip_first", "skip"):
                for bit in ("0", "1"):
                    if op == "GET":
                        t.add(op + ":" + state, "M", bit, op + ":skip", move=1)
                    else:
                        t.add(op + ":" + state, "M", bit, op + ":skipput" + bit, move=1)
            if op == "SET":
                for bit in ("0", "1"):
                    t.add(op + ":skipput" + bit, "R", "_", op + ":skip", bit, 1)
                t.add(op + ":skip", "M", "#", op + ":skiphash", move=1)
                t.add(op + ":skiphash", "R", "_", op + ":decstart", "#", 1)
            else:
                t.add(op + ":skip", "M", "#", op + ":decstart", move=1)
            t.add(op + ":decstart", "C", "^", op + ":borrow", move=1)
            t.add(op + ":borrow", "C", "0", op + ":borrow", "1", 1)
            t.add(op + ":borrow", "C", "1", op + ":decback", "0")
            t.rewind(op + ":decback", "C", op + ":czero")
            # GET copies only the addressed word. SET skips it and inserts V.
            for state in ("target_first", "target"):
                for bit in ("0", "1"):
                    if op == "GET":
                        t.add(op + ":" + state, "M", bit, op + ":targetput" + bit, move=1)
                    else:
                        t.add(op + ":" + state, "M", bit, op + ":target", move=1)
            if op == "GET":
                for bit in ("0", "1"):
                    t.add(op + ":targetput" + bit, "R", "_", op + ":target", bit, 1)
                t.add(op + ":target", "M", "#", "finish")
            else:
                t.add(op + ":target", "M", "#", "SET:vstart", move=1)
        elif op == "APPEND":
            t.add(op + ":cstart", "C", "_", "APPEND:mstart", "$")
            t.add("APPEND:mstart", "M", "^", "APPEND:tail_first", move=1)
        else:
            t.add("LEN:cstart", "C", "_", "LEN:cend", "0", 1)
            t.add("LEN:cend", "C", "_", "LEN:cback", "$")
            t.rewind("LEN:cback", "C", "LEN:mstart")
            t.add("LEN:mstart", "M", "^", "LEN:first", move=1)
            for state in ("first", "word"):
                for bit in ("0", "1"):
                    t.add("LEN:" + state, "M", bit, "LEN:word", move=1)
            t.add("LEN:word", "M", "#", "LEN:incstart", move=1)
            t.add("LEN:incstart", "C", "^", "LEN:carry", move=1)
            t.add("LEN:carry", "C", "1", "LEN:carry", "0", 1)
            t.add("LEN:carry", "C", "0", "LEN:incback", "1")
            t.add("LEN:carry", "C", "$", "LEN:extend", "1", 1)
            t.add("LEN:extend", "C", "_", "LEN:incback", "$")
            t.rewind("LEN:incback", "C", "LEN:first")
            t.add("LEN:first", "M", "$", "LEN:outstart")
            t.add("LEN:outstart", "C", "^", "LEN:out", move=1)
            for bit in ("0", "1"):
                t.add("LEN:out", "C", bit, "LEN:put" + bit, move=1)
                t.add("LEN:put" + bit, "R", "_", "LEN:out", bit, 1)
            t.add("LEN:out", "C", "$", "finish")
    # SET and APPEND insert an arbitrary-length V serially.
    for op in ("SET", "APPEND"):
        t.add(op + ":vstart", "V", "^", op + ":vfirst", move=1)
        for state in ("vfirst", "value"):
            for bit in ("0", "1"):
                t.add(op + ":" + state, "V", bit, op + ":vput" + bit, move=1)
        for bit in ("0", "1"):
            t.add(op + ":vput" + bit, "R", "_", op + ":value", bit, 1)
        t.add(op + ":value", "V", "$", op + ":vhash")
        t.add(op + ":vhash", "R", "_", "SET:tail_first" if op == "SET" else "finish", "#", 1)
        for state in ("tail_first", "tail"):
            for bit in ("0", "1"):
                t.add(op + ":" + state, "M", bit, op + ":tailput" + bit, move=1)
        for bit in ("0", "1"):
            t.add(op + ":tailput" + bit, "R", "_", op + ":tail", bit, 1)
        t.add(op + ":tail", "M", "#", op + ":tailhash", move=1)
        t.add(op + ":tailhash", "R", "_", op + ":tail_first", "#", 1)
        t.add(op + ":tail_first", "M", "$", "finish" if op == "SET" else "APPEND:vstart")
    t.add("finish", "R", "_", "reset:M", "$")
    for i, tape in enumerate(TAPES):
        t.rewind("reset:" + tape, tape, "ACCEPT" if i == len(TAPES) - 1 else "reset:" + TAPES[i + 1])
    return t


TABLE = build_table()
RULE_KEYS = sorted(TABLE.rules)
RULE_ID = {key: i for i, key in enumerate(RULE_KEYS)}


def build_copy_table():
    t = Table()
    for s in ("^", "0", "1", "#", "$"):
        t.add("COPY", "R", s, "copyput:" + s, move=0 if s == "$" else 1)
        t.add("copyput:" + s, "P", "_", "copyback:R" if s == "$" else "COPY", s, 0 if s == "$" else 1)
    t.rewind("copyback:R", "R", "copyback:P")
    t.rewind("copyback:P", "P", "COPYDONE")
    return t


COPY_TABLE = build_copy_table()


@dataclass
class Tape:
    cells: dict = field(default_factory=dict)
    head: int = 0

    def read(self):
        return self.cells.get(self.head, "_")

    def write(self, symbol):
        if symbol == "_":
            self.cells.pop(self.head, None)
        else:
            self.cells[self.head] = symbol

    @classmethod
    def from_text(cls, text):
        return cls(dict(enumerate(text)))


class Machine:
    """Base table plus finite-symbol history tape for injective execution.

    The history head always scans its next blank. A history symbol is one rule
    number from the fixed table, not an address, word, or saved Python object.
    reverse_step uses only the adjacent history symbol and the indicated local
    data cell. This is an explicit reversible simulation, still multitape.
    """
    def __init__(self, op, memory, address="^0$", value="^0$"):
        if op not in ("GET", "SET", "APPEND", "LEN"):
            raise ValueError(op)
        self.state = op
        self.tapes = {"M": Tape.from_text(memory), "A": Tape.from_text(address), "V": Tape.from_text(value), "C": Tape(), "R": Tape()}
        self.history = Tape()
        self.steps = 0

    def snapshot(self):
        return self.state, {k: (dict(t.cells), t.head) for k, t in self.tapes.items()}, dict(self.history.cells), self.history.head

    def step(self):
        if self.state in ("ACCEPT", "REJECT"):
            return False
        tape = self.tapes[TABLE.active[self.state]]
        key = self.state, tape.read()
        if key not in TABLE.rules:
            # Rejection is an external status, not a destructive state change.
            return False
        r = TABLE.rules[key]
        assert self.history.read() == "_"
        self.history.write(RULE_ID[key])
        self.history.head += 1
        tape.write(r.write)
        tape.head += r.move
        self.state = r.target
        self.steps += 1
        return True

    def reverse_step(self):
        if self.history.head == 0:
            return False
        self.history.head -= 1
        rule_id = self.history.read()
        q, old = RULE_KEYS[rule_id]
        r = TABLE.rules[q, old]
        assert self.state == r.target
        tape = self.tapes[r.tape]
        tape.head -= r.move
        assert tape.read() == r.write
        tape.write(old)
        self.state = q
        self.history.write("_")
        return True

    def run(self, max_steps=1_000_000):
        for _ in range(max_steps):
            if not self.step():
                return self.state == "ACCEPT"
        raise RuntimeError("transition budget exceeded")

    def rollback(self):
        while self.reverse_step():
            pass


def bits_frame(bits):
    bits = tuple(bits)
    if not bits or any(bit not in (0, 1) for bit in bits):
        raise ValueError("words are nonempty bit sequences")
    return "^" + "".join(map(str, bits)) + "$"


def memory_frame(words):
    return "^" + "".join(bits_frame(word)[1:-1] + "#" for word in words) + "$"


def clean_run(op, words, address=(0,), value=(0,)):
    """Copy the framed result serially, then reverse to clean all work tapes.

    Inputs are restored; P contains the answer/new array. This is Bennett
    compute-copy-uncompute, not in-place SET. Copy and reverse scheduling is
    given concretely in copy_table.json; the finite-symbol execution history
    supplies rollback. Spatial head-routing remains a separate compiler layer.
    """
    m = Machine(op, memory_frame(words), bits_frame(address), bits_frame(value))
    initial = m.snapshot()
    accepted = m.run()
    preserved = Tape()
    copy_steps = 0
    if accepted:
        r = m.tapes["R"]
        assert r.head == 0 and r.read() == "^"
        copy_tapes, state = {"R": r, "P": preserved}, "COPY"
        while state != "COPYDONE":
            tape = copy_tapes[COPY_TABLE.active[state]]
            rule = COPY_TABLE.rules[state, tape.read()]
            tape.write(rule.write)
            tape.head += rule.move
            state = rule.target
            copy_steps += 1
    m.rollback()
    assert m.snapshot() == initial
    out = ""
    if accepted:
        # Boundary adapter only, not used by the machine to address its input.
        while True:
            symbol = preserved.read()
            out += symbol
            if symbol == "$":
                break
            preserved.head += 1
        preserved.head = 0
    return {"accepted": accepted, "result": out, "forward_steps": m.steps,
            "copy_steps": copy_steps, "history_symbols": m.steps,
            "workspace_clean": m.snapshot() == initial}


if __name__ == "__main__":
    path = Path(__file__).with_name("transition_table.json")
    path.write_text(json.dumps(TABLE.payload(), indent=2) + "\n")
    copy_payload = COPY_TABLE.payload()
    copy_payload["initial_state"] = "COPY"
    path.with_name("copy_table.json").write_text(json.dumps(copy_payload, indent=2) + "\n")
    print(json.dumps({"states": len(TABLE.active), "rules": len(TABLE.rules), "path": str(path)}))
