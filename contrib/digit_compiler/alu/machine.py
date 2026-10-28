"""Fixed-control, asynchronous multi-tape digit arithmetic.

Construction uses Python to emit a finite transition table. Runtime execution
only reads one tape symbol, performs one table lookup, writes one symbol and
moves one head by at most one cell. No arithmetic value is decoded at runtime.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

ALPHABET = ("_", "P", "M", "0", "1", "R")
TAPE_COUNT = 8
OPERATIONS = ("ADD", "SUB", "MUL", "DIV", "MOD", "SHL", "SHR", "AND", "XOR", "EQ", "NE", "LT", "LE", "GT", "GE")
A, B, OUT, X, Y, Z, TMP, COUNT = range(TAPE_COUNT)


class Builder:
    def __init__(self):
        self.active = {}
        self.rules = {}
        self.number = 0

    def state(self, tape, name="s"):
        self.number += 1
        state = f"{name}:{self.number}"
        self.active[state] = tape
        return state

    def rule(self, state, read, nxt, write=None, move=0):
        assert (state, read) not in self.rules and read in ALPHABET
        assert move in (-1, 0, 1)
        self.rules[state, read] = (nxt, read if write is None else write, move)

    def action(self, tape, nxt, *, write=None, move=0, name="action"):
        state = self.state(tape, name)
        for symbol in ALPHABET:
            self.rule(state, symbol, nxt, write, move)
        return state

    def rewind(self, tape, nxt):
        state = self.state(tape, "rewind")
        for symbol in ALPHABET:
            self.rule(state, symbol, nxt if symbol in "PM" else state,
                      move=1 if symbol in "PM" else -1)
        return state

    def end(self, tape, nxt, move=0):
        state = self.state(tape, "seek_end")
        for symbol in ("0", "1", "P", "M"):
            self.rule(state, symbol, state, move=1)
        self.rule(state, "R", nxt, move=move)
        return self.rewind(tape, state)

    def clear(self, tape, nxt):
        # Clear an entire old framed word; return at its first digit.
        home = self.rewind(tape, nxt)
        tail = self.action(tape, home, write="R", move=-1)
        digit = self.action(tape, tail, write="0", move=1)
        back = self.rewind(tape, digit)
        scan = self.state(tape, "clear")
        for bit in "01":
            self.rule(scan, bit, scan, "_", 1)
        self.rule(scan, "R", back, "_", -1)
        return self.rewind(tape, scan)

    def normalize(self, tape, nxt):
        home = self.rewind(tape, nxt)
        end = self.action(tape, home, write="R", move=-1)
        zero_end = self.action(tape, home, write="R", move=-1)
        zero = self.action(tape, zero_end, write="0", move=1)
        scan = self.state(tape, "trim")
        self.rule(scan, "0", scan, "_", -1)
        self.rule(scan, "1", end, move=1)
        self.rule(scan, "P", zero, move=1)
        self.rule(scan, "M", zero, "P", 1)
        # Delete the old endpoint, then trim zeroes. Empty words are allowed here.
        delete_end = self.action(tape, scan, write="_", move=-1)
        return self.end(tape, delete_end)

    def sign(self, tape, negative, nxt):
        set_sign = self.state(tape, "set_sign")
        for symbol in ("P", "M"):
            self.rule(set_sign, symbol, nxt, "M" if negative else "P", 1)
        # Canonical zero must never acquire a negative sign.
        go_header = self.action(tape, set_sign, move=-1)
        return self.zero(tape, self.action(tape, set_sign, move=-1) if not negative else nxt,
                         go_header)

    def zero(self, tape, yes, no):
        scan = self.state(tape, "is_zero")
        self.rule(scan, "0", scan, move=1)
        self.rule(scan, "1", self.rewind(tape, no))
        self.rule(scan, "R", self.rewind(tape, yes))
        return self.rewind(tape, scan)

    def copy(self, src, dst, nxt):
        finish = self.rewind(src, self.rewind(dst, nxt))
        loop = self.state(src, "copy_read")
        for bit in "01":
            write = self.action(dst, loop, write=bit, move=1)
            self.rule(loop, bit, write, move=1)
        self.rule(loop, "R", self.action(dst, finish, write="R"))
        return self.clear(dst, self.rewind(src, loop))

    def increment(self, tape, nxt):
        done = self.rewind(tape, nxt)
        extend = self.action(tape, done, write="R", move=-1)
        loop = self.state(tape, "increment")
        self.rule(loop, "0", done, "1")
        self.rule(loop, "1", loop, "0", 1)
        self.rule(loop, "R", extend, "1", 1)
        return self.rewind(tape, loop)

    def decrement(self, tape, nxt):
        done = self.normalize(tape, nxt)
        loop = self.state(tape, "decrement")
        self.rule(loop, "0", loop, "1", 1)
        self.rule(loop, "1", done, "0")
        self.rule(loop, "R", "ERROR_UNDERFLOW")
        return self.rewind(tape, loop)

    def shift_left(self, tape, nxt):
        done = self.normalize(tape, nxt)
        tail = self.action(tape, done, write="R")
        loop = {b: self.state(tape, "shift_left") for b in "01"}
        for carry in "01":
            for bit in "01":
                self.rule(loop[carry], bit, loop[bit], carry, 1)
            self.rule(loop[carry], "R", tail, carry, 1)
        return self.rewind(tape, loop["0"])

    def shift_right(self, tape, zero, one):
        start = self.state(tape, "shift_right_first")
        for discarded, nxt in (("0", zero), ("1", one)):
            done = self.normalize(tape, nxt)
            put_end = self.action(tape, done, write="R")
            loop = self.state(tape, "shift_right_read")
            for bit in "01":
                advance = self.action(tape, loop, move=1)
                write = self.action(tape, advance, write=bit, move=1)
                self.rule(loop, bit, write, move=-1)
            self.rule(loop, "R", put_end, "_", -1)
            self.rule(start, discarded, loop, move=1)
        return self.rewind(tape, start)

    def compare(self, left, right, less, equal, greater):
        finishes = {r: self.rewind(left, self.rewind(right, out))
                    for r, out in ((-1, less), (0, equal), (1, greater))}
        loops = {r: self.state(left, "compare_left") for r in (-1, 0, 1)}
        for rel in (-1, 0, 1):
            for a in ("0", "1", "R"):
                read_b = self.state(right, "compare_right")
                self.rule(loops[rel], a, read_b, move=0 if a == "R" else 1)
                for b in ("0", "1", "R"):
                    if a == b == "R":
                        self.rule(read_b, b, finishes[rel])
                    else:
                        ai, bi = int(a == "1"), int(b == "1")
                        newrel = -1 if ai < bi else 1 if ai > bi else rel
                        self.rule(read_b, b, loops[newrel], move=0 if b == "R" else 1)
        return self.rewind(left, self.rewind(right, loops[0]))

    def arithmetic(self, left, right, dst, subtract, nxt):
        done = self.normalize(dst, self.rewind(left, self.rewind(right, nxt)))
        put_end = self.action(dst, done, write="R")
        carry_end = self.action(dst, put_end, write="1", move=1)
        loops = {carry: self.state(left, "arith_left") for carry in (0, 1)}
        for carry in (0, 1):
            for a in ("0", "1", "R"):
                read_b = self.state(right, "arith_right")
                self.rule(loops[carry], a, read_b, move=0 if a == "R" else 1)
                for b in ("0", "1", "R"):
                    if a == b == "R":
                        target = "ERROR_UNDERFLOW" if subtract and carry else carry_end if carry else put_end
                        self.rule(read_b, b, target)
                        continue
                    value = int(a == "1") + (-int(b == "1") - carry if subtract else int(b == "1") + carry)
                    bit, newcarry = value % 2, int(value < 0) if subtract else value // 2
                    write = self.action(dst, loops[newcarry], write=str(bit), move=1)
                    self.rule(read_b, b, write, move=0 if b == "R" else 1)
        return self.clear(dst, self.rewind(left, self.rewind(right, loops[0])))

    def bitwise(self, op, neg_a, neg_b, neg_out, nxt):
        finish = self.normalize(OUT, self.rewind(A, self.rewind(B, nxt)))
        put_end = self.action(OUT, finish, write="R")
        loops = {(ca, cb, co): self.state(A, "logic_left")
                 for ca in (0, 1) for cb in (0, 1) for co in (0, 1)}

        def complement(bit, negative, carry):
            if not negative:
                return bit, 0
            return (bit, int(bit == 0)) if carry else (1 - bit, 0)

        for (ca, cb, co), start in loops.items():
            for a in ("0", "1", "R"):
                av, ac = (int(neg_a), 0) if a == "R" else complement(int(a), neg_a, ca)
                rb = self.state(B, "logic_right")
                self.rule(start, a, rb, move=0 if a == "R" else 1)
                for b in ("0", "1", "R"):
                    bv, bc = (int(neg_b), 0) if b == "R" else complement(int(b), neg_b, cb)
                    val = av & bv if op == "AND" else av ^ bv
                    out, oc = complement(val, neg_out, co)
                    target = put_end if a == b == "R" else loops[ac, bc, oc]
                    write = self.action(OUT, target, write=str(out), move=1)
                    self.rule(rb, b, write, move=0 if b == "R" else 1)
        return self.clear(OUT, self.rewind(A, self.rewind(B, loops[int(neg_a), int(neg_b), int(neg_out)])))

    def multiply(self, nxt):
        loop = self.state(Y, "multiply_test_bit")
        again = self.zero(Y, nxt, self.rewind(Y, loop))
        shift_y = self.shift_right(Y, again, again)
        shift_x = self.shift_left(X, shift_y)
        add = self.arithmetic(OUT, X, TMP, False, self.copy(TMP, OUT, shift_x))
        self.rule(loop, "0", shift_x)
        self.rule(loop, "1", add)
        return self.clear(OUT, self.copy(A, X, self.copy(B, Y, again)))

    def divide(self, nxt):
        # Binary long division. Y=quotient, Z=remainder; X is a scanned copy
        # of the numerator whose head runs from its most significant bit left.
        scan = self.state(X, "division_numerator")
        next_digit = self.action(X, scan, move=-1)
        subtract = self.arithmetic(Z, B, TMP, True,
                                   self.copy(TMP, Z, self.increment(Y, next_digit)))
        compare = self.compare(Z, B, next_digit, subtract, subtract)
        shift_q = self.shift_left(Y, compare)
        for bit in "01":
            continuation = self.increment(Z, shift_q) if bit == "1" else shift_q
            self.rule(scan, bit, self.shift_left(Z, continuation))
        for header in "PM":
            self.rule(scan, header, nxt, move=1)
        run = self.copy(A, X, self.clear(Y, self.clear(Z, self.end(X, scan, move=-1))))
        return self.zero(B, "ERROR_DIV_ZERO", run)

    def constant(self, value, nxt):
        assert value in (0, 1)
        put = self.action(OUT, nxt, write=str(value))
        return self.clear(OUT, put)

    def signed_operation(self, op, na, nb):
        halt = "HALT"
        if op in ("ADD", "SUB"):
            nb_effective = nb ^ (op == "SUB")
            if na == nb_effective:
                return self.arithmetic(A, B, OUT, False, self.sign(OUT, na, halt))
            left = self.arithmetic(A, B, OUT, True, self.sign(OUT, na, halt))
            right = self.arithmetic(B, A, OUT, True, self.sign(OUT, nb_effective, halt))
            return self.compare(A, B, right, self.clear(OUT, halt), left)
        if op == "MUL":
            return self.multiply(self.sign(OUT, na ^ nb, halt))
        if op in ("DIV", "MOD"):
            if op == "DIV":
                finish = self.copy(Y, OUT, self.sign(OUT, na ^ nb, halt))
                correction = self.increment(Y, finish) if na != nb else finish
                finish = self.zero(Z, finish, correction)
            else:
                finish = self.copy(Z, OUT, self.sign(OUT, nb, halt))
                correction = self.arithmetic(B, Z, OUT, True, self.sign(OUT, nb, halt)) if na != nb else finish
                finish = self.zero(Z, finish, correction)
            return self.divide(finish)
        if op in ("SHL", "SHR"):
            if nb:
                return "ERROR_NEGATIVE_SHIFT"
            finish = self.sign(OUT, na, halt)
            loop = self.state(COUNT, "shift_loop")
            again = self.zero(COUNT, finish, loop)
            if op == "SHL":
                shift = self.shift_left(OUT, again)
            else:
                shift = self.shift_right(OUT, again, self.increment(OUT, again) if na else again)
            dec = self.decrement(COUNT, shift)
            for symbol in ALPHABET:
                self.rule(loop, symbol, dec)
            return self.copy(A, OUT, self.copy(B, COUNT, again))
        if op in ("AND", "XOR"):
            negative = na and nb if op == "AND" else na != nb
            return self.bitwise(op, na, nb, negative, self.sign(OUT, negative, halt))
        if op in ("EQ", "NE", "LT", "LE", "GT", "GE"):
            def truth(rel):
                result = {"EQ": rel == 0, "NE": rel != 0, "LT": rel < 0,
                          "LE": rel <= 0, "GT": rel > 0, "GE": rel >= 0}[op]
                return self.constant(int(result), halt)
            if na != nb:
                return truth(-1 if na else 1)
            return self.compare(A, B, truth(1 if na else -1), truth(0), truth(-1 if na else 1))
        raise ValueError(op)

    def operation(self, op):
        # The initial read of each sign occurs at cell 0. All other heads begin
        # at cell 1, and every continuation assumes that convention.
        start = self.state(A, op + "_sign_a")
        for na in (False, True):
            bsign = self.state(B, op + "_sign_b")
            self.rule(start, "M" if na else "P", bsign, move=1)
            for nb in (False, True):
                self.rule(bsign, "M" if nb else "P", self.signed_operation(op, na, nb), move=1)
        return start


def build():
    builder = Builder()
    entries = {op: builder.operation(op) for op in OPERATIONS}
    # Total local rule table: an unexpected symbol enters an explicit error
    # state without changing the symbol or head. This is not a proof that every
    # malformed infinite tape halts; the arithmetic theorem assumes framing.
    for state in builder.active:
        for symbol in ALPHABET:
            if (state, symbol) not in builder.rules:
                builder.rule(state, symbol, "ERROR_MALFORMED")
    return {"format": "asynchronous-multitape-tm-v1", "alphabet": list(ALPHABET),
            "blank": "_", "tape_count": TAPE_COUNT, "entries": entries,
            "initial_states": entries,
            "halting_states": ["HALT", "ERROR_UNDERFLOW", "ERROR_DIV_ZERO", "ERROR_NEGATIVE_SHIFT", "ERROR_MALFORMED"],
            "active_tape": builder.active, "transitions": builder.rules}


TABLE = build()


def machine_for(op):
    """Shared-interpreter export; all operations share the same finite table."""
    if op not in OPERATIONS:
        raise ValueError(op)
    return {**TABLE, "initial_state": TABLE["entries"][op]}


def framed(bits_le, negative=False):
    bits = tuple(bits_le)
    if not bits or any(bit not in (0, 1) for bit in bits):
        raise ValueError("A word is a nonempty tuple of bits")
    if len(bits) > 1 and bits[-1] != 1:
        raise ValueError("Words must have no leading magnitude zeroes")
    if bits == (0,) and negative:
        raise ValueError("Negative zero is not canonical")
    return ["M" if negative else "P", *map(str, bits), "R"]


@dataclass
class Configuration:
    state: str
    tapes: list[dict[int, str]]
    heads: list[int]
    steps: int = 0

    def snapshot(self):
        return (self.state, tuple(tuple(sorted(tape.items())) for tape in self.tapes), tuple(self.heads))


class Runner:
    """Reference table interpreter; a host counter measures work, not values."""
    def __init__(self, op, left_bits_le, right_bits_le, left_sign=0, right_sign=0, journal=False):
        values = [framed(left_bits_le, left_sign), framed(right_bits_le, right_sign)]
        values += [framed((0,)) for _ in range(TAPE_COUNT - 2)]
        self.config = Configuration(TABLE["entries"][op],
                                    [dict(enumerate(word)) for word in values],
                                    [0, 0] + [1] * (TAPE_COUNT - 2))
        self.history = [] if journal else None

    def step(self):
        c = self.config
        tape = TABLE["active_tape"][c.state]
        head = c.heads[tape]
        read = c.tapes[tape].get(head, "_")
        event = (c.state, read)
        nxt, write, move = TABLE["transitions"][event]
        if self.history is not None:
            self.history.append(event)
        if write == "_":
            c.tapes[tape].pop(head, None)
        else:
            c.tapes[tape][head] = write
        c.heads[tape] += move
        c.state = nxt
        c.steps += 1

    def undo(self):
        state, read = self.history.pop()
        nxt, write, move = TABLE["transitions"][state, read]
        c = self.config
        assert c.state == nxt
        tape = TABLE["active_tape"][state]
        c.heads[tape] -= move
        head = c.heads[tape]
        assert c.tapes[tape].get(head, "_") == write
        if read == "_":
            c.tapes[tape].pop(head, None)
        else:
            c.tapes[tape][head] = read
        c.state = state

    def run(self, limit=10_000_000):
        while self.config.state not in TABLE["halting_states"]:
            if self.config.steps >= limit:
                raise RuntimeError("Step limit exceeded")
            self.step()
        if self.config.state != "HALT":
            raise ArithmeticError(self.config.state)
        tape = self.config.tapes[OUT]
        bits, position = [], 1
        while tape.get(position) != "R":
            bits.append(int(tape[position]))
            position += 1
        return {"bits_le": tuple(bits), "sign": int(tape[0] == "M"),
                "steps": self.config.steps, "runner": self}


def run_binary(op, left_bits_le, right_bits_le, left_sign=0, right_sign=0, *, journal=False, limit=10_000_000):
    return Runner(op, left_bits_le, right_bits_le, left_sign, right_sign, journal).run(limit)


def write_table(path=None):
    payload = {k: v for k, v in TABLE.items() if k != "transitions"}
    payload["transitions"] = [[state, symbol, *rule] for (state, symbol), rule in TABLE["transitions"].items()]
    target = Path(path) if path else Path(__file__).with_name("transition_table.json")
    target.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    return target


if __name__ == "__main__":
    print(write_table())
