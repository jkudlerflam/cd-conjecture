#!/usr/bin/env python3
"""Full-space local gates of the integrated autonomous controller.

Basis order: head tag, neutral symbol, quantum record. Head tag zero means
no head; tag 1+2*control_index+phase has fresh phase0 or updated phase1.
The quantum record is a neutral qubit.  All columns, including zero/two
heads and malformed neutral symbols, are explicitly defined.
"""

from __future__ import annotations

from integrated_readout import IntegratedReadout


class WorkGates:
    def __init__(self, machine=None):
        self.machine = IntegratedReadout() if machine is None else machine
        self.dimension = 2 * self.machine.local_factor
        a, b = (1 + 1j) / 2, (1 - 1j) / 2
        self.oracles = {
            None: ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
            "CNOT": ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0)),
            "CV": ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, a, b), (0, 0, b, a)),
            "SWAP": ((1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 0), (0, 0, 0, 1)),
        }

    def encode(self, symbol, record, control=None, phase=0):
        if record not in (0, 1) or phase not in (0, 1):
            raise ValueError("Quantum record and head phase are binary")
        if control is None:
            if phase:
                raise ValueError("No-head basis states have canonical phase zero")
            head = 0
        else:
            head = 1 + 2 * self.machine.control_index[control] + phase
        return 2 * (head * self.machine.B + self.machine.symbol_index[symbol]) + record

    def decode(self, index):
        if not 0 <= index < self.dimension:
            raise ValueError("Work basis index out of range")
        rest, record = divmod(index, 2)
        head, neutral = divmod(rest, self.machine.B)
        if head == 0:
            return self.machine.symbols[neutral], record, None, 0
        control, phase = divmod(head - 1, 2)
        return self.machine.symbols[neutral], record, self.machine.controls[control].name, phase

    def f(self, index):
        symbol, record, control, phase = self.decode(index)
        if control is None:
            return index
        source = self.machine.encode(control, symbol)
        table = self.machine.collision if phase else self.machine.inverse
        target_control, target_symbol = self.machine.decode(table[source])
        return self.encode(target_symbol, record, target_control, 1 - phase)

    def seed(self, index):
        symbol, record, control, phase = self.decode(index)
        if symbol.static != "LEFT":
            return index
        if control is None:
            return self.encode(symbol, record, "F:ADVANCE", 0)
        if control == "F:ADVANCE" and phase == 0:
            return self.encode(symbol, record)
        return index

    def initial_allowed(self, index):
        symbol, record, control, phase = self.decode(index)
        return control is None and (symbol.bit, symbol.frontier, symbol.program, symbol.cursor) == (0, 0, 0, 0)

    def g(self, left, right):
        a, ra, qa, pa = self.decode(left)
        b, rb, qb, pb = self.decode(right)
        if (qa is None) == (qb is None):
            return {(left, right): 1+0j}
        on_left = qa is not None
        q, phase = (qa, pa) if on_left else (qb, pb)
        control = self.machine.control_by_name[q]
        forward = phase == 0 and ((control.direction == 1) == on_left)
        reverse = phase == 1 and ((control.direction == -1) == on_left)
        if not (forward or reverse):
            return {(left, right): 1+0j}
        matrix = self.oracles[control.oracle]
        adjoint = control.adjoint ^ reverse
        source_pair = 2 * ra + rb
        result = {}
        for target_pair in range(4):
            coefficient = (complex(matrix[source_pair][target_pair]).conjugate()
                           if adjoint else complex(matrix[target_pair][source_pair]))
            if coefficient == 0:
                continue
            ta, tb = divmod(target_pair, 2)
            if on_left:
                out = self.encode(a, ta), self.encode(b, tb, q, 1 - phase)
            else:
                out = self.encode(a, ta, q, 1 - phase), self.encode(b, tb)
            result[out] = coefficient
        return result


WORK = WorkGates()
WORK_DIM = WORK.dimension
encode = WORK.encode
decode = WORK.decode
initial_allowed = WORK.initial_allowed


def f_column(work):
    return {WORK.f(work): 1+0j}


def seed_column(work):
    return {WORK.seed(work): 1+0j}


def g_column(leftwork, rightwork):
    return WORK.g(leftwork, rightwork)
