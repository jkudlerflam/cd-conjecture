#!/usr/bin/env python3
"""Exact local parameter oracles and reversible extraction using one H helper.

This supplies whole-space permutations and an explicit nearest-neighbour gate
word.  The stopping rule is a uniform algorithm; its outer head/collision table
has NOT been compiled here.  No logical data symbol controls that algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass


def alphabet(q, m):
    return ([('D', s) for s in range(q)]
            + [('blank',), ('bit', 0), ('bit', 1)]
            + [('H', j, c) for j in range(m) for c in range(2)]
            + [('T', j, c, s) for j in range(m) for c in range(2)
               for s in range(q)])


def species(symbol):
    if symbol[0] == 'D':
        return symbol[1]
    if symbol[0] == 'T':
        return symbol[3]
    return None


@dataclass(frozen=True)
class Site:
    sigma: tuple
    digit: int = 0
    end: int = 0


def local(op, site, p, m):
    """Each gate is defined on every Sigma x digit x end basis state."""
    x, digit, end = site.sigma, site.digit, site.end
    if op in ('PARAM', 'END'):
        if x[0] == 'H':
            _, j, c = x
            value = ((p >> j) & 1) if op == 'PARAM' else int(j == m - 1)
            x = ('H', j, c ^ value)
    elif op in ('COPY_DIGIT', 'COPY_END'):
        if x[0] == 'H':
            if op == 'COPY_DIGIT':
                digit ^= x[2]
            else:
                end ^= x[2]
    elif op in ('SHIFT', 'SHIFT_DAG'):
        if x[0] in ('H', 'T'):
            delta = 1 if op == 'SHIFT' else -1
            x = (x[0], (x[1] + delta) % m, *x[2:])
    else:
        raise ValueError(op)
    return Site(x, digit, end)


def inverse_op(op):
    return {'SHIFT': 'SHIFT_DAG', 'SHIFT_DAG': 'SHIFT'}.get(op, op)


def apply_instruction(tape, instruction, p, m):
    op, at = instruction
    if op == 'SWAP':
        a, b = tape[at], tape[at + 1]
        # Only Sigma moves: trusted parameter tracks remain at physical sites.
        tape[at] = Site(b.sigma, a.digit, a.end)
        tape[at + 1] = Site(a.sigma, b.digit, b.end)
    else:
        tape[at] = local(op, tape[at], p, m)


def extraction_word(n, m):
    """Reference expansion of a compact walk, not an n-sized stored program.

    The physical helper begins at RIGHT.  At visited site j, PARAM writes
    digit j; END exposes whether its index is m-1.  SHIFT and a left SWAP
    advance unless END=1.  A right SWAP walk and inverse shifts restore it.
    """
    if not 1 <= m <= n:
        raise ValueError('Need at least m sites for m parameter bits')
    word = []
    at = n - 1
    for j in range(m):
        for op in ('PARAM', 'COPY_DIGIT', 'PARAM', 'END', 'COPY_END', 'END'):
            word.append((op, at))
        if j + 1 < m:
            word.extend([('SHIFT', at), ('SWAP', at - 1)])
            at -= 1
    while at < n - 1:
        word.extend([('SWAP', at), ('SHIFT_DAG', at + 1)])
        at += 1
    assert len(word) == 10 * m - 4
    return word


def run_word(initial, word, p, m, reverse=False):
    tape = list(initial)
    instructions = ([(inverse_op(op), at) for op, at in reversed(word)]
                    if reverse else word)
    for instruction in instructions:
        apply_instruction(tape, instruction, p, m)
    return tuple(tape)


def extract_until_end(initial, p, m):
    """Uniform stopping algorithm with a separately pinned H_(0,0) helper.

    Its loop never compares the growing index with m in classical control:
    it reads the END oracle's trusted, deterministic output.  The m argument
    is used only by the q-dependent local matrices, just as by SHIFT itself.
    This function is NOT a finite collision-table implementation.
    """
    assert initial[-1].sigma == ('H', 0, 0)
    assert all(a.digit == a.end == 0 for a in initial)
    tape = list(initial)
    at = len(tape) - 1
    trajectory = []
    while True:
        for op in ('PARAM', 'COPY_DIGIT', 'PARAM', 'END', 'COPY_END', 'END'):
            instruction = (op, at)
            apply_instruction(tape, instruction, p, m)
            trajectory.append(instruction)
        if tape[at].end:
            break
        if at == 0:
            raise ValueError('Physical tape shorter than indexed parameter width')
        for instruction in [('SHIFT', at), ('SWAP', at - 1)]:
            apply_instruction(tape, instruction, p, m)
            trajectory.append(instruction)
        at -= 1
    while at < len(tape) - 1:
        for instruction in [('SWAP', at), ('SHIFT_DAG', at + 1)]:
            apply_instruction(tape, instruction, p, m)
            trajectory.append(instruction)
        at += 1
    return tuple(tape), trajectory
