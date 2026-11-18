"""Literal nine-letter bond with a biased vacuum clock and L1 readout.

No marker is introduced. The boundary bias is -sum_(y not L) P_y P_W.
All changes are inside the same repeated Hermitian nearest-neighbour bond.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import scipy.linalg as la

LET = ['L0', 'L1', 'R0', 'R1', 'W', 'A', 'B', 'XA', 'XB']
IDX = {x: i for i, x in enumerate(LET)}
HEADS = ['W', 'A', 'B', 'XA', 'XB']
NOT_L = ['R0', 'R1']+HEADS
NOT_R = ['L0', 'L1']+HEADS
CHARGE = np.array([0, 2, 0, 2, 0, -1, 1, -1, 1])


def bond(p):
    """Return two-site part (including bias), on-site energies, bad projector."""
    W = np.zeros((81, 81), complex)

    def entry(out, inc, value):
        r = IDX[out[0]]*9+IDX[out[1]]
        c = IDX[inc[0]]*9+IDX[inc[1]]
        W[r, c] += value

    for a in (0, 1):
        entry((f'L{a}', 'W'), ('W', f'R{a}'), -p.get('q', 1.))
        for s in 'AB':
            z = (-1)**a if s == 'A' else 1
            phase = p.get('phi'+s, 0.)
            entry((f'L{a}', s), (s, f'R{a}'), -p.get('w', 1.)*z)
            entry((f'L{a}', 'X'+s), ('X'+s, f'R{a}'), -z*np.exp(1j*phase))
    for s in 'AB':
        for y in NOT_L:
            entry((y, 'X'+s), (y, s), -1.)
        for y in NOT_R:
            entry(('X'+s, y), (s, y), -p.get('td', 0.)*np.exp(1j*p.get('psi'+s, np.pi if s == 'A' else np.pi/2)))
    entry(('XA', 'B'), ('W', 'W'), p.get('lam', 0.)*np.exp(-1j*np.pi/4))
    entry(('XB', 'A'), ('W', 'W'), p.get('lam', 0.)*np.exp(1j*np.pi/4))
    entry(('B', 'R0'), ('A', 'R1'), p.get('g', 0.))
    W += W.conj().T
    for y in NOT_L:
        entry((y, 'W'), (y, 'W'), -p.get('bias', 1.))
    onsite = np.zeros(9)
    for h in HEADS:
        onsite[IDX[h]] = p.get('nu_head', 0.)
    for s in 'AB':
        onsite[IDX[s]] += p.get('Delta', 3.)
        onsite[IDX['X'+s]] += p.get('Delta', 3.)+p.get('D', 0.)
    onsite[IDX['L1']] += p.get('hL', 0.)
    bad = np.zeros((81, 81))
    for a in (0, 1):
        for b in (0, 1):
            i = IDX[f'L{a}']*9+IDX[f'R{b}']
            bad[i, i] = 1.
    return W, onsite, bad


def clock(n, q):
    H = -q*(np.eye(n+1, k=1)+np.eye(n+1, k=-1))
    H[0, 0] = -1.
    e, U = la.eigh(H)
    v = U[:, 0]
    if v[0] < 0:
        v = -v
    gamma = np.array([sum(v[j+1:]**2) for j in range(n)])
    return e[0], v, gamma


def cell_words(n):
    from itertools import product
    return [(tuple(IDX[x] for x in
                   [f'L{a}' for a in bits[:j]]+[h]+[f'R{a}' for a in bits[j:]]), h, j, bits)
            for h in HEADS for j in range(n+1) for bits in product((0, 1), repeat=n)]
