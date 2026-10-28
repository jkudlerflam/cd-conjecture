#!/usr/bin/env python3
"""Second-order logical bond of two 13-letter cells, and the exact two-cell check.

H(lam) = H0 (x) 1 + 1 (x) H0 + lam V,  V = sum (coef) opR(cell 1) opL(cell 2).
Second order: H2 = - P V (H0 (x) 1 + 1 (x) H0)^{-1} V P on the 16-dim code, computed in the product
eigenbasis of H0 (exact, no truncation).  Exact: des Cloizeaux bond from the 16 lowest states of the
two-cell Hamiltonian (dimension 104976, sparse Lanczos).
Logical operators: occupation qubits n_a of flavour a in each cell, basis order (n1, n2) per cell.
Claude (Opus 5.5), 2026-10-28.
"""
import numpy as np
from scipy import linalg
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh
import cell2 as C

I2 = np.eye(2)
SP, SM, NO = C.SP, C.SM, C.NOP


def lop(op, a, cell):
    """Logical operator op on flavour a (0,1) of cell (0,1) in the 16-dim two-cell code (cell 0 first)."""
    one = [I2, I2]; one[a] = op
    single = np.kron(one[0], one[1])
    return np.kron(single, np.eye(4)) if cell == 0 else np.kron(np.eye(4), single)


def basis_ops():
    """Named Hermitian operators used to read off the bond."""
    ops = {}
    for a in (0, 1):
        hop = lop(SP, a, 0) @ lop(SM, a, 1)
        ops[f'hop{a + 1}'] = hop + hop.conj().T                     # c+_a,1 c_a,2 + h.c.
        ops[f'hopi{a + 1}'] = 1j * (hop - hop.conj().T)             # imaginary part
        ops[f'n{a + 1}_1'] = lop(NO, a, 0)
        ops[f'n{a + 1}_2'] = lop(NO, a, 1)
        ops[f'nn{a + 1}'] = lop(NO, a, 0) @ lop(NO, a, 1)
    ops['n1n2_12'] = lop(NO, 0, 0) @ lop(NO, 1, 1)
    ops['n2n1_12'] = lop(NO, 1, 0) @ lop(NO, 0, 1)
    ops['n1n2_11'] = lop(NO, 0, 0) @ lop(NO, 1, 0)
    ops['n1n2_22'] = lop(NO, 0, 1) @ lop(NO, 1, 1)
    return ops


def fit(Heff):
    """Least-squares decomposition of a 16 x 16 Hermitian bond into the named operators plus identity."""
    ops = basis_ops()
    names = ['I'] + list(ops)
    A = np.array([np.eye(16).ravel()] + [ops[k].ravel() for k in ops]).T
    coef, *_ = np.linalg.lstsq(A, Heff.ravel(), rcond=None)
    resid = np.linalg.norm(A @ coef - Heff.ravel())
    return {k: float(np.real(v)) for k, v in zip(names, coef)}, float(resid)


def two_cell(par):
    H0, ports, code, qc = C.cell_ops(par['DN'], par['DC'], par['kap'], par['sg'], par['DCd'], 0.0)
    terms = C.bond_terms(par['Jp'], par['Jz'])
    return H0, ports, code, qc, terms


def second_order(par):
    H0, ports, code, qc, terms = two_cell(par)
    e, U = linalg.eigh(H0)
    G = np.kron(code, code).reshape(324, 324, 16)
    src = np.zeros((324, 324, 16), complex)
    for cf, r, l in terms:
        A, B = ports['R'][r], ports['L'][l]
        src += cf * np.einsum('ia,jb,abk->ijk', A, B, G, optimize=True)
    first = np.einsum('ijk,ijl->kl', G.conj(), src, optimize=True)
    s = np.einsum('ai,bj,abk->ijk', U.conj(), U.conj(), src, optimize=True)          # product eigenbasis
    den = e[:, None] + e[None, :]
    mask = den > 1e-9
    s2 = np.where(mask[:, :, None], s / np.where(mask, den, 1)[:, :, None], 0)
    leak = float(np.abs(np.where(mask[:, :, None], 0, s)).max())       # weight of the source on the kernel
    H2 = -np.einsum('ijk,ijl->kl', s.conj(), s2, optimize=True)
    return first, H2, leak, float(e[4])


def sparse_two_cell(par, lam):
    H0, ports, code, qc, terms = two_cell(par)
    H0s = sps.csr_matrix(np.where(np.abs(H0) > 1e-14, H0, 0))
    I = sps.identity(324, format='csr')
    H = sps.kron(H0s, I) + sps.kron(I, H0s)
    for cf, r, l in terms:
        A = sps.csr_matrix(np.where(np.abs(ports['R'][r]) > 1e-14, ports['R'][r], 0))
        B = sps.csr_matrix(np.where(np.abs(ports['L'][l]) > 1e-14, ports['L'][l], 0))
        H = H + lam * cf * sps.kron(A, B)
    return H.tocsr(), np.kron(code, code)


def exact_eff(par, lam, nev=20):
    H, G = sparse_two_cell(par, lam)
    w, V = eigsh(H, k=nev, which='SA', tol=1e-12)
    o = np.argsort(w); w, V = w[o], V[:, o]
    Phi = G.conj().T @ V[:, :16]
    S = Phi.conj().T @ Phi
    ev, Ev = linalg.eigh(S)
    Sm = Ev @ np.diag(ev ** -0.5) @ Ev.conj().T
    O = Phi @ Sm
    Heff = O @ np.diag(w[:16]) @ O.conj().T
    return Heff, w, float(np.sqrt(ev.min()))


DEFAULT = dict(DN=1.0, DC=1.0, kap=(0.6, 0.6), sg=-1.0, DCd=(0.0, -0.2), Jp=1.0, Jz=1.0)

if __name__ == '__main__':
    import sys, json
    par = dict(DEFAULT)
    first, H2, leak, gap = second_order(par)
    print('cell gap', gap, '|first order|', np.abs(first).max(), 'source weight on kernel', leak)
    c2, r2 = fit(H2)
    print('H2 fit residual', r2)
    print({k: round(v, 6) for k, v in c2.items()})
