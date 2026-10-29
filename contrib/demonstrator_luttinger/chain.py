#!/usr/bin/env python3
"""Open chains for the Luttinger test of the tuned ten-letter demonstrator.

micro : the full nearest-neighbour d = 10 bond of ../routeA_demonstrator/cell.py, restricted to its
        conserved period-3 type sector P B P P B P ... (types are exactly conserved by the bond, so this
        is an exact sector of the microscopic chain).  Sites have dimension 2 (P), 8 (B), 2 (P).  Every
        two-site term is the corresponding block of cell.bond(), unaltered.  M cells = 3M sites.
logical : the second-order logical chain, one qubit per cell,
        lam^2 [Jxy (s+s- + s-s+) + Jzz ZZ],  Jxy, Jzz from the closed forms (C47, collaborator README).
logical_x : the same with the lam-dependent couplings of the exact two-cell des Cloizeaux bond.

Collaborator's tuning: Jp, Jz free with D0 = D1 = 1, kap = 0.6.  We keep Jp Jz = -2 (so the hopping
Jxy = -9/32 is the same as in the original demonstrator) and vary rho = Jp/|Jz|:
    rho = 2   : Jp = 2,       Jz = -1        (original), Delta = 27/182
    rho = 1   : Jp = sqrt 2,  Jz = -sqrt 2,  Delta = 27/364
    rho = 1/2 : Jp = 1,       Jz = -2,       Delta = 27/728
No four-site (C53) term is used.  Charge U(1): micro site charges from cell.charge_site(); logical Z/2.
The sector of total charge Q is selected by a penalty w (Q - Q0)^2 in the MPO (Q is exactly conserved,
so eigenstates in that sector are unchanged; only other sectors are lifted).
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys
import numpy as np
from scipy import linalg

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'routeA_demonstrator'))
import cell as C  # noqa: E402

BASE = dict(D0=1.0, D1=1.0, kap=0.6, mu=0.55, K=2.0)
SP, SM, Z = C.SP, C.SM, C.Z


def couplings(rho):
    """(Jp, Jz) with Jp Jz = -2 and Jp/|Jz| = rho."""
    Jp = np.sqrt(2 * rho)
    return Jp, -2.0 / Jp


def second_order(Jp, Jz, D=1.0, kap=0.6):
    Jxy = Jp * Jz * kap ** 2 / (4 * D * (D ** 2 - kap ** 2))
    Jzz = 3 * Jp ** 2 * kap ** 4 / (32 * D * (D ** 2 - kap ** 2) * (4 * D ** 2 - kap ** 2))
    return Jxy, Jzz, 2 * Jzz / abs(Jxy)


def exact_two_cell(lam, Jp, Jz):
    import effective as E
    He, w, smin = E.exact_eff(lam, BASE['D0'], BASE['D1'], BASE['kap'], Jp, Jz)
    s = E.summary(He)
    return s['Jxy'], s['Jzz'], s


def two_site_terms(h2, d1, d2, tol=1e-13):
    """Operator Schmidt decomposition h2 = sum_k A_k (x) B_k."""
    T = h2.reshape(d1, d2, d1, d2).transpose(0, 2, 1, 3).reshape(d1 * d1, d2 * d2)
    U, S, Vh = linalg.svd(T, full_matrices=False)
    k = int((S > tol * max(S[0], 1e-300)).sum())
    return [(U[:, i].reshape(d1, d1) * S[i], Vh[i].reshape(d2, d2)) for i in range(k)]


class Chain:
    """Open chain: dims[s], onsite[s] (d x d), bonds[s] = list of (A, B) on (s, s+1), charge[s] (diag)."""

    def __init__(self, dims, onsite, bonds, charge, cell_len, meta):
        self.dims, self.onsite, self.bonds, self.charge = dims, onsite, bonds, charge
        self.L = len(dims)
        self.cell_len = cell_len
        self.meta = meta

    def mpo(self, Q0=None, w=0.0, const=0.0):
        """Upper-triangular MPO, channel 0 = start, D-1 = done, channels 1..nb = bond terms, nb+1 = charge."""
        nb = max(len(b) for b in self.bonds) if self.bonds else 0
        pen = Q0 is not None and w != 0.0
        D = 2 + nb + (1 if pen else 0)
        qch = 1 + nb
        Ws = []
        for s in range(self.L):
            d = self.dims[s]
            W = np.zeros((D, D, d, d))
            I = np.eye(d)
            q = np.diag(self.charge[s])
            W[0, 0] = I
            W[D - 1, D - 1] = I
            loc = self.onsite[s].copy()
            if pen:
                loc = loc + w * (q @ q - 2 * Q0 * q)
                if s == 0:
                    loc = loc + (w * Q0 * Q0 + const) * I
                W[0, qch] = q
                W[qch, qch] = I
                W[qch, D - 1] = 2 * w * q
            elif s == 0 and const:
                loc = loc + const * I
            W[0, D - 1] = loc
            if s < self.L - 1:
                for k, (A, B) in enumerate(self.bonds[s]):
                    W[0, 1 + k] = A
            if s > 0:
                for k, (A, B) in enumerate(self.bonds[s - 1]):
                    W[1 + k, D - 1] = B
            Ws.append(W)
        Ws[0] = Ws[0][:1]
        Ws[-1] = Ws[-1][:, -1:]
        return Ws


def micro_blocks(lam, Jp, Jz):
    h = C.bond(D0=BASE['D0'], D1=BASE['D1'], kap=BASE['kap'], lam=lam, Jp=Jp, Jz=Jz, mu=BASE['mu'], K=BASE['K'])
    iP, iB = np.arange(2), np.arange(2, 10)
    H = h.reshape(10, 10, 10, 10)
    blk = lambda a, b: H[np.ix_(a, b, a, b)].reshape(len(a) * len(b), len(a) * len(b))
    return blk(iP, iB), blk(iB, iP), blk(iP, iP), h


def micro(M, lam, Jp, Jz, field=0.0):
    """3M sites P B P ... ; field adds -field * Z_c on every B site (one-site term, not used by default)."""
    hPB, hBP, hPP, _ = micro_blocks(lam, Jp, Jz)
    q = C.charge_site()
    qP, qB = q[:2], q[2:]
    Zc = np.diag([1.0 if ((i >> 1) & 1) == 0 else -1.0 for i in range(8)])   # letter 4 sL + 2 c + sR
    dims, onsite, bonds, charge = [], [], [], []
    tPB, tBP, tPP = two_site_terms(hPB, 2, 8), two_site_terms(hBP, 8, 2), two_site_terms(hPP, 2, 2)
    for j in range(M):
        for t in 'PBP':
            dims.append(2 if t == 'P' else 8)
            charge.append(qP if t == 'P' else qB)
            onsite.append(-field * Zc if t == 'B' else np.zeros((2, 2)))
        bonds += [tPB, tBP]
        if j < M - 1:
            bonds.append(tPP)
    meta = dict(model='micro', M=M, lam=lam, Jp=Jp, Jz=Jz, field=field, **BASE)
    return Chain(dims, onsite, bonds, charge, 3, meta)


def logical(M, Jxy, Jzz, lam=1.0, field=0.0):
    h = lam ** 2 * (Jxy * (np.kron(SP, SM) + np.kron(SM, SP)) + Jzz * np.kron(Z, Z))
    t = two_site_terms(h, 2, 2)
    dims = [2] * M
    meta = dict(model='logical', M=M, lam=lam, Jxy=Jxy, Jzz=Jzz, field=field)
    return Chain(dims, [-field * Z] * M, [t] * (M - 1), [np.array([0.5, -0.5])] * M, 1, meta)


def encode(logical_mps):
    """Embed a logical MPS (tensors (a, c, b)) into the micro chain: Phi+ on (pL, sL) and (sR, pR)."""
    out = []
    r = 1 / np.sqrt(2)
    for A in logical_mps:
        a, _, b = A.shape
        pL = np.zeros((a, 2, 2 * a))
        for x in range(a):
            for p in range(2):
                pL[x, p, 2 * x + p] = r
        B = np.zeros((2 * a, 8, 2 * b))
        for x in range(a):
            for y in range(b):
                for sL in range(2):
                    for c in range(2):
                        for sR in range(2):
                            B[2 * x + sL, 4 * sL + 2 * c + sR, 2 * y + sR] = A[x, c, y] * r
        pR = np.zeros((2 * b, 2, b))
        for y in range(b):
            for p in range(2):
                pR[2 * y + p, p, y] = 1.0
        out += [pL, B, pR]
    return out


if __name__ == '__main__':
    # checks: blocks reproduce the cell model, charge conservation of every term
    for rho in (2.0, 1.0, 0.5):
        Jp, Jz = couplings(rho)
        print('rho', rho, 'Jp', Jp, 'Jz', Jz, 'second order (Jxy, Jzz, Delta)', second_order(Jp, Jz))
    lam, Jp, Jz = 0.1, *couplings(1.0)
    ch = micro(2, lam, Jp, Jz)
    # dense two-cell Hamiltonian from the chain terms vs cell.inter_cell
    L = ch.L
    dims = ch.dims
    Dt = int(np.prod(dims))

    def emb(ops):
        out = np.array([[1.0]])
        for s in range(L):
            out = np.kron(out, ops.get(s, np.eye(dims[s])))
        return out
    H = sum(emb({s: ch.onsite[s]}) for s in range(L))
    for s in range(L - 1):
        for A, B in ch.bonds[s]:
            H = H + emb({s: A, s + 1: B})
    H00, V, code, qc = C.inter_cell(lam=1.0, Jp=Jp, Jz=Jz, D0=1.0, D1=1.0, kap=0.6)
    Href = H00 + lam * V
    # cell ordering pL sL c sR pR == site ordering (pL)(sL c sR)(pR): identical kron order
    w1 = np.linalg.eigvalsh(H); w2 = np.linalg.eigvalsh(Href)
    shift = w1[0] - w2[0]
    print('two-cell spectrum match (up to the constant -mu per B):', np.abs(w1 - w2 - shift).max(), 'shift', shift)
    Qtot = sum(emb({s: np.diag(ch.charge[s])}) for s in range(L))
    print('charge commutator', np.abs(H @ Qtot - Qtot @ H).max())
