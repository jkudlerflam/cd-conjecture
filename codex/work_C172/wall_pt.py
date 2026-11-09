#!/usr/bin/env python3
"""C172 (Claude Opus 5.5 agent, 2026-11-08): first-order (Ising-limit) wall problem, PROOF section 4.3.

For a bond with J_zz = Dz >> all other entries, the Neel states along z are the unperturbed ground states and
the n-wall configurations are degenerate (energy 2 Dz per wall). First-order degenerate perturbation theory is
the projection P_n (H - Dz sum Z Z) P_n onto n-wall configurations. This script builds it exactly:
  n = 1 on odd rings (2L configurations), n = 2 on even rings (L(L-1) configurations),
and compares the bottom of the two-wall spectrum with twice the bottom of the one-wall band.
Binding (a kink-even state below the two-kink continuum) shows as  E2_min(L) - 2 E1_min  -> negative constant.

usage: python3 wall_pt.py 'g=0.5,G=0.1,...' [L2 [L1]]   (fam keywords of C165; Dz is irrelevant here)
"""
import sys
import numpy as np
from ed_kink import fam, bond_matrix


def configs_two(L):
    neel = sum(1 << j for j in range(1, L, 2))  # bit 1 = down; site j down for odd j
    out = {}
    for i in range(L):
        for r in range(1, L):
            m = 0
            for t in range(r):
                m |= 1 << ((i + t) % L)
            s = neel ^ m
            out.setdefault(s, len(out))
    return out


def configs_one(L):
    out = {}
    for b in range(L):          # ferro bond (b, b+1); spins alternate going right from b+1 around to b
        for up in (0, 1):
            s = 0
            spin = up
            j = (b + 1) % L
            for t in range(L):
                if spin:
                    s |= 1 << j
                j = (j + 1) % L
                spin ^= 1 if t < L - 1 else 0
            # now site b has same spin as site b+1 (one defect bond)
            out.setdefault(s, len(out))
    return out


def nwalls(s, L):
    return sum(1 for j in range(L) if ((s >> j) & 1) == ((s >> ((j + 1) % L)) & 1))


def projected(L, J, B, conf):
    Jp = J.copy(); Jp[2, 2] = 0.0          # drop the Ising part (constant in the sector)
    h = bond_matrix(Jp, B)
    n = len(conf)
    H = np.zeros((n, n), complex)
    for s, a in conf.items():
        for j in range(L):
            k = (j + 1) % L
            bj = (s >> j) & 1; bk = (s >> k) & 1
            i_in = 2 * bj + bk
            for o in range(4):
                amp = h[o, i_in]
                if abs(amp) < 1e-15:
                    continue
                s2 = s & ~((1 << j) | (1 << k)) | (((o >> 1) & 1) << j) | ((o & 1) << k)
                b = conf.get(s2)
                if b is not None:
                    H[b, a] += amp
    return H


if __name__ == '__main__':
    kw = eval('dict(' + sys.argv[1] + ')')
    L2 = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    L1 = int(sys.argv[3]) if len(sys.argv) > 3 else L2 + 1
    J, B = fam(**kw)
    c1 = configs_one(L1)
    assert all(nwalls(s, L1) == 1 for s in c1) and len(c1) == 2 * L1
    c2 = configs_two(L2)
    assert all(nwalls(s, L2) == 2 for s in c2)
    H1 = projected(L1, J, B, c1); H2 = projected(L2, J, B, c2)
    assert np.allclose(H1, H1.conj().T) and np.allclose(H2, H2.conj().T)
    e1 = np.linalg.eigvalsh(H1); e2 = np.linalg.eigvalsh(H2)
    print(f'{sys.argv[1]}  L1={L1} E1min={e1[0]:.6f}  L2={L2} E2min={e2[0]:.6f}  b=E2-2E1={e2[0]-2*e1[0]:+.6f}  '
          f'E2 lowest4={np.round(e2[:4],5)}')
