"""Check: spectrum of periodic Q-state Potts TL chain  -s*sum e_i  (L sites)
is contained in the union over S^z sectors m and real twists phi of the twisted XXZ ring
with N=2L sites and Delta = s*cos(gamma), sqrt(Q)=2cos(gamma).
Claude Fable 5.2, work_C201."""
import numpy as np, itertools, sys
from math import gcd
from functools import reduce

def potts_H(Q, L, s):
    dim = Q**L
    H = np.zeros((dim, dim))
    E1 = np.ones((Q, Q))/np.sqrt(Q)
    states = list(itertools.product(range(Q), repeat=L))
    idx = {st: i for i, st in enumerate(states)}
    for i, st in enumerate(states):
        # bond terms F = sqrt(Q) delta
        for j in range(L):
            if st[j] == st[(j+1) % L]:
                H[i, i] += np.sqrt(Q)
        # on-site E
        for j in range(L):
            for b in range(Q):
                st2 = list(st); st2[j] = b
                H[idx[tuple(st2)], i] += 1/np.sqrt(Q)
    return -s*H

def xxz_sector(N, Delta, phi, nup):
    """H_X = (1/2) sum [XX+YY+Delta(ZZ-1)], twist e^{i phi} on sigma^+_N sigma^-_1; sector with nup up spins."""
    confs = [c for c in itertools.combinations(range(N), nup)]
    idx = {c: i for i, c in enumerate(confs)}
    H = np.zeros((len(confs), len(confs)), dtype=complex)
    for c in confs:
        i = idx[c]; occ = set(c)
        for b in range(N):
            a, a2 = b, (b+1) % N
            za = 1 if a in occ else -1; zb = 1 if a2 in occ else -1
            H[i, i] += 0.5*Delta*(za*zb - 1)
            # (1/2)(XX+YY) = sigma^+_a sigma^-_a2 + h.c.
            if (a in occ) != (a2 in occ):
                new = set(occ)
                if a in occ:   # sigma^-_a sigma^+_{a2}
                    new.remove(a); new.add(a2); ph = np.exp(-1j*phi) if b == N-1 else 1.0
                else:          # sigma^+_a sigma^-_{a2}
                    new.remove(a2); new.add(a); ph = np.exp(1j*phi) if b == N-1 else 1.0
                H[idx[tuple(sorted(new))], i] += ph
    return H

def lcm(a, b): return a*b//gcd(a, b)

def run(Q, L, s):
    gam = np.arccos(np.sqrt(Q)/2)
    Ep = np.linalg.eigvalsh(potts_H(Q, L, s))
    N = 2*L
    step_den = 12*reduce(lcm, range(1, L+1))
    phis = [2*np.pi*k/step_den for k in range(step_den)]
    worst = 0.0; gs_info = None
    uniq = np.unique(np.round(Ep, 9))
    table = {}
    for m in range(0, L+1):
        nup = L - m
        for phi in phis:
            ev = np.linalg.eigvalsh(xxz_sector(N, s*np.cos(gam), phi, nup))
            table[(m, phi)] = ev
    allmin = min(ev[0] for ev in table.values())
    for E in uniq:
        best = (1e9, None)
        for key, ev in table.items():
            dd = np.min(np.abs(ev - E))
            if dd < best[0]: best = (dd, key)
        worst = max(worst, best[0])
        if abs(E - Ep[0]) < 1e-9: gs_info = best
    print(f"Q={Q} L={L} s={s:+d}: #distinct Potts levels={len(uniq)}, worst mismatch={worst:.2e}, "
          f"E_P={Ep[0]:.10f}, GS at (m,phi/pi)=({gs_info[1][0]},{gs_info[1][1]/np.pi:.4f}) [2gamma/pi={2*gam/np.pi:.4f}], "
          f"min over (m,phi grid) XXZ={allmin:.10f}")

if __name__ == "__main__":
    for Q, Lmax in ((2, 5), (3, 4), (4, 4)):
        for s in (+1, -1):
            for L in range(2, Lmax+1):
                run(Q, L, s)
