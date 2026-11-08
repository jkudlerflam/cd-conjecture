"""C164: sign of the mass-imbalance correction for an impenetrable two-species gas (the term C162 section 5 left
undetermined). Lattice hard-core particles on a ring of Lr sites that cannot pass each other, species hopping
t_A = 1 + eta, t_B = 1 - eta (masses 1/(2 t)), no other interaction. With the cyclic order fixed, the ground-state
energy depends on the species sequence only through the mass imbalance. Compare the segregated order AABB with
the alternating ABAB (and AAABBB / ABABAB / AABABB for six particles).
dE = E(alternating) - E(segregated) > 0 means the imbalance favours segregation (immiscible, ferromagnetic
pseudospin coupling); < 0 favours alternation.
usage: masscheck.py Lr N eta ...
"""
import sys, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def states(Lr, N):
    out = []
    for p1 in range(Lr):
        for rest in itertools.combinations(range(1, Lr), N - 1):
            out.append(tuple([p1] + [(p1 + r) % Lr for r in rest]))
    return out


def energy(Lr, seq, ts):
    N = len(seq)
    S = states(Lr, N)
    idx = {s: i for i, s in enumerate(S)}
    rows, cols, vals = [], [], []
    for i, s in enumerate(S):
        occ = set(s)
        for a in range(N):
            for dpos in (1, -1):
                q = (s[a] + dpos) % Lr
                if q in occ:
                    continue
                t = list(s); t[a] = q; t = tuple(t)
                j = idx.get(t)
                if j is None:
                    continue           # would change the cyclic labelling (cannot happen without passing)
                rows.append(j); cols.append(i); vals.append(-ts[seq[a]])
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(S), len(S)))
    w = spla.eigsh(H, k=1, which='SA', tol=1e-12)[0]
    return float(w[0])


if __name__ == '__main__':
    Lr, N = int(sys.argv[1]), int(sys.argv[2]); etas = [float(x) for x in sys.argv[3:]]
    if N == 4:
        orders = {'AABB': (0, 0, 1, 1), 'ABAB': (0, 1, 0, 1)}
    else:
        orders = {'AAABBB': (0, 0, 0, 1, 1, 1), 'AABABB': (0, 0, 1, 0, 1, 1), 'ABABAB': (0, 1, 0, 1, 0, 1)}
    for eta in etas:
        ts = (1 + eta, 1 - eta)
        E = {k: energy(Lr, s, ts) for k, s in orders.items()}
        ks = list(E)
        print(f'Lr={Lr} N={N} eta={eta}: ' + '  '.join(f'E({k})={E[k]:.8f}' for k in ks) +
              f'   E({ks[-1]}) - E({ks[0]}) = {E[ks[-1]] - E[ks[0]]:+.3e}', flush=True)
