"""C168, Lemma D check (a): ordering energy of an impenetrable two-species gas with unequal masses.
Ring of L sites, N_A = N_B = 2 hard-core particles, NN hopping t_A = 1 + eta, t_B = 1 - eta, no passing
(a hop onto an occupied site is forbidden, so the cyclic order of species is conserved).
Compares the ground energies of the cyclic orders AABB (segregated) and ABAB (alternating).
Prediction (C164 numerics, Lemma D): Delta = E_ABAB - E_AABB > 0, Delta ~ eta^2 / L^2 at low density.
Claude Opus 5.5, 2026-11-08."""
import itertools, sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla


def ground(L, pattern, tA, tB):
    # states: tuple of (position, species) for 4 particles; enumerate occupied sets and assign labels by
    # starting the cyclic pattern at any of its distinct rotations, at the smallest occupied site
    n = len(pattern)
    rots = sorted(set(tuple(pattern[(i + r) % n] for i in range(n)) for r in range(n)))
    states = []
    for occ in itertools.combinations(range(L), n):
        for lab in rots:
            states.append((occ, lab))
    index = {s: i for i, s in enumerate(states)}
    rows, cols, vals = [], [], []
    for i, (occ, lab) in enumerate(states):
        occset = set(occ)
        for p in range(n):
            x = occ[p]; t = tA if lab[p] == 'A' else tB
            for d in (1, -1):
                y = (x + d) % L
                if y in occset:
                    continue
                new = list(zip(occ, lab)); new[p] = (y, lab[p]); new.sort()
                no = tuple(q[0] for q in new); nl = tuple(q[1] for q in new)
                j = index[(no, nl)]
                rows.append(i); cols.append(j); vals.append(-t)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(states), len(states)))
    w = spla.eigsh(H, k=1, which='SA', tol=1e-12)[0]
    return float(w[0])


if __name__ == '__main__':
    for L in (16, 24, 32, 40):
        for eta in (0.0, 0.1, 0.2, 0.3):
            ea = ground(L, 'AABB', 1 + eta, 1 - eta); eb = ground(L, 'ABAB', 1 + eta, 1 - eta)
            # ground energy measured from the band bottom -2(tA+tB) per species pair, i.e. kinetic energy
            d = eb - ea
            print(f'L {L:3d} eta {eta:.2f}  E_AABB {ea:.10f}  E_ABAB {eb:.10f}  Delta {d:+.3e}  Delta*L^2 {d*L*L:+.4f}  Delta*L^2/eta^2 {d*L*L/eta**2 if eta else 0:+.4f}', flush=True)
