"""C168, Lemma D check (b): pseudospin anisotropy at a dilute two-species onset.
Ring of L sites, 2 A and 2 B hard-core particles, hoppings t_A = 1 + eta, t_B = 1 - eta, plus a swap term
-s (|A_j B_{j+1}> <B_j A_{j+1}| + h.c.) that lets unlike neighbours pass (the only exchange channel).
All cyclic orders mix. Output: ground energy and the ground-state weight w_alt on the alternating orders ABAB.
Prediction of Lemma D: J_perp ~ s n^3 while the mass-imbalance J_z ~ eta^2 n^2, so at fixed (eta, s) the
ground state becomes segregated (w_alt -> 0) as n = 4/L -> 0, while at eta = 0 the weight stays finite.
Claude Opus 5.5, 2026-11-08."""
import itertools, sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla


def build(L, tA, tB, s):
    states = []
    for occ in itertools.combinations(range(L), 4):
        for lab in set(itertools.permutations('AABB')):
            states.append((occ, lab))
    index = {st: i for i, st in enumerate(states)}
    rows, cols, vals = [], [], []
    for i, (occ, lab) in enumerate(states):
        pos = dict(zip(occ, lab))
        for x, a in pos.items():
            t = tA if a == 'A' else tB
            for d in (1, -1):
                y = (x + d) % L
                if y in pos:
                    if d == 1 and pos[y] != a:      # swap with the unlike right neighbour (h.c. from the other side)
                        new = dict(pos); new[x], new[y] = pos[y], pos[x]
                        o = tuple(sorted(new)); j = index[(o, tuple(new[q] for q in o))]
                        rows += [i, j]; cols += [j, i]; vals += [-s / 2, -s / 2]
                    continue
                new = dict(pos); del new[x]; new[y] = a
                o = tuple(sorted(new)); j = index[(o, tuple(new[q] for q in o))]
                rows.append(i); cols.append(j); vals.append(-t)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(states), len(states)))
    alt = np.array([1.0 if lab in (('A', 'B', 'A', 'B'), ('B', 'A', 'B', 'A')) else 0.0 for _, lab in states])
    return H, alt


if __name__ == '__main__':
    for eta in (0.0, 0.2):
        for s in (0.5,):
            for L in (12, 16, 20, 24, 28):
                H, alt = build(L, 1 + eta, 1 - eta, s)
                w, v = spla.eigsh(H, k=1, which='SA', tol=1e-12)
                p = np.abs(v[:, 0]) ** 2
                print(f'eta {eta:.2f} s {s:.2f} L {L:3d} E0 {w[0]:.8f} w_alt {np.sum(p*alt):.4f}  (random-order value 1/3)', flush=True)
