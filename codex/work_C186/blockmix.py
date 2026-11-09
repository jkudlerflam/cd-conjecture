"""C186: DMRG of the mixed second-order cusp model (PROOF eq. 3.1) on open chains, 3-site blocks (local dim 5).

Extends C173's block_dmrg.py (copied logic, C173's dmrg_d.py copied unmodified) by the number-changing vertex
    r sum_j [ b+_j (n_{j+2} - n_{j-2}) + h.c. ],
which has range 3 and therefore spans at most two blocks. Spin language (hard-core bosons). Blockade exact inside a
block, a penalty U n n across block boundaries (as in C173). No particle-number symmetry is used.
"""
import numpy as np, json, sys, time
from dmrg_d import dmrg

LOC = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1)]


def op_matrix(states, nsites, w, t2, V2, V3, mu, U, r, split=None, mu_sites=()):
    idx = {st: i for i, st in enumerate(states)}
    H = np.zeros((len(states), len(states)), dtype=complex)
    keep = (lambda i, j: True) if split is None else (lambda i, j: min(i, j) < split <= max(i, j))
    for a, st in enumerate(states):
        occ = list(st)
        d = 0.0
        for i in range(nsites):
            if i in mu_sites:
                d += -mu * occ[i]
            if i + 2 < nsites and keep(i, i + 2):
                d += V2 * occ[i] * occ[i + 2]
            if i + 3 < nsites and keep(i, i + 3):
                d += V3 * occ[i] * occ[i + 3]
            if i + 1 < nsites and keep(i, i + 1):
                d += U * occ[i] * occ[i + 1]
        H[a, a] += d
        for dist, amp in ((1, w), (2, t2)):
            for i in range(nsites - dist):
                j = i + dist
                if not keep(i, j):
                    continue
                for src, dst, am in ((i, j, amp), (j, i, np.conj(amp))):
                    if occ[src] == 1 and occ[dst] == 0:
                        if dist == 2 and occ[i + 1] == 1:
                            continue
                        new = occ.copy(); new[src] = 0; new[dst] = 1
                        t = tuple(new)
                        if t in idx:
                            H[idx[t], a] += am
        # number-changing vertex: create/annihilate at j conditioned on a particle at j +- 2
        if r != 0:
            for j in range(nsites):
                for part, sgn in ((j + 2, +1.0), (j - 2, -1.0)):
                    if not (0 <= part < nsites) or not keep(j, part) or occ[part] != 1:
                        continue
                    new = occ.copy()
                    if occ[j] == 0:
                        new[j] = 1; am = sgn * r
                    else:
                        new[j] = 0; am = sgn * np.conj(r)
                    t = tuple(new)
                    if t in idx:
                        H[idx[t], a] += am
    return H


def block_model(w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, mu=0.0, U=40.0, r=0.0):
    pair_states = [a + b for a in LOC for b in LOC]
    hb = op_matrix(pair_states, 6, w, t2, V2, V3, mu, U, r, split=3)
    hs = op_matrix(LOC, 3, w, t2, V2, V3, mu, U, r, split=None, mu_sites=(0, 1, 2))
    assert np.abs(hb - hb.conj().T).max() < 1e-12 and np.abs(hs - hs.conj().T).max() < 1e-12
    return hb, hs


def density_profile(M, nloc):
    L = len(M)
    out = np.zeros(L)
    Rs = [None] * (L + 1); Rs[L] = np.ones((1, 1), dtype=complex)
    for j in range(L - 1, -1, -1):
        A = M[j]
        Rs[j] = np.einsum('asb,bc,dsc->ad', A, Rs[j + 1], A.conj())
    for j in range(L):
        A = M[j]
        out[j] = np.real(np.einsum('asb,s,bc,asc->', A, nloc, Rs[j + 1], A.conj()))
    return out


if __name__ == "__main__":
    tag = sys.argv[1]
    cfgs = json.loads(sys.argv[2])
    f = open(f"bm_{tag}.jsonl", "a")
    nloc = np.array([sum(s) for s in LOC], float)
    for cfg in cfgs:
        t0 = time.time()
        w = cfg['u'] * np.exp(1j * cfg['phi'])
        v = cfg['v']; r = cfg.get('r', np.sqrt(v))
        hb, hs = block_model(w, V2=-2.0 + v, mu=cfg['mu'], U=cfg.get('U', 40.0), r=r)
        out = dmrg(hb, cfg['Lb'], chis=tuple(cfg.get('chis', (16, 32, 48, 64, 64))), cutoff=1e-11,
                   seed=cfg.get('seed', 0), onsite=[hs] * cfg['Lb'])
        S = out['S']
        dens = density_profile(out['M'], nloc)
        Lb = cfg['Lb']
        rec = dict(cfg, r=r, E=out['E'], n_mid=float(dens[Lb // 4: 3 * Lb // 4].mean() / 3),
                   Smid=float(np.mean(S[len(S) // 2 - 1: len(S) // 2 + 2])), S=[float(x) for x in S],
                   hist=out['hist'], secs=time.time() - t0)
        f.write(json.dumps(rec) + "\n"); f.flush()
        print(json.dumps({k: rec[k] for k in rec if k not in ('S', 'hist')}), flush=True)
