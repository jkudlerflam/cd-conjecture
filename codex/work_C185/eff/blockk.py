"""C185: the transverse-DM cusp effective model (PROOF 3.3), including its dN = 1 term, by DMRG with 3-site blocks.

Real-site model (units |t2| = 1, phi = pi/2 after the leading hop shift is cancelled):
  H = P sum_j [ w b+_{j+1} b_j + t2 b+_{j+2} b_j + h.c. + V2 n_j n_{j+2} + V3 n_j n_{j+3} - mu n_j ] P
    + kappa P sum_s [ b+_s (n_{s+2} - n_{s-2}) + h.c. ] P,
with V2 = -2 + v, V3 = -2/3, and at leading order v = kappa^2 (kappa = 4 lam/Delta, v = 16 lam^2/Delta^2).
Blocks of three sites, local basis 000,100,010,001,101 (blockade exact inside a block, penalty U across, as in
C173's block_dmrg.py, which this extends). The kappa term creates or removes a particle at s only when s-1 and s+1
are empty (explicit projector), conditioned on n_{s+2} or n_{s-2}. Every term spans at most two blocks.
N is not conserved. Entanglement at block cuts.
"""
import numpy as np, json, sys, time
from dmrg_d import dmrg

LOC = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1)]


def op_matrix(states, nsites, w, t2, V2, V3, mu, U, kappa, across=None, mu_sites=(), left_end=False, right_end=False):
    idx = {st: i for i, st in enumerate(states)}
    H = np.zeros((len(states),) * 2, dtype=complex)

    def keep(lo, hi):
        return True if across is None else (lo < across <= hi)

    for a, st in enumerate(states):
        occ = list(st); d = 0.0
        for i in range(nsites):
            if i in mu_sites:
                d += -mu * occ[i]
            if i + 2 < nsites and keep(i, i + 2): d += V2 * occ[i] * occ[i + 2]
            if i + 3 < nsites and keep(i, i + 3): d += V3 * occ[i] * occ[i + 3]
            if i + 1 < nsites and keep(i, i + 1): d += U * occ[i] * occ[i + 1]
        H[a, a] += d
        for dist, amp in ((1, w), (2, t2)):
            for i in range(nsites - dist):
                j = i + dist
                if not keep(i, j): continue
                for src, dst, am in ((i, j, amp), (j, i, np.conj(amp))):
                    if occ[src] == 1 and occ[dst] == 0:
                        if dist == 2 and occ[i + 1] == 1: continue
                        new = occ.copy(); new[src] = 0; new[dst] = 1
                        t = tuple(new)
                        if t in idx: H[idx[t], a] += am
        if kappa != 0:
            # creation at s conditioned on s+c occupied (c = +2: sign +, c = -2: sign -); support s-1 .. s+2 or s-2 .. s+1
            for s in range(nsites):
                for c, sg in ((2, 1.0), (-2, -1.0)):
                    q = s + c
                    lo, hi = min(s - 1, q), max(s + 1, q)
                    if q < 0 or q >= nsites: continue
                    # neighbours of s that exist in this window; at a window edge the term belongs to another window
                    # a neighbour outside the window is allowed only past a chain end
                    if s - 1 < 0 and not left_end: continue
                    if s + 1 >= nsites and not right_end: continue
                    lo_c, hi_c = max(lo, 0), min(hi, nsites - 1)
                    if not keep(lo_c, hi_c): continue
                    if occ[q] != 1: continue
                    nb_empty = all(occ[x] == 0 for x in (s - 1, s + 1) if 0 <= x < nsites)
                    if not nb_empty: continue
                    if occ[s] == 0:
                        new = occ.copy(); new[s] = 1
                        t = tuple(new)
                        if t in idx:
                            H[idx[t], a] += sg * kappa
                            H[a, idx[t]] += sg * np.conj(kappa)
    return H


def block_model(w, kappa, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, mu=0.0, U=40.0):
    pair_states = [a + b for a in LOC for b in LOC]
    hb = op_matrix(pair_states, 6, w, t2, V2, V3, mu, U, kappa, across=3)
    hs = op_matrix(LOC, 3, w, t2, V2, V3, mu, U, kappa, across=None, mu_sites=(0, 1, 2))
    hl = op_matrix(LOC, 3, w, t2, V2, V3, mu, U, kappa, across=None, mu_sites=(0, 1, 2), left_end=True)
    hr = op_matrix(LOC, 3, w, t2, V2, V3, mu, U, kappa, across=None, mu_sites=(0, 1, 2), right_end=True)
    return hb, hs, hl, hr


if __name__ == "__main__":
    tag = sys.argv[1]; a = sys.argv[2]
    cfgs = json.load(open(a)) if a.endswith('.json') else json.loads(a)
    f = open(f"bk_{tag}.jsonl", "a")
    for cfg in cfgs:
        t0 = time.time()
        w = cfg['wa'] * np.exp(1j * cfg['phi'])
        hb, hs, hl, hr = block_model(w, cfg['kappa'], V2=cfg['V2'], V3=cfg.get('V3', -2 / 3), mu=cfg['mu'], U=cfg.get('U', 40.0))
        Lb = cfg['Lb']
        out = dmrg(hb, Lb, chis=tuple(cfg['chis']), cutoff=1e-11, seed=cfg.get('seed', 0), onsite=[hl] + [hs] * (Lb - 2) + [hr])
        M = out['M']; S = np.array(out['S'])
        nloc = np.array([sum(s) for s in LOC], float)
        R = [None] * (Lb + 1); R[Lb] = np.ones((1, 1), complex)
        for j in range(Lb - 1, -1, -1):
            R[j] = np.einsum('asb,bc,dsc->ad', M[j], R[j + 1], M[j].conj())
        dens = np.array([np.real(np.einsum('asb,s,bc,asc->', M[j], nloc, R[j + 1], M[j].conj())) for j in range(Lb)])
        Lr = 3 * Lb; m = np.arange(1, Lb); l = 3 * m
        x = np.log((2 * Lr / np.pi) * np.sin(np.pi * l / Lr)); sel = (l >= Lr / 6) & (l <= 5 * Lr / 6)
        c_cc = 6 * np.linalg.lstsq(np.vstack([x[sel], np.ones(sel.sum())]).T, S[sel], rcond=None)[0][0]
        rec = dict(cfg, E=out['E'], c_cc=float(c_cc), n_mid=float(dens[Lb // 4: 3 * Lb // 4].mean() / 3),
                   Smid=float(S[len(S) // 2]), S=[float(v) for v in S], hist=out['hist'], secs=time.time() - t0)
        f.write(json.dumps(rec) + "\n"); f.flush()
        print(json.dumps({k: rec[k] for k in rec if k != 'S'}), flush=True)
