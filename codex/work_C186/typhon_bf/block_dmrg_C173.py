"""DMRG of the blockaded second-order model on an open chain, using 3-site blocks (local dimension 5).

Real-site model (units |t2| = 1):
  H = P sum_j [ w b+_{j+1} b_j + t2 b+_{j+2} b_j + h.c. + V2 n_j n_{j+2} + V3 n_j n_{j+3} - mu n_j ] P
On an open chain particles never pass, so fermions = hard-core bosons here. Blockade inside a block is exact
(local basis 000, 100, 010, 001, 101); across a block boundary it is a penalty U n n (U = 40 by default; a 1e4 penalty stalls Lanczos).
Terms of range <= 3 span at most two consecutive blocks.
Entanglement is read at block boundaries, i.e. real cuts l = 3m.
"""
import numpy as np, json, sys, time
from dmrg_d import dmrg

LOC = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1)]
LOC8 = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]


def terms(nsites, w, t2, V2, V3, mu, U, mu_sites):
    """List of operators on nsites real sites as functions on occupation tuples."""
    return dict(w=w, t2=t2, V2=V2, V3=V3, mu=mu, U=U, mu_sites=mu_sites)


def op_matrix(states, nsites, w, t2, V2, V3, mu, U, pairs_only_across=None, mu_sites=(), allow_pass=False):
    """Matrix of the real-site Hamiltonian restricted to the given product states (tuples of occupations).
    pairs_only_across: if set to an integer s, keep only terms touching both sides of the split at s
    (sites < s and >= s); used for the block-bond part."""
    idx = {st: i for i, st in enumerate(states)}
    n = len(states)
    H = np.zeros((n, n), dtype=complex)

    def keep(i, j):
        if pairs_only_across is None:
            return True
        s = pairs_only_across
        return (min(i, j) < s) and (max(i, j) >= s)

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
                        if dist == 2 and occ[i + 1] == 1 and not allow_pass:
                            continue
                        new = occ.copy(); new[src] = 0; new[dst] = 1
                        t = tuple(new)
                        if t in idx:
                            H[idx[t], a] += am
    return H


def block_model(w, t2=-1.0, V2=-2.0, V3=-2.0 / 3.0, mu=0.0, U=40.0, control=False):
    """control=True: no blockade, U = 0, hard-core bosons whose NNN hop may pass (w = 0 gives two
    decoupled XX legs, c = 2 exactly; used to calibrate the c = 2 reading)."""
    loc = LOC8 if control else LOC
    if control:
        U = 0.0
    pair_states = [a + b for a in loc for b in loc]
    # bond: terms crossing the block boundary at real site 3
    hb = op_matrix(pair_states, 6, w, t2, V2, V3, mu, U, pairs_only_across=3, allow_pass=control)
    # onsite: terms inside one block, plus -mu n
    hs = op_matrix(loc, 3, w, t2, V2, V3, mu, U, pairs_only_across=None, mu_sites=(0, 1, 2), allow_pass=control)
    return hb, hs


def run(Lb, w, mu, chis=(16, 32, 48, 64, 64), seed=0, **kw):
    hb, hs = block_model(w, mu=mu, **kw)
    out = dmrg(hb, Lb, chis=chis, cutoff=1e-11, seed=seed, onsite=[hs] * Lb)
    M = out['M']
    nloc = np.array([sum(s) for s in (LOC8 if kw.get('control') else LOC)], dtype=float)
    return out, nloc


def density_profile(M, nloc):
    """<n> per block for an MPS that is left-canonical with the centre at the last site."""
    L = len(M)
    out = np.zeros(L)
    R = np.ones((1, 1), dtype=complex)
    Rs = [None] * (L + 1); Rs[L] = R
    for j in range(L - 1, -1, -1):
        A = M[j]
        Rs[j] = np.einsum('asb,bc,dsc->ad', A, Rs[j + 1], A.conj())
    for j in range(L):
        A = M[j]
        out[j] = np.real(np.einsum('asb,s,bc,asc->', A, nloc, Rs[j + 1], A.conj()))
    return out


def cc_open(S, Lreal):
    """S at block cuts m = 1..Lb-1 (real cut l = 3m). Fit S = (c/6) ln chord + a on the middle 2/3."""
    Lb = len(S) + 1
    m = np.arange(1, Lb)
    l = 3 * m
    x = np.log((2 * Lreal / np.pi) * np.sin(np.pi * l / Lreal))
    sel = (l >= Lreal / 6) & (l <= 5 * Lreal / 6)
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    c6 = np.linalg.lstsq(A, S[sel], rcond=None)[0][0]
    return 6 * c6


if __name__ == "__main__":
    tag = sys.argv[1]
    a = sys.argv[2]
    cfgs = json.load(open(a)) if a.endswith('.json') else json.loads(a)
    f = open(f"bd_{tag}.jsonl", "a")
    for cfg in cfgs:
        t0 = time.time()
        w = cfg['wa'] * np.exp(1j * cfg['phi'])
        kw = {k: cfg[k] for k in ('t2', 'V2', 'V3', 'control') if k in cfg}
        out, nloc = run(cfg['Lb'], w, cfg['mu'], chis=tuple(cfg.get('chis', (16, 32, 48, 64, 64))), **kw)
        S = out['S']
        c = cc_open(S, 3 * cfg['Lb'])
        dens = density_profile(out['M'], nloc)
        rec = dict(cfg, E=out['E'], c_cc=c, n_mid=float(dens[len(dens) // 4: 3 * len(dens) // 4].mean() / 3), Smid=float(S[len(S) // 2]), S=[float(x) for x in S],
                   hist=out['hist'], secs=time.time() - t0)
        f.write(json.dumps(rec) + "\n"); f.flush()
        print(json.dumps({k: rec[k] for k in rec if k != 'S'}), flush=True)
