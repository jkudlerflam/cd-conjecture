"""C168: detection efficiency of the C155 ring screen (rings 12, 16) for c = 2 floating liquids of a qubit chain.

Proxy family (range 3, one qubit per site, outside B_2): Jordan-Wigner free fermions with complex NN and NNN
hopping and a chemical potential,
    H0 = sum_j [ -t1 e^{i a1} s+_j s-_{j+1} - t2 e^{i a2} s+_j (-Z_{j+1}) s-_{j+2} + h.c. ] - (mu/2) sum_j Z_j,
band eps(k) = -2 t1 cos(k - a1) - 2 t2 cos(2k - a2) - mu (k measured with c_k ~ e^{ikj}),
plus a random generic nearest-neighbour two-site perturbation W of norm eps_w (all 15 Pauli strings on (j, j+1)),
which removes the U(1) and every on-site symmetry. Off resonance the number of Fermi seas (one or two) fixes
c = 1 or 2 for small eps_w (C154 section 4.1: every charged operator carries incommensurate momentum).
Output: calib_c2.jsonl (one line per sample: parameters, number of seas, sea widths, C155 estimator at 12/16).
Claude Opus 5.5, 2026-11-08.  Usage: python3 calib_c2.py NSAMP NPROC SEED EPSW
"""
import json, sys, numpy as np
from multiprocessing import Pool
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]], complex); Z = np.diag([1., -1.]).astype(complex)
SP = np.array([[0, 1], [0, 0]], complex); SM = SP.T.copy()
P = [X, Y, Z]


def k3(a, b, c):
    return np.kron(np.kron(a, b), c)


def seas(t1, a1, t2, a2, mu, n=20000):
    k = np.linspace(-np.pi, np.pi, n, endpoint=False)
    e = -2 * t1 * np.cos(k - a1) - 2 * t2 * np.cos(2 * k - a2) - mu
    occ = e < 0
    if occ.all() or not occ.any():
        return 0, []
    # count runs on the circle
    d = np.diff(np.concatenate([occ, occ[:1]]).astype(int))
    starts = np.nonzero(d == 1)[0]; ends = np.nonzero(d == -1)[0]
    nr = len(starts)
    widths = []
    for s in starts:
        e_after = ends[ends > s]
        en = e_after[0] if len(e_after) else ends[0] + n
        widths.append((en - s) * 2 * np.pi / n)
    return nr, widths


def bond(t1, a1, t2, a2, mu, W):
    h = -t1 * np.exp(1j * a1) * k3(SP, SM, I2) - t2 * np.exp(1j * a2) * k3(SP, -Z, SM)
    h = h + h.conj().T
    h = h - (mu / 2) * k3(Z, I2, I2)
    h = h + np.kron(W, I2)
    return h


def sample(args):
    i, seed, epsw = args
    import est_c155 as est
    rng = np.random.default_rng(seed * 100003 + i)
    t1 = 1.0; t2 = rng.uniform(0.5, 1.5); a1 = rng.uniform(0, 2 * np.pi); a2 = rng.uniform(0, 2 * np.pi)
    mu = rng.uniform(-2 * (t1 + t2), 2 * (t1 + t2))
    nr, w = seas(t1, a1, t2, a2, mu)
    c = rng.normal(size=(4, 4)); c[0, 0] = 0
    W = sum(c[a, b] * np.kron(([I2] + P)[a], ([I2] + P)[b]) for a in range(4) for b in range(4))
    W = W * epsw / np.linalg.norm(c)
    h = bond(t1, a1, t2, a2, mu, W)
    try:
        e = est.estimate(h, Ls=(12, 16), r=3, nev=2)
        e = {k: e[k] for k in ('c_S', 'c_fit', 'c_cas', 'c_obj', 'crit', 'ndeg', 'm0')}
    except Exception as ex:
        e = dict(error=str(ex))
    return dict(i=i, t2=t2, a1=a1, a2=a2, mu=mu, nseas=nr, widths=w, epsw=epsw, est=e)


if __name__ == '__main__':
    ns, npr, seed, epsw = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    with Pool(npr) as p, open(f'calib_c2_s{seed}_e{epsw}.jsonl', 'w') as f:
        for r in p.imap_unordered(sample, [(i, seed, epsw) for i in range(ns)]):
            f.write(json.dumps(r) + '\n'); f.flush()
