"""Full-spectrum block ED for the coarse seam-pinching functionals of C121.

Models (bond h on sites i,i+1, ring of n sites):
  'xxz' : h = Sx Sx + Sy Sy + Delta Sz Sz + (hz/2)(Sz_i + Sz_{i+1}),  d = 2
  'su3' : h = P (swap),                                             d = 3
Pinching is in the computational (Sz or color) basis, which is conserved
blockwise because the projectors commute with every U(1) charge used below.

Every eigenvalue in every block is kept. Symmetry multiplicities (spin flip
at hz = 0, color permutations for su3) are used only for quantities that are
invariant sums over records (log Z, nonselective channel traces), never for
history-resolved weights.

Conventions follow codex/work_C121 and codex/work_C123:
  F = exp(-tau H/4), E = exp(-tau H/2), tau = beta/N,
  W_s = F P_{s_N} E ... E P_{s_1} F,  q_s = Tr W_s W_s^dagger,
  Q   = sum_s q_s  (nonselective channel trace),
  S   = sum_s q_s^2,  Q_D = channel trace on the 2n ring with sites 0, n pinched.
All traces carry the factor exp(beta E0) of their own ring.
"""
from itertools import product, permutations
from collections import Counter
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp


def _digits(n, d):
    powers = d ** np.arange(n, dtype=np.int64)
    basis = np.arange(d ** n, dtype=np.int64)
    return powers, (basis[:, None] // powers[None, :]) % d


def blocks(n, model, delta=1.0, hz=0.0, reduce=True):
    d = 3 if model == 'su3' else 2
    powers, digits = _digits(n, d)
    counts = np.stack([(digits == c).sum(axis=1) for c in range(d)], axis=1)
    groups = {}
    for s, c in enumerate(counts):
        groups.setdefault(tuple(int(x) for x in c), []).append(s)
    out = []
    for key, states in groups.items():
        mult = 1
        if reduce and model == 'xxz' and hz == 0.0:
            if key[0] < key[1]:
                continue
            mult = 1 if key[0] == key[1] else 2
        if reduce and model == 'su3':
            if tuple(sorted(key, reverse=True)) != key:
                continue
            mult = 6
            for v in Counter(key).values():
                mult //= [1, 1, 2, 6][v]
        states = np.array(states)
        lookup = np.full(d ** n, -1, dtype=np.int64)
        lookup[states] = np.arange(len(states))
        dd = digits[states]
        m = len(states)
        h = np.zeros((m, m))
        col = np.arange(m)
        for i in range(n):
            j = (i + 1) % n
            x, y = dd[:, i], dd[:, j]
            changed = states + (y - x) * powers[i] + (x - y) * powers[j]
            if model == 'su3':
                h[lookup[changed], col] += 1.0
            else:
                sx, sy = 0.5 - x, 0.5 - y          # digit 0 = spin up
                h[col, col] += delta * sx * sy + 0.5 * hz * (sx + sy)
                keep = x != y
                h[lookup[changed[keep]], col[keep]] += 0.5
        ev, vec = eigh(h, check_finite=False, driver='evd')
        out.append(dict(ev=ev, vec=vec, mult=mult, digits=dd, key=key))
    total = sum(len(b['ev']) * b['mult'] for b in out)
    assert total == d ** n, (total, d ** n)
    E0 = min(float(b['ev'][0]) for b in out)
    return dict(blocks=out, E0=E0, n=n, d=d)


def logW(spec, beta):
    """log Tr exp(-beta (H - E0))."""
    return float(logsumexp(np.concatenate(
        [-beta * (b['ev'] - spec['E0']) + np.log(b['mult']) for b in spec['blocks']])))


def marginal(spec, beta, site=0):
    """Thermal one-site marginal in the pinching basis (needs reduce=False
    unless the model is symmetric, in which case it is uniform)."""
    d = spec['d']
    p = np.zeros(d)
    for b in spec['blocks']:
        w = np.exp(-beta * (b['ev'] - spec['E0']))
        diag = (b['vec'] ** 2) @ w
        for s in range(d):
            p[s] += b['mult'] * diag[b['digits'][:, site] == s].sum()
    return p / p.sum()


def channel(spec, beta, N, sites):
    """Nonselective pinched trace Q (times exp(beta E0)), pinching all listed
    sites at each of the N insertion times."""
    total = 0.0
    d = spec['d']
    for b in spec['blocks']:
        ev, vec = b['ev'], b['vec']
        e = (vec * np.exp(-beta * (ev - spec['E0']) / (2 * N))) @ vec.T
        lab = sum(d ** j * b['digits'][:, s] for j, s in enumerate(sites))
        mask = lab[:, None] == lab[None, :]
        rho = e.copy()
        for j in range(N):
            rho *= mask
            if j < N - 1:
                rho = e @ rho @ e
        total += b['mult'] * float(np.sum(rho * e.T))
    return total


def histories(spec, beta, N, site=0):
    """All d^N record weights q_s (times exp(beta E0)). Requires reduce=False."""
    d = spec['d']
    labels = list(product(range(d), repeat=N))
    q = np.zeros(len(labels))
    for b in spec['blocks']:
        assert b['mult'] == 1
        ev, vec = b['ev'], b['vec']
        f = (vec * np.exp(-beta * (ev - spec['E0']) / (4 * N))) @ vec.T
        e = f @ f
        col = b['digits'][:, site]
        for k, hist in enumerate(labels):
            w = f.copy()
            for j, lab in enumerate(hist):
                w = (col == lab)[:, None] * w
                w = (f if j == N - 1 else e) @ w
            q[k] += float(np.sum(w * w))
    return q, labels


def decomposition(small, big, beta, N, small_full=None):
    """Return R, G and the exact pieces of G = R - 2 eps_L + eps_D - H2(pi).

    small_full: unreduced spec of the small ring (for histories/marginal).
    """
    sf = small_full if small_full is not None else small
    L, d = small['n'], small['d']
    defect = big['E0'] - 2 * small['E0']
    lwL, lwD = logW(small, beta), logW(big, beta)
    q, _ = histories(sf, beta, N)
    QL = q.sum()
    QL_chan = channel(small, beta, N, [0])
    QD = channel(big, beta, N, [0, L])
    R = beta * defect + 2 * lwL - lwD
    G = beta * defect + np.log(q @ q) - np.log(QD)
    epsL = lwL - np.log(QL)
    epsD = lwD - np.log(QD)
    pi = q / QL
    H2 = -np.log(pi @ pi)
    ident = R - 2 * epsL + epsD - H2
    p = marginal(sf, beta)
    s2, s4 = p @ p, (p ** 2) @ (p ** 2)
    return dict(L=L, d=d, beta=beta, N=N, R=float(R), G=float(G),
                epsL=float(epsL), epsD=float(epsD), H2=float(H2),
                identity_residual=float(G - ident),
                chan_vs_hist=float(abs(QL - QL_chan) / QL_chan),
                p=[float(x) for x in p],
                cluster_pred_GmR=float(N * np.log(s4 / s2 ** 2)),
                cluster_pred_epsL=float(-N * np.log(s2)),
                GmR=float(G - R), LDelta=float(L * defect))
