"""C164 Task 1: classical lattice screen for two-valley soliton bands (Claude Opus agent, 2026-11-08).

For a bond v (est.py coordinates, normalized to the unit 8-sphere as in est.py):
  1. uniform background: global minimum on S^2, LSW positivity, no periodic texture lower (rings 2..Lmax);
  2. isolated localized textures on an open chain (N free sites, both ends clamped to bg): multi-start relaxation
     from 2 pi loops about random axes with random widths and tilts; keep distinct local minima with a large
     excursion (max |n - bg| > 1, i.e. a "2 pi twist" class texture, not a small magnon lump);
     each carries energy E (relative to bg), solid angle Omega and Haldane momentum P = Omega/2 mod 2 pi;
  3. two-valley flag: two distinct soliton minima A, B with |P_A - P_B| (mod 2 pi) > 0.1, both with E > 0.
Pair interactions (like and unlike valleys, both orders) are computed by pair_interactions() at fixed core
positions (one core site of each soliton clamped to its isolated value).
"""
import sys, json, time
import numpy as np
import clat

TWO_PI = 2 * np.pi


def norm_v(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def pdist(a, b):
    d = abs((a - b) % TWO_PI)
    return min(d, TWO_PI - d)


def find_solitons(J, f, bg, N=40, nax=6, seed=0, widths=(1.5, 4.0), tilts=(-0.9, 0.0, 0.9)):
    rng = np.random.default_rng(seed)
    ch = clat.Chain(J, f, bg, N)
    sols = []
    u, w = clat.frame(bg)
    axes = [np.cos(t) * u + np.sin(t) * w for t in np.linspace(0, np.pi, nax, endpoint=False)]
    for ax in axes:
        for wd in widths:
            for tl in tilts:
                for turns in (1.0, -1.0):
                    x0 = ch.loop_start(ax + 0.05 * rng.normal(size=3), wd, (N + 1) / 2, tilt=tl, turns=turns)
                    E, n = ch.relax(x0, gtol=1e-9, maxiter=4000)
                    exc = np.linalg.norm(n - bg, axis=1).max()
                    if exc < 1.0:
                        continue
                    # localized: the texture should be back to bg well inside the chain
                    edge = max(np.linalg.norm(n[:3] - bg, axis=1).max(), np.linalg.norm(n[-3:] - bg, axis=1).max())
                    Om = ch.omega(n.ravel())
                    P = clat.momentum(Om)
                    dup = False
                    for s in sols:
                        if abs(s['E'] - E) < 1e-4 * max(1.0, abs(E)) and pdist(s['P'], P) < 0.02:
                            dup = True; break
                    if dup:
                        continue
                    # polish and test that it is a genuine local minimum (tangent Hessian >= 0)
                    E, n = ch.relax(n.ravel(), gtol=1e-11, maxiter=30000)
                    H, _ = ch.hessian_tangent(n)
                    hmin = float(np.linalg.eigvalsh(H)[0])
                    Om = ch.omega(n.ravel()); P = clat.momentum(Om)
                    edge = max(np.linalg.norm(n[:3] - bg, axis=1).max(), np.linalg.norm(n[-3:] - bg, axis=1).max())
                    exc = np.linalg.norm(n - bg, axis=1).max()
                    rec = dict(E=float(E), Om=float(Om), P=float(P), exc=float(exc), edge=float(edge), hmin=hmin, n=n)
                    if exc < 1.0 or hmin < -1e-7:
                        continue
                    dup = False
                    for s in sols:
                        if abs(s['E'] - E) < 1e-4 * max(1.0, abs(E)) and pdist(s['P'], P) < 0.02:
                            dup = True; break
                    if not dup:
                        sols.append(rec)
    sols.sort(key=lambda s: s['E'])
    return sols


def two_valley(sols, dP=0.1, edge_max=0.05):
    good = [s for s in sols if s['edge'] < edge_max]
    if len(good) < 2:
        return None
    A = good[0]
    for B in good[1:]:
        if pdist(A['P'], B['P']) > dP:
            return A, B
    return None


def analyze(v, N=40, Lmax=16, seed=0):
    v = norm_v(v)
    J, f = clat.J_f(v)
    st = clat.stable_uniform(J, f, Lmax=Lmax, nstart=3, seed=seed)
    out = dict(v=v.tolist(), stable=st['stable'], lsw_min=st['lsw_min'], om_min=st['om_min'],
               texture_gain=st['texture_gain'], bg=st['n'].tolist())
    if not st['stable']:
        return out
    bg = st['n']
    sols = find_solitons(J, f, bg, N=N, seed=seed)
    out['nsol'] = len(sols)
    out['sols'] = [dict(E=s['E'], P=s['P'], exc=s['exc'], edge=s['edge'], hmin=s['hmin']) for s in sols[:6]]
    tv = two_valley(sols)
    out['two_valley'] = tv is not None
    if tv is not None:
        A, B = tv
        out['A'] = dict(E=A['E'], P=A['P']); out['B'] = dict(E=B['E'], P=B['P'])
        out['dE_AB'] = B['E'] - A['E']
        out['gap_ratio'] = min(A['E'], B['E']) / max(st['om_min'], 1e-9)
        # classical candidate: soliton A is the lowest excitation (below the magnon gap), E_A > 0 (background is
        # the ground state), and the second valley is within 30 percent of it (tunable to degeneracy)
        out['cand'] = bool(A['E'] > 0 and A['E'] < st['om_min'] and (B['E'] - A['E']) < 0.3 * A['E'])
    return out


# ------------------------------------------------------------------------------------------ pair interactions
def core_site(n, bg):
    return int(np.argmax(np.linalg.norm(n - bg, axis=1)))


def pair_energy(J, f, bg, nA, nB, R, margin=14, relax_iter=20000):
    """Two solitons with core sites R apart on an open clamped chain; both core spins clamped.
    Returns U = E_pair - E_A - E_B (E_A, E_B isolated, same clamp), final texture."""
    cA, cB = core_site(nA, bg), core_site(nB, bg)
    hA = [nA[i] for i in range(len(nA))]
    # crop cores to +-margin around the core site
    def crop(n, c):
        lo, hi = max(0, c - margin), min(len(n), c + margin + 1)
        return n[lo:hi], c - lo
    a, ca = crop(nA, cA); b, cb = crop(nB, cB)
    Ntot = 2 * margin + R + 2 * margin + 4
    x = np.tile(bg, (Ntot, 1)).astype(float)
    p1 = margin + 2; p2 = p1 + R
    # paste A then B; in the overlap region combine by choosing the larger deviation
    for src, c, p in ((a, ca, p1), (b, cb, p2)):
        for i in range(len(src)):
            j = p - c + i
            if 0 <= j < Ntot and np.linalg.norm(src[i] - bg) > np.linalg.norm(x[j] - bg):
                x[j] = src[i]
    clampA, clampB = nA[cA].copy(), nB[cB].copy()
    E = clamped_relax(J, f, bg, x, {p1: clampA, p2: clampB}, relax_iter)
    EA = clamped_relax(J, f, bg, place(a, ca, Ntot, p1, bg), {p1: clampA}, relax_iter)
    EB = clamped_relax(J, f, bg, place(b, cb, Ntot, p2, bg), {p2: clampB}, relax_iter)
    return E[0] - EA[0] - EB[0], E[1]


def place(src, c, Ntot, p, bg):
    x = np.tile(bg, (Ntot, 1)).astype(float)
    for i in range(len(src)):
        j = p - c + i
        if 0 <= j < Ntot:
            x[j] = src[i]
    return x


def clamped_relax(J, f, bg, x0, clamps, maxiter):
    """Open chain with bg at both ends and some interior spins clamped. Energy relative to all-bg."""
    from scipy.optimize import minimize
    Ntot = len(x0)
    free = np.array([i for i in range(Ntot) if i not in clamps])
    E_ref = (Ntot + 1) * (bg @ J @ bg) + Ntot * (f @ bg)
    def fun(xf):
        V = xf.reshape(-1, 3); r = np.linalg.norm(V, axis=1); nf = V / r[:, None]
        n = np.empty((Ntot, 3)); n[free] = nf
        for i, c in clamps.items():
            n[i] = c
        full = np.vstack([bg, n, bg])
        E = np.einsum('ia,ab,ib->', full[:-1], J, full[1:]) + np.sum(n @ f) - E_ref
        G = full[2:] @ J.T + full[:-2] @ J + f
        G = G[free]
        G = (G - np.sum(G * nf, 1)[:, None] * nf) / r[:, None]
        return E, G.ravel()
    from scipy.optimize import minimize
    res = minimize(fun, x0[free].ravel(), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, gtol=1e-10, ftol=1e-16))
    V = res.x.reshape(-1, 3); nf = V / np.linalg.norm(V, axis=1)[:, None]
    n = np.empty((Ntot, 3)); n[free] = nf
    for i, c in clamps.items():
        n[i] = c
    return res.fun, n


def pair_interactions(J, f, bg, nA, nB, Rs=(3, 4, 5, 6, 8, 10, 13)):
    out = {}
    for lab, (X, Y) in dict(AA=(nA, nA), BB=(nB, nB), AB=(nA, nB), BA=(nB, nA)).items():
        out[lab] = [float(pair_energy(J, f, bg, X, Y, R)[0]) for R in Rs]
    out['R'] = list(Rs)
    return out


# ------------------------------------------------------------------------------------------ sampling
def sample_near_hF(rng, eps, kind='N'):
    """J = -1 + eps (S0 + A(D)), B = eps b with (S0, D, b) uniform on S^10; kind 'T' sets b = 0 (S^7)."""
    S = rng.normal(size=(3, 3)); S = (S + S.T) / 2; S -= np.trace(S) / 3 * np.eye(3)
    D = rng.normal(size=3); b = rng.normal(size=3) if kind == 'N' else np.zeros(3)
    nrm = np.sqrt(np.sum(S * S) + D @ D + b @ b)
    S, D, b = S / nrm, D / nrm, b / nrm
    J = -np.eye(3) + eps * S
    J[0, 1] += eps * D[2]; J[1, 0] -= eps * D[2]
    J[1, 2] += eps * D[0]; J[2, 1] -= eps * D[0]
    J[2, 0] += eps * D[1]; J[0, 2] -= eps * D[1]
    # back to est.py coordinates in the principal frame of the symmetric part
    Js = (J + J.T) / 2; Ja = (J - J.T) / 2
    w, R = np.linalg.eigh(Js)
    Ja2 = R.T @ Ja @ R
    f2 = R.T @ (eps * b)
    v = np.array([w[0], w[1], w[2], Ja2[1, 2], Ja2[2, 0], Ja2[0, 1], *f2])
    return v


if __name__ == '__main__':
    # usage: solscreen.py seed count eps kind out.jsonl   (eps = 0 means uniform on S^8)
    seed, count, eps, kind, fn = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4], sys.argv[5]
    rng = np.random.default_rng(seed)
    t0 = time.time()
    with open(fn, 'a') as fo:
        for i in range(count):
            if kind in ('dm', 'hF_q', 'hFT_q'):
                import qscreen
                v, meta = qscreen.sample(rng, {'dm': 'dm', 'hF_q': 'hF', 'hFT_q': 'hFT'}[kind])
                eps = meta.get('eps', eps)
            elif eps > 0:
                v = sample_near_hF(rng, eps, kind)
            else:
                v = rng.normal(size=9)
                if kind == 'T':
                    v[6:] = 0
            try:
                o = analyze(v, seed=seed * 100000 + i)
            except Exception as ex:
                o = dict(v=list(map(float, v)), error=repr(ex))
            o.update(seed=seed, idx=i, eps=eps, kind=kind)
            fo.write(json.dumps(o) + '\n'); fo.flush()
    print('done', fn, time.time() - t0)
