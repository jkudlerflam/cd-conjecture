"""C186: binding of two dressed quasiparticles in the mixed cusp model (3.1), on rings with a particle-number cap.

1. Tune the hop phase phi so that the two dressed valleys are degenerate, and mu so that the dressed one-particle
   minimum sits at zero (infinite line, mixed_fewbody.py, N_max_f and box R1; the dressing cloud is compact).
2. Ring of L sites, particle number <= Nmax (ring_mixed.Ring(L, Nmax)). For every grid momentum q the dressed
   one-particle level E1(q) = lowest level with <N> < 1.7. Threshold T(K) = min_q E1(q) + E1(K - q).
3. B(L) = T(K*) - E_low(K*), K* the grid momentum closest to k_a + k_b. Bound pairs give B -> E_b > 0.
"""
import numpy as np, json, sys, time
from ring_mixed import Ring, model
import mixed_fewbody as mf


def tune(u, v, rmode, Nf=4, R1=30):
    r = 0.0 if rmode == "zero" else None
    m = mf.Model(u, np.pi / 2, v, r=r)
    if rmode == "zero":
        phi = np.pi / 2
        mf.tune_mu(m, 1, 4)
    else:
        phi = mf.tune_phi(m, Nf, R1, phi0=np.pi / 2, span=0.2)
    (ka, ea), (kb, eb) = mf.threshold(m, Nf if rmode != "zero" else 1, R1 if rmode != "zero" else 4)
    return phi, m.mu, ka, kb, ea, eb


def run(u, v, rmode, L, Nmax, tuned=None):
    t0 = time.time()
    phi, mu, ka, kb, ea, eb = tuned if tuned is not None else tune(u, v, rmode)
    ring = Ring(L, Nmax)
    m = model(u, phi, v, mu, r=(0.0 if rmode == "zero" else None))
    qs = 2 * np.pi * np.arange(L) / L
    E1 = np.full(L, np.inf); low = {}
    for i, q in enumerate(qs):
        e, Nx = ring.lowest(m, q, k=4)
        low[i] = (e, Nx)
        for ei, ni in zip(e, Nx):
            if 0.5 < ni < 1.7:
                E1[i] = ei
                break
    Kt = np.angle(np.exp(1j * (ka + kb))) % (2 * np.pi)
    iK = int(round(Kt * L / (2 * np.pi))) % L
    T = min(E1[j] + E1[(iK - j) % L] for j in range(L))
    e, Nx = low[iK]
    out = dict(u=u, v=v, rmode=rmode, L=L, Nmax=Nmax, phi=phi, mu=mu, ka=ka, kb=kb, ea=ea, eb=eb,
               K=2 * np.pi * iK / L, T=T, E_low=float(e[0]), N_low=float(Nx[0]), B=float(T - e[0]),
               E1min=float(E1.min()), q1min=float(qs[np.argmin(E1)]),
               levels=[float(x) for x in e], Ns=[float(x) for x in Nx], sec=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    return out


if __name__ == "__main__":
    u = float(sys.argv[1]); rmode = sys.argv[2]
    vs = [float(x) for x in sys.argv[3].split(',')]
    Ls = [int(x) for x in sys.argv[4].split(',')]
    Nmaxs = [int(x) for x in sys.argv[5].split(',')]
    for v in vs:
        tuned = tune(u, v, rmode)
        print(json.dumps(dict(tuned=dict(u=u, v=v, rmode=rmode, phi=tuned[0], mu=tuned[1], ka=tuned[2],
                                         kb=tuned[3], ea=tuned[4], eb=tuned[5]))), flush=True)
        for L in Ls:
            for Nmax in Nmaxs:
                run(u, v, rmode, L, Nmax, tuned)
