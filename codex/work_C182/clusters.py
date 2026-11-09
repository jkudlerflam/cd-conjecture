"""C182: cluster ladder E_N^min (N = 1..Nmax) of the Fibonacci-cusp model at one (|w|, phi).
Usage: python3 clusters.py wa phi Nmax R2 R3 R4 R5 [R6] > out.json
For N = 2 the band is scanned over all K (pairband.analyse). For N >= 3 the cluster band is evaluated at the
valley-sum momenta K = n_a k_a + n_b k_b (all splits n_a + n_b = N) and refined by a golden search over |dK| <= 0.15
(N = 3 also gets a coarse full-K scan). Every energy is a Dirichlet upper bound; the box size R is reported, and
for N = 3 the minimum is recomputed at 1.5 R.
"""
import numpy as np, json, sys, time
from pairband import analyse, valleys
import nbody

def golden(f, a, b, it=18):
    g = (np.sqrt(5)-1)/2
    c, d = b-g*(b-a), a+g*(b-a); fc, fd = f(c), f(d)
    for _ in range(it):
        if fc < fd: b, d, fd = d, c, fc; c = b-g*(b-a); fc = f(c)
        else: a, c, fc = c, d, fd; d = a+g*(b-a); fd = f(d)
    return (c, fc) if fc < fd else (d, fd)

def phi_for_delta(wa, delta):
    if delta <= 0: return np.pi/2
    f = lambda ph: (lambda v: (v[1][1]-v[0][1]) if len(v) > 1 else np.inf)(valleys(wa*np.exp(1j*ph)))
    lo, hi = 1.0, np.pi/2
    for _ in range(50):
        mid = 0.5*(lo+hi)
        if f(mid) > delta: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def run(wa, phi, Nmax, Rs):
    w = wa*np.exp(1j*phi)
    v = valleys(w)
    (ka, ea) = v[0]; (kb, eb) = v[1] if len(v) > 1 else v[0]
    out = dict(wa=wa, phi=phi, ka=ka, kb=kb, ea=ea, eb=eb, delta=eb-ea, clusters={})
    t0 = time.time()
    p = analyse(wa, phi, R=Rs[2], nK=181)
    out['pair_mins'] = p['mins']
    out['clusters']['2'] = dict(E=p['mins'][0]['E'], K=p['mins'][0]['K'], R=Rs[2])
    for N in range(3, Nmax+1):
        R = Rs[N]
        cache = {}
        def f(K):
            K = float(np.angle(np.exp(1j*K)))
            key = round(K, 10)
            if key not in cache:
                cache[key] = nbody.lowest(N, K, w, R)[0][0]
            return cache[key]
        cands = sorted(set(round(float(np.angle(np.exp(1j*(na*ka+(N-na)*kb)))), 6) for na in range(N+1)))
        if N == 3:
            cands = sorted(set(cands) | set(np.round(np.linspace(-np.pi, np.pi, 24, endpoint=False), 6)))
        res = []
        for K0 in cands:
            K, E = golden(f, K0-0.15, K0+0.15, it=14 if N <= 4 else 10)
            res.append(dict(K0=K0, K=float(K), E=float(E)))
        res.sort(key=lambda r: r['E'])
        best = res[0]
        rec = dict(E=best['E'], K=best['K'], R=R, local=res[:6])
        if N == 3:
            R2 = int(1.5*R)
            rec['E_1p5R'] = float(nbody.lowest(N, best['K'], w, R2)[0][0]); rec['R2'] = R2
        out['clusters'][str(N)] = rec
        out[f't{N}'] = time.time()-t0
        print(json.dumps(dict(N=N, **rec)), file=sys.stderr, flush=True)
    out['per_particle'] = {N: (out['clusters'][N]['E'])/int(N) - ea for N in out['clusters']}
    return out

if __name__ == "__main__":
    wa = float(sys.argv[1]); a2 = sys.argv[2]
    if a2.startswith('d'):            # d<x>: valley offset delta = x * E_b(phi = pi/2)
        x = float(a2[1:])
        from pairband import Eb_exact
        Eb = Eb_exact(wa)
        phi = phi_for_delta(wa, x*Eb)
    else:
        phi = float(a2)
    Nmax = int(sys.argv[3])
    Rs = {N: int(r) for N, r in zip(range(2, Nmax+1), sys.argv[4:])}
    out = run(wa, phi, Nmax, Rs)
    out['arg'] = a2
    print(json.dumps(out))
