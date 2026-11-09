"""C186: N_max convergence of pairbind2.py. Starts the ring tuning from the (phi, mu) found at a lower N_max and
refines with one phi-mu iteration in a narrow window, then measures B exactly as pairbind2.run."""
import numpy as np, json, sys, time
from scipy.optimize import brentq
import pairbind2 as pb
from ring_mixed import model

src, v, Nmax = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
d0 = [json.loads(l) for l in open(src) if json.loads(l)['v'] == v and json.loads(l)['rmode'] == 'locked']
d0 = max(d0, key=lambda d: d['Nmax'])
t0 = time.time()
tu = pb.Tuner(d0['L'], Nmax, d0['u'], v, 'locked')
phi, mu = d0['phi'], d0['mu']
for it in range(2):
    f = lambda ph: tu.valley(tu.ia0, ph, mu)[0] - tu.valley(tu.ib0, ph, mu)[0]
    span = 0.006
    while f(phi - span) * f(phi + span) > 0:
        span *= 2
    phi = brentq(f, phi - span, phi + span, xtol=1e-6)
    g = lambda m_: min(tu.valley(tu.ia0, phi, m_)[0], tu.valley(tu.ib0, phi, m_)[0])
    lo, hi = mu - 0.02, mu + 0.02
    while g(lo) < 0: lo -= 0.02
    while g(hi) > 0: hi += 0.02
    mu = brentq(g, lo, hi, xtol=1e-7)
L = d0['L']
(ea, ia), (eb, ib) = tu.valley(tu.ia0, phi, mu), tu.valley(tu.ib0, phi, mu)
E1 = np.array([tu.E1(i, phi, mu) for i in range(L)])
iK = (ia + ib) % L
T = min(E1[j] + E1[(iK - j) % L] for j in range(L))
e, Nx = tu.ring.lowest(model(d0['u'], phi, v, mu), 2 * np.pi * iK / L, k=4)
out = dict(L=L, Nmax=Nmax, u=d0['u'], v=v, rmode='locked', r=np.sqrt(v), phi=phi, mu=mu, qa=2*np.pi*ia/L, qb=2*np.pi*ib/L,
           ea=ea, eb=eb, K=2*np.pi*iK/L, T=float(T), E_low=float(e[0]), N_low=float(Nx[0]), B=float(T - e[0]),
           levels=[float(x) for x in e], Ns=[float(x) for x in Nx], sec=round(time.time() - t0, 1))
print(json.dumps(out), flush=True)
