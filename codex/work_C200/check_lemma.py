"""C200 checks of the half-filling lemma (Astra 6.2 proof, PROOF.md Part 1). Floating point (mpmath for the delicate ones).
(a) inequality (4): sinh^2(a u) + sin^2 eta >= 2 a sinh^2 u, a = 2 - 2 eta/pi, eta in (0, pi/4]
(b) h(x) <= sin(2eta)/(4pi), h unimodal
(c) identity -D_y(m) = int_{z-d}^{z+d} h, by quadrature against the direct formula
(d) |J(y)| <= sin(2 eta)|cos 2 alpha|/(2 pi) by direct evaluation of J for random y, L"""
import numpy as np
from scipy.integrate import quad
rng = np.random.default_rng(7)
# (a)
worst = np.inf
for eta in np.concatenate([np.geomspace(1e-8, np.pi/4, 400), [np.pi/4]]):
    a = 2-2*eta/np.pi
    u = np.concatenate([np.geomspace(1e-9, 6, 4000)])
    # ratio (sinh^2(au)+sin^2 eta)/(2 a sinh^2 u) computed stably
    r = (np.sinh(a*u)**2+np.sin(eta)**2)/(2*a*np.sinh(u)**2)
    worst = min(worst, r.min())
print('(a) min of [sinh^2(au)+sin^2 eta]/[2a sinh^2 u] =', worst, '(must be >= 1)')
# (b)
wb = 0; nonuni = 0
for eta in np.geomspace(1e-6, np.pi/4, 200):
    a = 2-2*eta/np.pi
    u = np.geomspace(1e-10, 30, 20000)
    h = a*np.sin(2*eta)/(2*np.pi)*np.sinh(u)**2/(np.sinh(a*u)**2+np.sin(eta)**2)
    wb = max(wb, h.max()/(np.sin(2*eta)/(4*np.pi)))
    dh = np.diff(h); sg = np.sign(dh[np.abs(dh) > 1e-18*h.max()])
    if (sg[1:] != sg[:-1]).sum() > 1: nonuni += 1
print('(b) max h/(sin 2eta/4pi) =', wb, ' non-unimodal cases:', nonuni)
# (c)
def Theta(w, eta): return np.arctan(np.tanh(w)/np.tan(eta))/np.pi
errc = 0
for _ in range(300):
    eta = np.exp(rng.uniform(np.log(1e-4), np.log(np.pi/4))); b = np.pi/(np.pi-eta); a = 2/b
    alpha = rng.uniform(0.01, np.pi/4); phi = rng.uniform(0.001, np.pi/2-0.001)
    wp = np.log(np.tan(alpha+phi)/np.tan(alpha))/b if phi < np.pi/2-alpha else np.inf
    wm = np.log(np.tan(alpha)/np.tan(alpha-phi))/b if phi < alpha else np.inf
    D = Theta(wp, eta)-Theta(wm, eta)
    z = np.sin(2*alpha)/np.tan(phi); d = np.cos(2*alpha)
    def h(x):
        if x <= 1: return 0.0
        u = np.arctanh(1/x)
        return a*np.sin(2*eta)/(2*np.pi)*np.sinh(u)**2/(np.sinh(a*u)**2+np.sin(eta)**2)
    pts = [p for p in [1.0, 1.0+eta, 1+10*eta] if z-d < p < z+d]
    val = quad(h, z-d, z+d, points=pts or None, limit=400, epsabs=1e-13, epsrel=1e-11)[0]
    errc = max(errc, abs(val+D)/eta)
print('(c) max |int h + D|/eta =', errc)
# (d) direct J(y) = int s(L(N(v)-N(y))) T(y-v) dv, cell by cell with the antiderivative Theta
def J_direct(eta, y, L):
    b = np.pi/(np.pi-eta)
    N = lambda v: np.arctan(np.exp(b*v))/np.pi
    Ny = N(y)
    # breakpoints: L(N(v)-N(y)) = k, k integer, plus y itself
    kmin = int(np.ceil(-L*Ny+1e-12)); kmax = int(np.floor(L*(0.5-Ny)-1e-12))
    ks = np.arange(kmin, kmax+1)
    lev = Ny+ks/L
    vk = np.log(np.tan(np.pi*lev))/b
    # J = sum_cells int T(y-v)(tau - k - 1/2) dv, tau = L(N(v)-N(y)); integrate by parts:
    # int T(y-v) tau dv over R minus sum over cells of (k+1/2)*mass_k, mass_k = Theta(y-v_lo)-Theta(y-v_hi)
    edges = np.concatenate([[-np.inf], vk, [np.inf]])
    kk = np.concatenate([[kmin-1], ks])            # floor(tau) on each cell
    Th = lambda w: np.arctan(np.tanh(w)/np.tan(eta))/np.pi
    mass = Th(y-edges[:-1])-Th(y-edges[1:])
    # int T(y-v) L (N(v)-N(y)) dv by quadrature in w with matched pairs
    f = lambda w: (np.sin(2*eta)/(np.pi*(np.cosh(2*w)-np.cos(2*eta))))*L*(N(y+w)+N(y-w)-2*Ny)
    I1 = quad(f, 0, 60, points=[eta, 10*eta, 1.0], limit=500, epsabs=1e-14, epsrel=1e-12)[0]
    return I1-np.sum((kk+0.5)*mass)
wd = 0
for _ in range(2000):
    eta = np.exp(rng.uniform(np.log(1e-3), np.log(np.pi/4))); b = np.pi/(np.pi-eta)
    L = np.exp(rng.uniform(np.log(2), np.log(3000)))
    alpha = np.exp(rng.uniform(np.log(1e-4), np.log(np.pi/2-1e-4)))
    y = np.log(np.tan(alpha))/b
    Jv = J_direct(eta, y, L)
    bound = np.sin(2*eta)*abs(np.cos(2*alpha))/(2*np.pi)
    wd = max(wd, abs(Jv)/bound if bound > 1e-12 else 0)
print('(d) max |J|/bound =', wd, '(must be <= 1)')
