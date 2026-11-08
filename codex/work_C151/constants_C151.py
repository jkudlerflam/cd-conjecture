"""C151 constants (PROOF Parts 1-2, 6). Floating point; each analytic majorant is checked on a grid."""
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

out = {}
hH = 2*np.log(1.5)-0.5                 # h(1/2), h(t) = 2 log(1+t) - t
alpha = 2*np.pi*hH                      # Harnack factor in Psi >= alpha eta L sigma
out['h_half'] = hH; out['alpha_H'] = alpha
out['conn_rate'] = 4*np.pi*hH           # Psi_s >= conn_rate * L sigma * s for s <= eta/2
G = lambda t: -np.log1p(-np.exp(-t))+np.arcsin(np.exp(-t))
IG = quad(G, 0, np.inf, limit=200)[0]
out['I_G'] = IG; out['I_G_closed'] = np.pi**2/6+np.pi/2*np.log(2)

# kappa_T(eta/2), 1/cos^2(eta), kappa_eps; grid check of the sup over u of |T(u+is)|/T(u) and |a1'(u+is)|/a1(u)
def kT(eta, s): return np.sin(eta)**2/(np.sin(eta)**2-np.sin(s)**2)
etas = np.linspace(1e-4, np.pi/4, 400)
ke = np.maximum(kT(etas, etas/2), 1/np.cos(etas)**2)
out['kappa_eps_sup_pi4'] = float(ke.max()); out['kappa_eps_eta0'] = float(ke[0])
u = np.linspace(-6, 6, 24001)
worstT = 0; worsta = 0
for eta in [1e-3, 0.05, 0.2, 0.5, np.pi/4]:
    s = eta/2
    T = lambda z: np.sin(2*eta)/(np.pi*(np.cosh(2*z)-np.cos(2*eta)))
    a1 = lambda z: np.sin(eta)/(np.pi*(np.cosh(2*z)+np.cos(eta)))
    da1 = lambda z: -2*np.sin(eta)*np.sinh(2*z)/(np.pi*(np.cosh(2*z)+np.cos(eta))**2)
    rT = np.max(np.abs(T(u+1j*s))/T(u))/kT(eta, s)
    ra = np.max(np.abs(da1(u+1j*s))/a1(u))/(2/np.cos(2*s)**2)
    worstT = max(worstT, rT); worsta = max(worsta, ra)
    # positivity of Re sigma pieces on |Im| < eta
    for ss in np.linspace(0, 0.999*eta, 7):
        assert np.min(T(u+1j*ss).real) > 0 and np.min(a1(u+1j*ss).real) > 0
out['grid_ratio_kT'] = worstT; out['grid_ratio_a1'] = worsta    # both should be <= 1

# tail and band lengths, b >= 1, kappa(d) = 1 - 4 e^{-2d}
kap = lambda d: 1-4*np.exp(-2*d)
At = minimize_scalar(lambda d: d+np.pi/(2*kap(d)), bounds=(0.8, 5), method='bounded')
out['A_tail'] = At.fun; out['A_tail_d'] = At.x
lB = minimize_scalar(lambda d: d+np.log(4/kap(d)), bounds=(0.8, 5), method='bounded')
out['ell_B'] = lB.fun; out['ell_B_d'] = lB.x
Ce = 8.43                                 # C146 Theorem 1

def K_of(M0, mfun, kap_e, amax=np.pi*np.sin(np.pi/4), cmin=0.75):
    core = 4*kap_e*out['ell_B']*sum(2**(k+1)*M0*G(alpha*mfun(2**k*M0)) for k in range(80))
    tail = 2*np.pi*M0*out['A_tail']
    strad = 2*amax/cmin
    conn = 2*2*amax*kap_e*IG/(out['conn_rate']*np.pi*cmin)
    mism = np.pi**2*Ce/8
    K = core+tail+strad+conn+mism
    return dict(K=K, C=5*K+np.pi**2/4, core=core, tail=tail, straddle=strad, connectors=conn, mismatch=mism)

kap_e = out['kappa_eps_sup_pi4']
out['Gc_M0_2'] = K_of(2.0, lambda l: 0.75*l, kap_e)
best = min(((K_of(M0, lambda l: 0.75*l, kap_e)['K'], M0) for M0 in np.linspace(0.5, 4, 351)))
out['Gc_best'] = dict(M0=best[1], **K_of(best[1], lambda l: 0.75*l, kap_e))
# general comparability constant c on the shell (deep core keeps 3/4 >= c from C142 II.6)
for c in [0.75, 0.5, 0.25, 0.1, 0.05]:
    bestc = min(((K_of(M0, lambda l, c=c: c*l, kap_e, cmin=c)['K'], M0) for M0 in np.linspace(0.5, 12, 461)))
    out[f'c_{c}'] = dict(M0=bestc[1], **K_of(bestc[1], lambda l, c=c: c*l, kap_e, cmin=c))
# unconditional: K_0 = {lambda >= th*C_D(eta)} where C142 II.6 gives sigma >= (1-1/th) r; the rest by the sawtooth
def uncond(eta, th):
    CD = (4/np.pi)*(np.log(1/eta)+np.pi/2)
    c = 1-1/th
    return K_of(th*CD, lambda l: c*l, kap_e, cmin=c)['K']
out['uncond_K'] = {f'{e:g}': min(uncond(e, th) for th in np.linspace(1.05, 3, 196)) for e in [1e-1, 1e-2, 1e-4, 1e-6, 1e-12]}
json.dump(out, open('constants_C151.json', 'w'), indent=1, default=float)
for k, v in out.items(): print(k, v)
