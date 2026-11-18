"""C200: sector constant at half filling from Theorem H (eta L |sigma_L - r_0| <= C#, C# = 2/pi) inserted in
C151 Theorem 2.1 (formula copied from codex/work_C151/constants_C151.py, K_of).  Mismatch term is zero at half filling.
Floating point evaluation of explicit expressions."""
import json, numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
hH = 2*np.log(1.5)-0.5; alpha = 2*np.pi*hH; conn_rate = 4*np.pi*hH
G = lambda t: -np.log1p(-np.exp(-t))+np.arcsin(np.exp(-t))
IG = np.pi**2/6+np.pi/2*np.log(2)
kap = lambda d: 1-4*np.exp(-2*d)
A_tail = minimize_scalar(lambda d: d+np.pi/(2*kap(d)), bounds=(0.8, 5), method='bounded').fun
ell_B = minimize_scalar(lambda d: d+np.log(4/kap(d)), bounds=(0.8, 5), method='bounded').fun
kap_e = 2.0; amax = np.pi*np.sin(np.pi/4); Csharp = 2/np.pi
def K_of(M0, mismatch):
    c = 1-Csharp/M0
    core = 4*kap_e*ell_B*sum(2**(k+1)*M0*G(alpha*(2**k*M0-Csharp)) for k in range(80))
    tail = 2*np.pi*M0*A_tail
    strad = 2*amax/c
    conn = 2*2*amax*kap_e*IG/(conn_rate*np.pi*c)
    mism = np.pi**2*8.43/8 if mismatch else 0.0
    K = core+tail+strad+conn+mism
    return dict(M0=M0, c=c, K=K, C=5*K+np.pi**2/4, core=core, tail=tail, straddle=strad, connectors=conn, mismatch=mism)
out = {}
for tag, mm in [('half_filling_no_mismatch', False), ('with_C151_mismatch_term', True)]:
    best = min((K_of(M0, mm) for M0 in np.linspace(0.8, 6, 521)), key=lambda r: r['K'])
    out[tag] = best; out[tag+'_M0_2'] = K_of(2.0, mm)
json.dump(out, open('constants_C200.json', 'w'), indent=1, default=float)
for k, v in out.items(): print(k, {a: round(float(b), 3) for a, b in v.items()})

# near half filling: h = ceil(L/2) - n >= 0, eta L sigma_L >= lambda - (h+2)/pi on {lambda >= M0} (PROOF Corollary 4.9)
def K_h(h, M0):
    Cs = (h+2)/np.pi; c = 1-Cs/M0
    if c <= 0.05: return dict(K=np.inf)
    core = 4*kap_e*ell_B*sum(2**(k+1)*M0*G(alpha*(2**k*M0-Cs)) for k in range(80))
    K = core+2*np.pi*M0*A_tail+2*amax/c+2*2*amax*kap_e*IG/(conn_rate*np.pi*c)+np.pi**2*8.43/8
    return dict(h=h, M0=M0, c=c, K=K, C=5*K+np.pi**2/4)
nh = {}
for h in [0, 1, 2, 3, 5, 10, 20]:
    best = min((K_h(h, M0) for M0 in np.linspace(0.8, 40, 1961)), key=lambda r: r['K'])
    nh[str(h)] = best; print('near half filling h =', h, {a: round(float(b), 3) for a, b in best.items()})
out['near_half_filling'] = nh
json.dump(out, open('constants_C200.json', 'w'), indent=1, default=float)
