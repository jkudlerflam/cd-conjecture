"""C149 constants. Floating point evaluation of explicit one-variable functions on fine grids.
Every quantity is a closed form in eta; sups are over a geometric grid of 2e5 points in [1e-9, pi/2 - 1e-9]
plus the endpoint limits.  Output constants_C149.json."""
import numpy as np, json
from scipy.optimize import brentq

A_core, A_tail = 2.71765, 3.98873          # C145 / C147 closed forms
def I_of_b(b):                              # C146 Lemma 2.3
    f = lambda e: b*b*e*(4+e)-4*(1-e)**2
    es = brentq(f, 1e-15, 1-1e-15); ss = 0.5*np.log(4/es)
    return (1/b)*(4*ss+(b*b/2)*(5/(1-es)-5+np.log(1-es)))

def pieces(eta, k0):
    g = np.pi-eta; b = np.pi/g; a = np.pi*np.sin(eta); q = 1-2*eta/np.pi; B = np.pi/(2*eta)
    ep = np.pi*eta/g
    TA0 = 1/(np.tan(ep/2)*g); TB0 = 1/(np.tan(eta/2)*np.pi)
    CD = (2*k0+1)*eta/np.tan(eta)/np.pi
    CD1 = (2*k0+1)*2*eta**2*TA0*TB0
    M = 4*CD
    VW = a*(4*I_of_b(b)-8)+12*a*B                     # C146 Theorem 2 second line
    core = (2*a/9)*(1+CD1*A_core/(eta*CD))
    var = VW/6
    tail = 4*a*M*A_tail/eta
    strad = 16*a/3
    outer = 4*a*q*(k0+0.5)
    bdry = 4*a*q*q*(k0+0.5)**2
    K = core+var+tail+strad+outer+bdry
    return dict(CD=CD, CD1=CD1, ratio=CD1/CD, core=core, var=var, tail=tail, strad=strad, outer=outer, bdry=bdry, K=K)

def C1_const(eta):
    """sup_d eta T(d) / int_d^inf T  (transfer lemma)."""
    te = np.tan(eta)
    d = np.concatenate([[0.0], np.geomspace(1e-6*eta, 40, 20000)])
    t = np.tanh(d)
    T = te*(1-t*t)/(np.pi*(t*t+te*te))
    omt = 2/(np.exp(2*d)+1)
    tail = 2*np.arctan(te*omt/(t+te*te))/(2*np.pi)   # = (phi(d) - 2 eta)/(2 pi) = int_d^inf T, stable form
    return float(np.max(eta*T/tail))

out = {}
etas = np.geomspace(1e-9, np.pi/2-1e-9, 20000)
for k0 in [0.3, 0.35, 0.4, 0.45]:
    rows = [pieces(e, k0) for e in etas]
    keys = rows[0].keys()
    sup = {k: float(max(r[k] for r in rows)) for k in keys}
    i = int(np.argmax([r['K'] for r in rows]))
    lim = pieces(1e-9, k0)
    out[str(k0)] = dict(sup=sup, eta_at_supK=float(etas[i]), at_supK=rows[i], limit_eta0={k: float(v) for k, v in lim.items()},
                        C_strip=5*sup['K']+np.pi**2/4, C_strip_eta_to_0=5*lim['K']+np.pi**2/4)
    print(k0, json.dumps(out[str(k0)], indent=None)[:900])
out['C1_transfer'] = {str(e): C1_const(e) for e in [1e-6, 1e-3, 0.1, 0.5, 1.0, 1.5]}
print(out['C1_transfer'])
json.dump(out, open('constants_C149.json', 'w'), indent=1)

# transfer lemma (PROOF Lemma 4.5): C_2(eta) = max(2 eta T(0)/int_eta^inf T, C_1(eta) e^(kappa_eta eta)),
# kappa_eta = sqrt(2) max(1/sin eta, 2) (C146 Prop 3.1); Lambda_e = q (k0 + 1/2) e^(kappa_eta eta).
def transfer(eta, k0=0.4):
    te = np.tan(eta); q = 1-2*eta/np.pi
    T0 = 1/(np.pi*te)
    t = np.tanh(eta); omt = 2/(np.exp(2*eta)+1)
    tail_eta = 2*np.arctan(te*omt/(t+te*te))/(2*np.pi)
    ke = np.sqrt(2)*max(1/np.sin(eta), 2.0)
    C1 = C1_const(eta)
    C2 = max(2*eta*T0/tail_eta, C1*np.exp(ke*eta))
    Le = q*(k0+0.5)*np.exp(ke*eta)
    # c'(r) with delta/eta <= Le/lambda_* = r:  c' = e^(-ke eta r)/(1 + C2 r); need c' - (2k0+1) eta cot(eta)/(2 pi) >= 1/2 c'
    return dict(C1=C1, C2=C2, Lambda_e=Le, kappa_eta_eta=ke*eta)
tr = {str(e): transfer(e) for e in [1e-6, 1e-3, 1e-2, 0.1, 0.3, 0.5]}
print(tr)
out['transfer'] = tr
json.dump(out, open('constants_C149.json', 'w'), indent=1)
