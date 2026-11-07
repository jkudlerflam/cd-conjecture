"""C146 constants. Evaluation of explicit analytic formulas (floating point; grid sups use monotone-piece
upper bounds, so they are rigorous up to rounding). J=1."""
import json
import numpy as np

# (W3'') edge bound: eta Z^2 >= pi Q/(4(1+kQ)), k = 16/(3 pi^2)
k = 16/(3*np.pi**2)
# c_Z = inf_Q pi Q/(4(1+kQ) tanh(2Q)); clean split bound and monotone-piece grid bound
cZ_clean = (np.pi/8)/(1+k/2)          # Q <= 1/2 uses tanh(2Q) <= 2Q, Q >= 1/2 uses tanh <= 1
Qg = np.geomspace(1e-6, 5, 400001)  # Q < 1e-6: tanh(2Q) <= 2Q gives c >= (pi/8)/(1+k 1e-6), h <= 2(1+k 1e-6)
lo, hi = Qg[:-1], Qg[1:]
# on [lo,hi]: Q/(1+kQ) increasing, tanh(2Q) increasing -> lower bound pi lo/(4(1+k lo) tanh(2 hi))
cZ_grid_lower = float(np.min(np.pi*lo/(4*(1+k*lo)*np.tanh(2*hi))))
# tail Q >= 5: pi Q/(4(1+kQ)) >= pi*5/(4(1+5k)) and tanh <= 1
cZ_tail = np.pi*5/(4*(1+5*k))
cZ = min(cZ_grid_lower, cZ_tail, (np.pi/8)/(1+k*1e-6))
# C_e = 4 sup_Q h(Q), h = tanh(2Q)(1+kQ)/Q; upper bound per piece tanh(2hi)(1+k hi)/lo
h_up = float(np.max(np.tanh(2*hi)*(1+k*hi)/lo))
h_tail = 1*(1+k*5)/5 + 0  # for Q >= 5: h <= (1+kQ)/Q <= 1/5 + k
Ce = 4*max(h_up, 1/5+k, 2*(1+k*1e-6))
Ce_clean = 4*2*(1+k/2)

# (W1) replacement: I(b) = int_0^inf min(4, b^2 e(4+e)/(1-e)^2) dd, e = 4 exp(-2bd)
def I_of_b(b):
    # crossover e* solves b^2 e(4+e) = 4(1-e)^2
    A, B_, C = b**2-4, 4*b**2+8, -4
    es = (-B_+np.sqrt(B_**2-4*A*C))/(2*A) if abs(A) > 1e-14 else -C/B_
    ss = 0.5*np.log(4/es)
    return (4*ss + (b**2/2)*(5/(1-es)-5+np.log(1-es)))/b
bs = np.linspace(1, 2, 2001)
Ib = np.array([I_of_b(b) for b in bs])
I_sup = float(Ib.max())
# check one value by quadrature
from scipy.integrate import quad
def I_quad(b):
    f = lambda d: min(4.0, b**2*(4*np.exp(-2*b*d))*(4+4*np.exp(-2*b*d))/(1-4*np.exp(-2*b*d))**2) if 4*np.exp(-2*b*d) < 1 else 4.0
    return quad(f, 0, 40, limit=400, points=[0.5*np.log(4)/b, 2/b])[0]
# V_W(eta) = a(4 I(b) - 8) + 12 a B, a = pi sin eta, B = pi/(2 eta), b = pi/(pi - eta)
etas = np.linspace(1e-6, np.pi/2, 200001)
bb = np.pi/(np.pi-etas)
VW = np.pi*np.sin(etas)*(4*np.array([I_of_b(x) for x in bb[::100]]).repeat(100)[:len(etas)]-8) + 12*np.pi**2*np.sin(etas)/(2*etas)
# conservative: use I_sup everywhere
VW_cons = np.pi*np.sin(etas)*(4*I_sup-8) + 6*np.pi**2*np.sin(etas)/etas
VW_sup = float(VW_cons.max())
VW_eta0 = 6*np.pi**2 + 0.0

# Theorem A.2 constants (C145): C_G1 analytic, A_tail, M0 = 2
kap = lambda d: 1-4*np.exp(-2*d)
ds = np.linspace(0.75, 3, 2251)
Atail = float(np.min(ds+2/kap(ds))); Acore = float(np.min(ds+1/kap(ds)))
CG1 = (8/9)*(np.pi/2+(3.61+4/np.pi)*Acore)
def K(CGS=1.0, VW=VW_sup, Ce=Ce, M0=2.0):
    terms = [np.pi*(CG1+CGS)/4, VW/6, 4*np.pi*M0*Atail, 16*np.pi/3, np.pi**2*Ce/8]
    return float(sum(terms)), [float(t) for t in terms]
Kf, terms = K()
out = dict(k=k, cZ_clean=cZ_clean, cZ_grid=cZ, Ce_grid=Ce, Ce_clean=Ce_clean, h_sup_upper=h_up,
           I_sup=I_sup, I_at_b1=float(I_of_b(1.0)), I_at_b2=float(I_of_b(2.0)), I_quad_b1=I_quad(1.0), I_quad_b15=I_quad(1.5), I_formula_b15=float(I_of_b(1.5)),
           VW_sup=VW_sup, VW_eta_to_0=VW_eta0, eta_argmax_VW=float(etas[np.argmax(VW_cons)]),
           A_tail=Atail, A_core=Acore, CG1_analytic=CG1,
           K_strip=Kf, K_terms=terms, C_strip=5*Kf+np.pi**2/4)
for CGS in [1.0, 5.0]:
    kk, _ = K(CGS=CGS); out['C_strip_CGS_%g' % CGS] = 5*kk+np.pi**2/4
json.dump(out, open('constants_C146.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

# Proposition G.1: (G) derivative half from the comparability half alone, with
# C_GS(eta) = (4/3)[2 eta log(M1/M0) + eta kappa |S|],  |S| <= 2 min_d [d + (1/b) log(2 M1/(kappa(d) M0))],
# eta kappa <= sqrt(2) max(eta/sin eta, 2 eta),  M1 = 4 C_D(eta) <= (16/pi)(log(1/eta) + pi/2) (C142 analytic majorant).
rows = []
dd = np.linspace(0.7, 4, 3301)
for eta in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9, 1e-12]:
    b = np.pi/(np.pi-eta); M0 = 2.0; M1 = (16/np.pi)*(np.log(1/eta)+np.pi/2)
    S = 2*np.min(dd+np.log(2*M1/(kap(dd)*M0))/b)
    ek = np.sqrt(2)*max(eta/np.sin(eta), 2*eta)
    CGS = (4/3)*(2*eta*np.log(M1/M0)+ek*S)
    kk, _ = K(CGS=CGS)
    rows.append(dict(eta=eta, M1=M1, S_len_bound=float(S), eta_kappa=float(ek), CGS_loglog=float(CGS), C_strip=float(5*kk+np.pi**2/4)))
out['PropG1_table'] = rows
json.dump(out, open('constants_C146.json', 'w'), indent=1)
for r in rows: print(r)
