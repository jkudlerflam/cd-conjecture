"""C145 constants. Floating point evaluation of explicit analytic formulas. J=1."""
import json
import numpy as np
from check_attractive import bounds

kap = lambda d: 1-4*np.exp(-2*d)            # proved (W2) at distance d, b >= 1
ds = np.linspace(0.75, 3, 2251)
dcore = ds[np.argmin(ds+1/kap(ds))]; dtail = ds[np.argmin(ds+2/kap(ds))]
Acore = dcore+1/kap(dcore); Atail = dtail+2/kap(dtail)

# (1) C138 Theorem B with proved (W2) at optimized split distances, C138 values for the rest
CD, CD1, CW, Ce = 0.5, 2.0, 8.0, 1.3
KB = (2*np.pi/(9*CD))*(CD+CD1*Acore)+np.pi*CW/12+16*np.pi*CD*Atail+16*np.pi/3+np.pi**2*Ce/8
KB_c138 = (2*np.pi/(9*CD))*(CD+CD1*(1+1/0.7))+np.pi*CW/12+16*np.pi*CD*(1+2/0.7)+16*np.pi/3+np.pi**2*Ce/8

# (2) proved part of (G) on K1 = {eta L r >= 4 C_D(eta)}: eta Var_K1(r/sigma) <= (8/(9 C_D))[eta C_D + C_D1 Acore]
etas = np.geomspace(1e-7, np.pi/2-1e-3, 400)
G1 = []
for eta in etas:
    bd = bounds(eta)
    G1.append((8/(9*bd['CD']))*(eta*bd['CD']+bd['CD1']*Acore))
G1 = np.array(G1)
nan_count = int(np.isnan(G1).sum()); CG1 = float(np.nanmax(G1)); eta_argmax = float(etas[np.nanargmax(G1)])

def Kstar_G(CGS, M0=2.0, CW=8.0, Ce=1.3):
    return np.pi*(CG1+CGS)/4+np.pi*CW/12+4*np.pi*M0*Atail+16*np.pi/3+np.pi**2*Ce/8

out = dict(d_core=float(dcore), A_core=float(Acore), d_tail=float(dtail), A_tail=float(Atail),
           kappa_d1=float(kap(1.0)), KB_with_proved_W2=float(KB), KB_C138=float(KB_c138),
           C_U1_B=float(5*KB+np.pi**2/4),
           CG1_sup=CG1, CG1_eta_argmax=eta_argmax, CG1_nan_points=nan_count, CG1_at_eta_min=float(G1[0]), CG1_at_eta_max=float(G1[-1]))
for CGS in [1.0, 2.0, 5.0]:
    k = Kstar_G(CGS); out['Kstar_G_CGS_%g' % CGS] = float(k); out['C_U1_G_CGS_%g' % CGS] = float(5*k+np.pi**2/4)
json.dump(out, open('constants_C145.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

# (3) final assembly: analytic C_G1 = 13.18, (W1') -> C_W = 4 pi, (W3') with c_Z = 3/4 -> C_e = 4 pi/3, C_GS = 1, M0 = 2
CG1_an = (8/9)*(np.pi/2+(3.61+4/np.pi)*Acore)
def K_final(CGS=1.0, M0=2.0, CW=4*np.pi, Ce=4*np.pi/3, CG1=CG1_an):
    terms = [np.pi*(CG1+CGS)/4, np.pi*CW/12, 4*np.pi*M0*Atail, 16*np.pi/3, np.pi**2*Ce/8]
    return float(sum(terms)), [float(t) for t in terms]
K, terms = K_final()
out.update(CG1_analytic=float(CG1_an), K_final=K, K_final_terms=terms, C_strip_final=5*K+np.pi**2/4,
           C_strip_final_numerical_inputs=5*K_final(CGS=0.002, CW=8.0, Ce=np.pi/0.808, CG1=CG1)[0]+np.pi**2/4)
json.dump(out, open('constants_C145.json', 'w'), indent=1)
print('final', json.dumps({k: out[k] for k in ['CG1_analytic', 'K_final', 'K_final_terms', 'C_strip_final', 'C_strip_final_numerical_inputs']}))
