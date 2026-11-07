"""C147 referee: recompute every constant of C142 Part I and C145 Parts A, C, E. Floating point."""
import json
import numpy as np
from scipy.integrate import quad
pi = np.pi; out = {}
# ---- C142 Lemma I.1: f(t) = atan(t)(1+t^2)/t increasing, f(1)=pi/2
t = np.linspace(1e-6, 1, 100001); f = np.arctan(t)*(1+t**2)/t
out['I1_f_monotone'] = bool(np.all(np.diff(f) > 0)); out['I1_f1'] = float(f[-1])
for Q, bb in [(0.0, 3.), (0.7, 3.), (2.0, 10.)]:
    lhs = quad(lambda x: 1/np.cosh(bb*x), Q, np.inf)[0]; rhs = (pi/(2*bb))/np.cosh(bb*Q)
    out['I1_check_Q%g_b%g' % (Q, bb)] = lhs/rhs
# ---- C142 Part I constants
gs = np.linspace(1e-7, pi/2, 400001)
q = 1-2*gs/pi; qJ = (pi-2*gs)/(2*(pi-gs)); s = 1-gs/pi; ab = pi**2*np.sin(gs)/gs
gG = (2+q)*(1+qJ)/(1-qJ)*s*(1+qJ*q/2)
D = 2*gG/(pi*gs); M1 = (2/pi)*(0.325/np.sin(gs)**2+0.25/np.sin(gs)); D1 = 4*M1*gG
cr = 1+4/pi; Vf = ab*(4+4*D1/((pi/gs)*D)) + 2*(2+40/pi)*ab
core = Vf/8; ends = 2*s*ab*gG*(2*cr/pi+2); curv = ab*(1+q)**2/2*(qJ*s)**2/2
K = core+ends+curv
theta1 = 1-2*(2-np.log(2))*(pi-gs)/pi**2
out.update(g_sup=float(gG.max()), g_monotone_decreasing=bool(np.all(np.diff(gG) <= 1e-12)),
           core_sup=float(core.max()), core_at_gamma_argmax=float(gs[core.argmax()]), core_at_0=float(core[0]),
           ends_sup=float(ends.max()), curv_sup=float(curv.max()), K_sup=float(K.max()),
           K_argmax=float(gs[K.argmax()]), theta1_min=float(theta1.min()),
           L0_max=float(np.ceil((4*gs*D).max())), gsinM1_max=float((gs*np.sin(gs)*M1).max()),
           ratio_I2=float(((1/(2*(pi-gs)))/((2-np.log(2))/pi**2)).min()))
Kan = (8*pi**2+80*pi+8*pi**2*0.575)/8 + 2*pi**2*(45/4)*(2*cr/pi+2) + pi**2/4
out['K_analytic'] = Kan; out['C_rep_analytic'] = 5*815+pi**2/4; out['C_rep_grid'] = 5*float(K.max())+pi**2/4
out['small_L_seam'] = 4*28
# ---- C145 constants
kap = lambda d: 1-4*np.exp(-2*d)
u = 2-np.sqrt(3); dc = 0.5*np.log(4/u); Ac = dc+1/kap(dc)
u = 3-2*np.sqrt(2); dt = 0.5*np.log(4/u); At = dt+2/kap(dt)
et = np.linspace(1e-6, pi/2, 200001)
r1 = 4*(9/(16*np.sqrt(3))*et/np.sin(et)+0.25*et)
ratio_an = float(r1.max()) + 4/pi
CG1_an = (8/9)*(pi/2+ratio_an*Ac)
out.update(A_core=Ac, d_core=dc, A_tail=At, d_tail=dt, kappa1=kap(1.0), ratio_CD1_CD_analytic=ratio_an, CG1_analytic=CG1_an)
def KG(CG1, CGS=1.0, M0=2.0, CW=4*pi, Ce=4*pi/3):
    T = [pi*(CG1+CGS)/4, pi*CW/12, 4*pi*M0*At, 16*pi/3, pi**2*Ce/8]
    return sum(T), T
k, T = KG(CG1_an); out.update(K_G=k, K_G_terms=T, C_strip=5*k+pi**2/4)
k13, _ = KG(13.19); out['K_G_with_13.19'] = k13; out['C_strip_with_13.19'] = 5*k13+pi**2/4
# C138 (4.3)
def KB(CD, CD1, kc, kt, CW, Ce):
    return [(2*pi/(9*CD))*(CD+CD1*kc), pi*CW/12, 16*pi*CD*kt, 16*pi/3, pi**2*Ce/8]
T = KB(0.5, 2.0, 1+1/0.7, 1+2/0.7, 8.0, 1.3); out.update(C138_terms=T, C138_K=sum(T), C138_C=5*125+pi**2/4)
T = KB(0.5, 2.0, Ac, At, 8.0, 1.3); out.update(C138_provedW2_K=sum(T), C138_provedW2_C=5*sum(T)+pi**2/4)
# C142 II.4 tail coefficient and Theorem II.6 constants
out['C142_II4_coef'] = 16*pi*(1+2/0.7); out['C142_II6_CD1_const'] = float(((8/pi)*(0.32476*(et/np.sin(et))**2+0.25*et**2/np.sin(et))).max())
json.dump(out, open('check_constants.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
