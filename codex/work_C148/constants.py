"""PROOF Section 7 constants as functions of lambda (J=1), exact expressions on a grid,
and the analytic majorant. Output constants.json."""
import numpy as np, json
def sval(lam, nl=400):
    if lam >= 1.36: return 0.5
    t = lambda k: np.exp(-2*np.abs(k)*lam); k = np.arange(1, int(80/lam)+20); Jk = t(k)/(1+t(k))
    ell = np.linspace(1e-4, np.pi, nl)
    I0 = lambda y: (y/2+(Jk[None,:]*np.sin(2*np.outer(y,k))/k[None,:]).sum(1))/np.pi  # int_0^y J
    ratio = 2*I0(ell)/(1-2*I0(ell/2))
    return 0.5+0.5*float(ratio.max())
def consts(lam, s=None):
    s = sval(lam) if s is None else s
    a, b = np.pi*np.sinh(lam), np.pi/lam; ab = a*b; M0 = (1+np.cosh(lam))/4
    g = 13.5*s
    K0 = 1/(np.pi*np.tanh(lam)); Kp1 = 4/(np.pi*np.sinh(2*lam))
    C = np.cosh(2*lam); c = (-C+np.sqrt(C*C+8))/2
    M1 = 2*np.sinh(2*lam)/np.pi*np.sqrt(1-c*c)/(C-c)**2
    D = g*Kp1+K0*s/2; D1 = M1*(4*g+s/2)
    cp = 2+4/np.expm1(np.pi**2/lam)
    kA = min(np.pi**2/(2*lam), cp); kw = min(np.pi**2/(2*lam), 1+np.log(2*cp))
    Vv = ab*(4+12*kw); Varg = 4*kA*(1+D1/(b*D)); Vf = ab*Varg+2*Vv
    ends = 2*s*ab*(2*g+4*min(np.pi/(2*b), np.pi/2)*D)
    curv = ab*s*s/4
    K = Vf/8+ends+curv
    return dict(lam=lam, Delta=float(np.cosh(lam)), s=s, g=g, D=D, D1=D1, kA=kA, kw=kw,
                core=Vf/8/M0, ends=ends/M0, curv=curv/M0, K_over_M0=K/M0, ab_over_M0=ab/M0)
lams = np.r_[np.geomspace(1e-4, 0.05, 30), np.linspace(0.05, 6, 600), np.geomspace(6, 200, 40)]
rows = [consts(l) for l in lams]
rows32 = [consts(l, 1.5) for l in lams]
best = max(rows, key=lambda r: r['K_over_M0']); best32 = max(rows32, key=lambda r: r['K_over_M0'])
out = dict(sup_exact=best, sup_with_s_3_2=best32,
           C_exact=2.5*best['K_over_M0']+np.pi**2/4, C_s32=2.5*best32['K_over_M0']+np.pi**2/4,
           s_max=max(r['s'] for r in rows), samples={f"{r['Delta']:.6g}": r['K_over_M0'] for r in rows[::40]})
json.dump(out, open('constants.json','w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ['sup_exact','C_exact','C_s32','s_max']}, indent=1, default=float))
