"""Explicit constants of the C139 Theorem A (J=1) as functions of gamma, Delta=cos gamma in [0,1)."""
import numpy as np, json
cr = 1+4/np.pi
def parts(g):
    B = np.pi/(2*g); s = np.sin(g); a = np.pi*s; b = np.pi/g; q = 1-2*g/np.pi
    ab = a*b; cot = np.cos(g)/s
    M1 = (2/np.pi)*(0.325/s**2+0.25/s)
    mu_m = 3*B/4
    D0 = B*cot/np.pi
    D1_0 = 2*M1*B
    Dp = D0 + (2/g)*cr*mu_m
    D1p = D1_0 + b*(2/g)*cr*mu_m
    Vf = 8*ab + 80*ab/np.pi + 4*a*D1p/Dp
    Eend = 2*(1.5*B)*ab*cr*(1.5*B + 4*cr*Dp/b)
    Ecurv = ab*(1+q)**2/2*mu_m**2/2
    Klarge = Vf/8 + Eend + Ecurv
    Lcut = max(4*g*Dp*np.cosh(1), 5*B)
    Ksmall = (Lcut+1)**2
    return dict(gamma=g, Delta=np.cos(g), B=B, Dp=Dp, D1p=D1p, Vf=Vf, Eend=Eend, Ecurv=Ecurv,
                Klarge=Klarge, Lcut=Lcut, Ksmall=Ksmall, K=max(Ksmall, Klarge))
def K(g): return parts(g)['K']
if __name__ == '__main__':
    out = {}
    for D in [0, 0.1, 0.5, 0.9, 0.99, 0.999]:
        p = parts(np.arccos(D)); out[str(D)] = p
        print(D, {k: round(v, 3) for k, v in p.items()})
    # uniform majorant on [0,1-delta]: sup_gamma K(gamma) * delta
    worst = 0
    for delta in np.logspace(-6, 0, 300):
        g0 = np.arccos(1-delta)
        gs = np.linspace(g0, np.pi/2, 400)
        Kmax = max(K(g) for g in gs)
        worst = max(worst, Kmax*delta)
    print('max over delta of delta*sup K =', worst)
    out['max_delta_times_supK'] = worst
    json.dump(out, open('constants.json', 'w'), indent=1)
