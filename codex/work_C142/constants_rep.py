"""Constants of PROOF.md Part I (repulsive side, Theorem I.9), as explicit functions of gamma in (0,pi/2].
Every line is a formula from PROOF.md; the script only evaluates them and takes sup over a fine gamma grid,
and also checks the analytic majorants used in the text. J=1."""
import json
import numpy as np

def consts(g):
    pi = np.pi
    q = 1-2*g/pi
    qJ = (pi-2*g)/(2*(pi-g))
    s = 1-g/pi                               # sawtooth constant, Lemma I.5
    a = pi*np.sin(g); b = pi/g; ab = a*b
    inv = (1+qJ)/(1-qJ)                      # restricted inverse, C139 (3.1)
    gG = (2+q)*inv*s*(1+qJ*q/2)              # L||G||, Lemma I.6
    D = 2*gG/(pi*g)                          # density error bound (>= 2 cot(g) g/pi), Lemma I.7
    sn = np.sin(g)
    M1 = (2/pi)*(0.325/sn**2+0.25/sn)        # max|K'|, C139 Lemma 5.0
    D1 = 4*M1*gG                             # derivative error, Lemma I.7
    cr = 1+4/pi
    Vv = (2+40/pi)*ab                        # Var(v_Q), C139 (3.2)
    Vf = ab*(4+4*D1/(b*D)) + 2*Vv            # Lemma I.8 (core)
    core = Vf/8
    ends = 2*s*ab*(cr*pi*D/b + 2*gG)         # Lemma I.8 (ends)
    a2 = ab*(1+q)**2/2
    curv = a2*(qJ*s)**2/2                    # mismatch curvature
    K = core+ends+curv
    L0 = int(np.ceil(4*g*D))                 # A exists for L >= L0
    return dict(q=q, qJ=qJ, s=s, g=gG, D=D, D1=D1, Vf=Vf, core=core, ends=ends, curv=curv,
                K_large=K, L0=L0, theta1=1-2*(2-np.log(2))*(pi-g)/pi**2,
                lam=1/(2*(pi-g)), lam_needed=(2-np.log(2))/pi**2)

if __name__ == '__main__':
    gs = np.linspace(1e-6, np.pi/2, 200001)
    rows = [consts(g) for g in gs]
    keys = ['g', 'core', 'ends', 'curv', 'K_large', 'L0', 'theta1']
    sup = {k: float(max(r[k] for r in rows)) for k in keys}
    sup['theta1_min'] = float(min(r['theta1'] for r in rows))
    sup['lam_min_over_needed'] = float(min(r['lam']/r['lam_needed'] for r in rows))
    sup['gamma_sin_M1_max'] = float(max(g*np.sin(g)*(2/np.pi)*(0.325/np.sin(g)**2+0.25/np.sin(g)) for g in gs))
    sup['L0_max'] = int(max(r['L0'] for r in rows))
    sup['sector_doubling_over_M0'] = max(5*sup['K_large'], 4*(sup['L0_max']-1))
    sup['C_total'] = sup['sector_doubling_over_M0'] + np.pi**2/4
    for gg in [np.pi/2, 1.0, 0.5, 0.177, 0.05, 1e-3, 1e-6]:
        c = consts(gg); sup[f'K_large_at_gamma_{gg:.4g}'] = c['K_large']
    print(json.dumps(sup, indent=1))
    json.dump(sup, open('constants_rep.json', 'w'), indent=1)
