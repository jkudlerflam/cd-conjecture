"""Constants of the C150 doubling bound, as functions of eta in (0, eta0].

All sups are of explicit one-variable functions; they are evaluated on a
log-uniform grid in s (from 1e-9 eta to 40) whose resolution is checked by
doubling it. Asymptotic values at s -> infinity are included analytically.
"""
import json
import numpy as np


def grid(eta, N):
    return np.concatenate([[0.0], np.geomspace(1e-9*eta, 40.0, N)])


def ratios_T(eta, s):
    c = np.cos(2*eta); C = np.cosh(2*s); S = np.sinh(2*s)
    d = C-c
    r1 = -2*S/d
    r2 = 4*(-(d)*C+2*S*S)/d**2
    return r1, r2


def ratios_a1(eta, s):
    # a1 ~ 1/(cosh 2x + cos eta): same algebra with c -> -cos eta
    c = -np.cos(eta); C = np.cosh(2*s); S = np.sinh(2*s)
    d = C-c
    r1 = -2*S/d
    r2 = 4*(-d*C+2*S*S)/d**2
    # third derivative: f = 1/d, f' = -2S/d^2, use exact formula via derivatives of 1/d
    # (1/d)''' = -(d'''/d^2) + 6 d' d''/d^3 - 6 d'^3/d^4, with d' = 2S, d'' = 4C, d''' = 8S
    r3 = (-(8*S)/d**2+6*(2*S)*(4*C)/d**3-6*(2*S)**3/d**4)*d
    return r1, r2, r3


def consts(eta, N=400000):
    s = grid(eta, N)
    r1, r2 = ratios_T(eta, s)
    tau1 = eta*max(np.max(np.abs(r1)), 2.0)
    tau1_exact = 2*eta/np.sin(2*eta)
    tau2 = eta**2*max(np.max(np.abs(r2)), 4.0)
    q1, q2, q3 = ratios_a1(eta, s)
    al1 = max(np.max(np.abs(q1)), 2.0)
    al2 = max(np.max(np.abs(q2)), 4.0)
    al3 = max(np.max(np.abs(q3)), 8.0)
    s1 = max(tau1, eta*al1); s2 = max(tau2, eta**2*al2)
    C3 = tau2+3*tau1*s1+s2+3*s1**2
    C3e = eta**2*al3+3*eta*al2*s1+al1*(s2+3*s1**2)
    # pair half-width bound delta: tanh(delta) = tan(eta)
    delta = np.arctanh(np.tan(eta))
    # H = (1/2)[1 + sup_{s>=0} T((s-delta)_+)/T(s)]
    c = np.cos(2*eta)
    ss = np.concatenate([np.linspace(0, 20*delta, 400001), np.linspace(20*delta, delta+12.0, 100001)])
    ratio = (np.cosh(2*ss)-c)/(np.cosh(2*np.maximum(ss-delta, 0))-c)
    H = 0.5*(1+np.max(ratio))
    return dict(eta=eta, tau1=tau1, tau1_exact=tau1_exact, tau2=tau2, al1=al1, al2=al2, al3=al3,
                s1=s1, s2=s2, C3=C3, C3e=C3e, delta=delta, delta_over_eta=delta/eta, H=H)


def vnorm():
    # v(tau) = int_{-1/2}^tau (w - 1/16), w = (1/4-|tau|)_+ ; sup |v|
    t = np.linspace(-0.5, 0.5, 2000001)
    w = np.maximum(0.25-np.abs(t), 0)-1/16
    v = np.concatenate([[0], np.cumsum((w[1:]+w[:-1])/2*np.diff(t))])
    return float(np.max(np.abs(v))), float(abs(v[-1]))


if __name__ == '__main__':
    out = {}
    vn, vend = vnorm()
    out['vnorm'] = vn; out['vnorm_exact'] = 9/512; out['v_end'] = vend
    rows = []
    for eta in [1e-6, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5]:
        r = consts(eta)
        r2 = consts(eta, N=800000)
        r['grid_change_C3'] = abs(r2['C3']-r['C3'])
        rows.append(r)
        print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}))
    out['rows'] = rows
    out['sup_over_eta_le_0.5'] = {k: max(r[k] for r in rows) for k in ['tau1', 'tau2', 's1', 's2', 'C3', 'C3e', 'H', 'delta_over_eta']}
    print(vn, 9/512, out['sup_over_eta_le_0.5'])
    json.dump(out, open('constants_C150.json', 'w'), indent=1)
