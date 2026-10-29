"""Small analytic-formula and quadrature checks for C106, no many-body numerics."""
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def entropy_kernel(x):
    y = np.abs(x)
    z = np.exp(-y)
    return np.log1p(z) + y * z / (1.0 + z)


def a_for_ratio(delta, r, t=1.0):
    target = r*r*(4*t*t+delta*delta)
    f = lambda a: 4*t*t*(1-1/a**2)-delta*delta*(a*a-1)**2/a**2-target
    a0 = 1/np.sqrt(1-r*r)
    return brentq(f, a0, a0+0.1, xtol=1e-14)


def entropy(delta, temp, a, t=1.0):
    # sin(k)=temp*x/(2t), using all four quarter periods.
    limit = 2*t/temp
    # Boltzmann tails beyond 80 + |mu/T| are negligible here.
    top = min(limit, 80 + abs(a*delta/temp))
    def integrand(x):
        root = np.sqrt(x*x+(delta/temp)**2)
        return (entropy_kernel(root-a*delta/temp)
                + entropy_kernel(root+a*delta/temp)) / np.sqrt(1-(x/limit)**2)
    kf_scaled = delta*np.sqrt(a*a-1)/temp
    value, error = quad(integrand, 0, top,
                        points=[kf_scaled] if kf_scaled < top else None,
                        epsabs=1e-10, epsrel=1e-10, limit=200)
    return temp*value/(np.pi*t)


def fermion_ops(n):
    z = np.diag([1., -1.])
    f = np.array([[0., 1.], [0., 0.]])
    out = []
    for j in range(n):
        op = np.array([[1.]])
        for k in range(n):
            op = np.kron(op, z if k < j else f if k == j else np.eye(2))
        out.append(op)
    return out


def reduced(vector, left):
    mat = vector.reshape(4,4)
    return mat @ mat.conj().T if left else mat.T @ mat.conj()


lines = ['# C106 checks', '', '[Numerical] Small matrix identities and one-dimensional quadratures only. These check the analytic proof; no infrared extrapolation of an interacting chain is used.', '']
t = 1.
r = 0.25
fs = fermion_ops(4)
sx = np.array([[0.,1.],[1.,0.]])
sz = np.diag([1.,-1.])
for delta in [1e-1, 1e-2, 1e-3]:
    a = a_for_ratio(delta,r,t)
    X = delta*(sx-a*np.eye(2))
    B = -1j*t*sz
    one = np.block([[X/2,B],[B.conj().T,X/2]])
    h = sum(one[i,j]*fs[i].conj().T@fs[j] for i in range(4) for j in range(4))
    vals, vecs = np.linalg.eigh(h)
    M = np.sqrt(4*t*t+delta*delta)
    residual = abs((vals[-1]-vals[0])/2-M)
    marginal = max(np.linalg.norm(reduced(vecs[:,j],True)-reduced(vecs[:,j],False)) for j in [0,-1])
    velocity = 2*t*np.sqrt(1-1/a**2)*np.sqrt(1-delta*delta*(a*a-1)/(4*t*t))
    assert residual < 1e-12 and marginal < 1e-12
    assert abs(velocity/M-r) < 1e-12
    lines.append(f'[Numerical] delta={delta:g}: bond half-spread error {residual:.3g}, extremal-state equal-marginal residual {marginal:.3g}, v/M={velocity/M:.12g}.')
lines += ['', '| delta | T/delta | s(T) / [2 pi T/(3v)] |', '| :-- | --: | --: |']
for delta in [1e-2,1e-3,1e-4]:
    a = a_for_ratio(delta,r,t)
    v = r*np.sqrt(4*t*t+delta*delta)
    for factor in [1e-4,1e-3,0.01,1.,10.]:
        T=delta*factor
        ratio=entropy(delta,T,a,t)/(2*np.pi*T/(3*v))
        lines.append(f'| {delta:g} | {factor:g} | {ratio:.10f} |')
        if factor == 1e-4:
            assert abs(ratio-1) < 1e-3
        if factor == 10.:
            assert ratio < 0.5
lines += ['', '[Numerical] Additional one-temperature and entropy-increment checks at r=0.025, epsilon=1e-8.']
r_small = 0.025
delta_small = 1e-8
a_small = a_for_ratio(delta_small, r_small, t)
M_small = np.sqrt(4*t*t+delta_small*delta_small)
v_small = r_small*M_small
w_small = 0.75*v_small
for kappa in [0.1, 1., 10.]:
    for p in [1., 1.5, 2.]:
        T = kappa*M_small*(w_small/M_small)**p
        quotient = entropy(delta_small,T,a_small,t)/(2*np.pi*T/(3*v_small))
        assert quotient < 0.5
    T = kappa*w_small*w_small/M_small
    increment = entropy(delta_small,T,a_small,t)-entropy(delta_small,T/2,a_small,t)
    quotient = increment/(2*np.pi/(3*v_small)*(T/2))
    assert 0 < quotient < 0.5
    lines.append(f'[Numerical] kappa={kappa:g}: all tested E_p ratios below 1/2; increment ratio {quotient:.10f}, also below 1/2.')
lines += ['', '[Proved check] In the main r=1/4 family A/M = 2r = 1/2, so the free bound A <= 2M is respected.', '', '[Numerical] All assertions passed.']
Path('CHECKS.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
