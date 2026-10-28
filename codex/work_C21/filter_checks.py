"""C21 delayed Jacobi filters, exact moments and direct weight contraction.

Extends C19's formulas, rather than re-running its reflection-cell ED.
Only standard library, NumPy and SciPy; output stays beside this file.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
from fractions import Fraction
from math import comb, factorial, sqrt, exp, pi, cos, sin
import json
import numpy as np
import scipy.linalg as la

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def kron(xs):
    out = np.ones((1, 1), dtype=complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def filters(n, s):
    k = n+s
    raw, hs = [], []
    for j in range(n):
        q = [0]*(j+1)
        for l in range(j+1):
            a = comb(j+2*s, l)*comb(j+1, j-l)
            for r in range(j-l+1):
                q[l+r] += a*comb(j-l, r)*(-1)**(j-l-r)
        f = [0]*(k+1)
        for l, a in enumerate(q):
            for r in range(s+1):
                f[1+l+r] += a*comb(s, r)*(-1)**r
        raw.append(f[1:])
        hs.append(Fraction(j+1, (2*j+2*s+2)*(j+2*s+1)))
    coeff = np.array(raw, float).T/np.sqrt(np.array(hs, float))[None, :]
    return coeff, raw, hs


def exact_checks(n, s):
    _, raw, hs = filters(n, s)
    k = n+s
    for a in range(n):
        for r in range(s):
            assert sum((e+1)**r*raw[a][e] for e in range(k)) == 0
        assert sum((e+1)**s*raw[a][e] for e in range(k)) == (
            (-1)**s*factorial(s)*comb(a+2*s, a))
        for b in range(n):
            value = sum((Fraction(raw[a][e]*raw[b][f], e+f+2)
                         for e in range(k) for f in range(k)), Fraction())
            assert value == (hs[a] if a == b else 0)
    d = [(-1)**r*comb(s, r) for r in range(s+1)] + [0]*(n-1)
    g = sum((Fraction(d[e], k*(e+f+2))
             for e in range(k) for f in range(k)), Fraction())
    beta2s = Fraction(factorial(s-1), factorial(s+1))
    betaks = Fraction(factorial(k+1)*factorial(s-1), factorial(k+s+1))
    assert g == (beta2s-betaks)/k and g > 0
    return g


def weights(n, s):
    k, N = n+s, 2**n
    fs, _, _ = filters(n, s)
    gamma = []
    for a in range(n):
        gamma += [kron([Z]*a+[X]+[I2]*(n-a-1)),
                  kron([Z]*a+[Y]+[I2]*(n-a-1))]
    P = kron([Z]*n)
    eta = [1j*g@P for g in gamma]
    A1 = np.zeros((3, 3), complex)
    A1[0, 2] = A1[2, 0] = 1/sqrt(2)
    A2 = np.zeros((3, 3), complex)
    A2[0, 2], A2[2, 0] = 1j/sqrt(2), -1j/sqrt(2)
    alpha = exp(-1)
    js, mus = [], []
    for a in range(1, n+1):
        pl, pr = pi*alpha**(2*a), pi*alpha**(2*a-1)
        t, theta = 1+alpha**(2*a), -(pl+pr)/2
        js.append(-t/2*np.array([[sin(theta), cos(theta)],
                                [-cos(theta), sin(theta)]]))
        mus.append(2*t*cos((pr-pl)/2))
    F = sum((-mus[a]/2*kron([I2]*a+[Z]+[I2]*(n-a-1))
             for a in range(n)), np.zeros((N, N), complex))
    dr, dl = [], []
    for e in range(k):
        right = np.zeros((3*N, 3*N), complex)
        left = right.copy()
        for a in range(n):
            for u, A in enumerate((A1, A2)):
                right += fs[e, a]*np.kron(A, eta[2*a+u])
                op = sum((js[a][u, v]*gamma[2*a+v] for v in range(2)),
                         np.zeros((N, N), complex))
                left -= fs[e, a]*np.kron(A.T, op)
        dr.append(right)
        dl.append(left)
    d = np.array([(-1)**r*comb(s, r) for r in range(s+1)]+[0]*(n-1), float)
    g = float(exact_checks(n, s))
    rmax = sqrt(2*np.max(np.sum(fs**2, axis=1)))
    rho = min(1/(4*k*rmax), sqrt(g/(4*k*np.max(abs(d))*la.norm(F, 2))))
    MR = [np.eye(3*N)/k+rho*a for a in dr]
    ML = [np.eye(3*N)/k+rho*a for a in dl]
    for e in range(k):
        ML[e][N:2*N, N:2*N] -= rho*rho/g*d[e]*F
    K = np.zeros((N*N, N*N), complex)
    for a in range(n):
        for u in range(2):
            for v in range(2):
                K += js[a][u, v]*np.kron(eta[2*a+u], gamma[2*a+v])
    es = np.arange(1, k+1)
    C = 1/(es[:, None]+es[None, :])
    scalar = -3*np.sum(C)/k**2
    H2 = np.zeros_like(K)
    for e in range(k):
        for f in range(k):
            for a in range(3):
                for b in range(3):
                    H2 -= C[e, f]*np.kron(MR[e][a*N:(a+1)*N, b*N:(b+1)*N],
                                          ML[f][a*N:(a+1)*N, b*N:(b+1)*N])
    pred = scalar*np.eye(N*N)+rho*rho*(K+np.kron(np.eye(N), F))
    grade = np.kron(np.diag([-1, -1, 1]), P)
    # Hermitian off-diagonal syndrome components of the first nonzero moment.
    mr = sum((float(e+1)**s*MR[e] for e in range(k)), np.zeros_like(MR[0]))
    ml = sum((float(e+1)**s*ML[e] for e in range(k)), np.zeros_like(ML[0]))
    R1 = (mr[:N, 2*N:]+mr[2*N:, :N])/sqrt(2)
    L2 = (ml[:N, 2*N:]-ml[2*N:, :N])/(1j*sqrt(2))
    comm = R1@L2-L2@R1
    row = dict(n=n, s=s, poles=k, rho2=rho*rho,
               min_weight=min(la.eigvalsh(m)[0] for m in ML+MR),
               weight_sum_error=max(la.norm(sum(ms)-np.eye(3*N), 2) for ms in (ML, MR)),
               parity_error=max(la.norm(grade@m-m@grade, 2) for m in ML+MR),
               h2_error=la.norm(H2-pred, 2),
               h2_error_over_rho2=la.norm(H2-pred, 2)/(rho*rho),
               end_moment_commutator=la.norm(comm, 2))
    assert row['h2_error'] < 1e-12
    assert row['parity_error'] < 1e-12
    assert row['weight_sum_error'] < 1e-12
    assert row['min_weight'] >= 3/(4*k)-1e-12
    assert row['end_moment_commutator'] > 1e-12
    return row


def main():
    exact = [(n, s) for n in (1, 2, 3, 5, 8) for s in (1, 2, 4, 8)]
    for n, s in exact:
        exact_checks(n, s)
    rows = [weights(n, s) for n in (1, 2, 3) for s in (1, 2, 3)]
    result = dict(exact_cases=exact, numerical_weights=rows)
    Path(__file__).with_name('filter_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
