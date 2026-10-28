"""[numerical] Small independent controls for the C46 proofs.
Uses only numpy, scipy and the standard library. No wrapper simulation.
"""
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import logsumexp


def embed(op, sites, count, local_d):
    dim = local_d ** count
    out = np.zeros((dim, dim), complex)
    labels = [np.unravel_index(i, (local_d,) * count) for i in range(dim)]
    rest = [j for j in range(count) if j not in sites]
    for i, a in enumerate(labels):
        ai = np.ravel_multi_index(tuple(a[j] for j in sites), (local_d,) * len(sites))
        for j, b in enumerate(labels):
            if all(a[k] == b[k] for k in rest):
                bi = np.ravel_multi_index(tuple(b[k] for k in sites), (local_d,) * len(sites))
                out[i, j] = op[ai, bi]
    return out


def cp_extend_two(op):
    units = {}
    for a in range(2):
        for b in range(2):
            x = np.zeros((3, 3), complex)
            x[a, b] = 1
            if a == b:
                x[2, 2] = .5
            units[a, b] = x
    out = np.zeros((9, 9), complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    out += op[2*a+c, 2*b+d] * np.kron(units[a,b], units[c,d])
    return out


def pressure_checks():
    q = np.diag([0., 0., 1.])
    p = np.eye(3) - q
    pp = np.kron(p, p)
    qq = np.eye(9) - pp
    raw = np.sin(np.arange(81).reshape(9, 9) * .37)
    raw = (raw + raw.T) / 2
    w = pp @ raw @ pp + qq @ raw @ qq
    w *= .07 / np.linalg.norm(w, 2)
    # Nonzero transitions between states with a different number of defects.
    nq2 = np.kron(q, np.eye(3)) + np.kron(np.eye(3), q)
    changing = np.linalg.norm(w @ nq2 - nq2 @ w, 2)
    code_indices = [0, 1, 3, 4]
    wp = w[np.ix_(code_indices, code_indices)]
    ew = cp_extend_two(wp)
    rows = []
    delta0 = 1.7
    for m in (2, 3, 4):
        bonds = [(j, j+1) for j in range(m-1)]
        if m > 2:
            bonds.append((m-1, 0))
        degrees = [sum(j in b for b in bonds) for j in range(m)]
        local_j = max(degrees) * np.linalg.norm(w, 2)
        delta = delta0 - 2*local_j
        nq = sum(embed(q, [j], m, 3) for j in range(m))
        D = delta0*nq + sum(embed(w, b, m, 3) for b in bonds)
        EH = sum(embed(ew, b, m, 3) for b in bonds)
        HP = sum(embed(wp, b, m, 2) for b in bonds)
        lower = np.linalg.eigvalsh(D-EH-delta*nq).min()
        ed = np.linalg.eigvalsh(D)
        ep = np.linalg.eigvalsh(HP)
        for beta in (.2, 1., 4.):
            ratio_log = float(logsumexp(-beta*ed)-logsumexp(-beta*ep))
            upper = m*math.log1p(.5*math.exp(-beta*delta))
            assert lower >= -1e-12
            assert -1e-12 <= ratio_log <= upper+1e-12
            rows.append(dict(M=m, beta=beta, log_ratio=ratio_log,
                             log_upper=upper, operator_slack=float(lower)))
    assert changing > 1e-3
    return dict(defect_number_commutator=float(changing), cases=rows)


def thermal_counterexample():
    beta, delta, u = 4., 1., 2.
    a = 1 + math.exp(-beta*delta)
    b = 1 + math.exp(-beta*(delta+u))
    cov = a*(b-a)/(3*a+b)**2
    with localcontext() as ctx:
        ctx.prec = 120
        r = Decimal(-4).exp()
        alpha = Decimal(-8)
        p = [1+r]
        for m in range(1, 81):
            p.append(r * alpha**m / Decimal(math.factorial(m)))
        c = [Decimal(0)]
        for m in range(1, 81):
            c.append((Decimal(m)*p[m] - sum(Decimal(k)*c[k]*p[m-k]
                       for k in range(1,m))) / (Decimal(m)*p[0]))
        partial = {str(N): str(sum(abs(x) for x in c[1:N+1])) for N in (10,20,40,80)}
    return dict(covariance=cov, scalar_taylor_radius=math.sqrt(1+(math.pi/4)**2)/2,
                absolute_coefficient_partial_sums=partial,
                vacuum_excursion_coefficients=[2**m for m in range(1,9)])


def free_thermal_check():
    A = np.array([[1., .18], [.18, .35]])
    B = np.array([[.12, .025], [.025, .015]])
    C = np.diag([.07, .013])
    def symbol(k):
        return -2*A*np.cos(k)+2*B*np.sin(k)+C
    def derivative(k):
        return 2*A*np.sin(k)+2*B*np.cos(k)
    grid = np.linspace(-math.pi, math.pi, 1001)
    roots=[]
    for band in range(2):
        f=lambda k: np.linalg.eigvalsh(symbol(k))[band]
        for a,b in zip(grid[:-1], grid[1:]):
            if f(a)*f(b)<0:
                k=brentq(f,a,b,xtol=1e-14)
                vals, vecs=np.linalg.eigh(symbol(k))
                z=vecs[:,np.argmin(abs(vals))]
                vel=float(z @ derivative(k) @ z)
                roots.append((k,vel))
    roots.sort()
    assert len(roots)==4
    coefficient=math.pi/12*sum(1/abs(v) for _,v in roots)
    points=[k for k,_ in roots]
    rows=[]
    for temperature in (.02,.01,.005,.0025):
        def f(k):
            return np.logaddexp(0,-abs(np.linalg.eigvalsh(symbol(k)))/temperature).sum()
        integ=quad(f,-math.pi,math.pi,points=points,epsabs=1e-12,epsrel=1e-11,limit=200)[0]
        measured=integ/(2*math.pi*temperature)
        rows.append(dict(T=temperature, minus_free_energy_change_over_T2=measured,
                         predicted=coefficient, error=measured-coefficient))
    assert abs(rows[-1]['error'])<2e-5
    return dict(roots_and_signed_velocities=roots, asymptotics=rows)


if __name__ == '__main__':
    result={'status':'numerical, independent toy controls only',
            'pressure':pressure_checks(),
            'thermal_counterexample':thermal_counterexample(),
            'free_thermal':free_thermal_check()}
    Path('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
