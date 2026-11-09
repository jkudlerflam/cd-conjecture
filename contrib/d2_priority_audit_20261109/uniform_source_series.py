"""Canonical (des Cloizeaux) effective Hamiltonian power series.

Use H/V = H0 + eps B + eps**2 C, eps=1/V. B contains real
pairing p and antisymmetric facilitated amplitude f. C contains hopping
a+ib and chemical potential m. No unpaired onsite field is present.
Recursively solve the invariant-graph Riccati equation, then conjugate
the Bloch Hamiltonian by the square root of its graph metric.
This is a finite periodic-ring algebra check, not an IR phase test.
"""
from pathlib import Path
import sys, json, argparse
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'codex/work_C186'))
import sw_general as sw


def multiply(A, B, order):
    out = [np.zeros((A[0].shape[0], B[0].shape[1]), complex) for _ in range(order+1)]
    for n in range(order+1):
        for k in range(n+1):
            out[n] += A[k] @ B[n-k]
    return out


def canonical_series(L, p, f, a, b, m, order=6):
    P = sw.pspace(L, True)
    Q = [s for s in range(2**L) if s not in P]
    B = sw.build_full(L, sw.bond_matrix(0, 0, p, -f, f), sw.site_matrix(0, 0), True)
    C = sw.build_full(L, sw.bond_matrix(0, a+1j*b, 0, 0, 0), sw.site_matrix(m, 0), True)
    bpq, bqp, bqq = B[np.ix_(P,Q)], B[np.ix_(Q,P)], B[np.ix_(Q,Q)]
    cpp, cpq, cqp, cqq = C[np.ix_(P,P)], C[np.ix_(P,Q)], C[np.ix_(Q,P)], C[np.ix_(Q,Q)]
    inv = 1 / np.array([sw.adj_count(s, L, True) for s in Q], float)
    X = [np.zeros((len(Q),len(P)),complex) for _ in range(order+1)]
    for n in range(1,order+1):
        res = bqp.copy() if n==1 else cqp.copy() if n==2 else np.zeros_like(bqp)
        res += bqq @ X[n-1]
        if n>=2:
            res += cqq @ X[n-2] - X[n-2] @ cpp
        for i in range(1,n):
            j=n-1-i
            if j>=1:
                res -= X[i] @ bpq @ X[j]
            j=n-2-i
            if j>=1:
                res -= X[i] @ cpq @ X[j]
        X[n] = -inv[:,None]*res
    HB = [np.zeros((len(P),len(P)),complex) for _ in range(order+1)]
    HB[2] += cpp
    for n in range(2,order+1):
        HB[n] += bpq @ X[n-1] + cpq @ X[n-2]
    # Graph metric M = 1 + X^dagger X. Its nonconstant part starts at eps^2.
    U = multiply([x.conj().T for x in X], X, order)
    identity = [np.zeros_like(U[0]) for _ in range(order+1)]
    identity[0] = np.eye(len(P), dtype=complex)
    powers=[identity,U]
    for k in range(2,order//2+1):
        powers.append(multiply(powers[-1],U,order))
    def binomial_power(alpha):
        out=[np.zeros_like(U[0]) for _ in range(order+1)]
        coef=1.0
        for k, power in enumerate(powers):
            if k: coef *= (alpha-k+1)/k
            for n in range(order+1): out[n] += coef*power[n]
        return out
    HC=multiply(multiply(binomial_power(.5),HB,order),binomial_power(-.5),order)
    pos={s:i for i,s in enumerate(P)}
    return HC,pos


def source(L=8,p=1.,f=.27386127875,a=-.075,b=.5175,m=1.675):
    HC,pos=canonical_series(L,p,f,a,b,m)
    z=[h[pos[1<<2],pos[0]] for h in HC]
    return {
        'L':L,'p':p,'f':f,'a':a,'b':b,'m':m,
        'scaled_coefficients':[[float(q.real),float(q.imag)] for q in z],
        'max_hermiticity_residual':max(float(np.max(np.abs(h-h.conj().T))) for h in HC),
    }


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path); args=ap.parse_args()
    cases=[dict(L=L) for L in (6,7,8)]
    cases += [dict(L=7,p=p,f=f,a=a,b=b,m=m) for p,f,a,b,m in [
        (1,.3,0,.5,0),(1,.3,.2,.5,0),(1,.3,0,.5,2),
        (2,.3,0,.5,0),(1,.6,0,.5,0),(1,.3,0,1,0),
        (1,.3,.2,.5,2),(1,.3,-.2,.5,0),(1,.3,.2,-.5,0),
        (1,0,.2,.5,2),(0,.3,.2,.5,2)]]
    rows=[source(**case) for case in cases]
    out={'description':'Coefficient of eps^n in H_eff/V; physical order is V^(1-n).', 'rows':rows}
    if args.out: args.out.write_text(json.dumps(out,indent=2)+'\n')
    else: print(json.dumps(out,indent=2))
