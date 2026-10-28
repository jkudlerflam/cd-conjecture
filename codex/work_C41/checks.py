"""Independent C41 checks. Python/numpy/scipy only; all output stays here."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import itertools
import json
from math import comb, pi, sin, cos
from pathlib import Path
import numpy as np
from scipy.linalg import eigh


def positions(L, M):
    xs = list(itertools.combinations(range(L), M))
    return xs, {x: i for i, x in enumerate(xs)}


def charge(L, M, theta):
    xs, ids = positions(L, M)
    H = np.zeros((len(xs), len(xs)), complex)
    for col, x in enumerate(xs):
        for a in x:
            for step in (-1, 1):
                c = (a + step) % L
                if c in x:
                    continue
                y = tuple(sorted((set(x) - {a}) | {c}))
                winding = 1 if a == L-1 and c == 0 else -1 if a == 0 and c == L-1 else 0
                H[ids[y], col] -= np.exp(1j * theta * winding)
    return H


def offset_test(L, M, b):
    # Distinct ordered colors 0,...,b-1, repeated. The label offset is r.
    xs, ids = positions(L, M)
    C = len(xs)
    H = np.zeros((C*b, C*b))
    for col, x in enumerate(xs):
        for r in range(b):
            for a in x:
                for step in (-1, 1):
                    c = (a+step) % L
                    if c in x:
                        continue
                    y = tuple(sorted((set(x)-{a}) | {c}))
                    winding = 1 if a == L-1 and c == 0 else -1 if a == 0 and c == L-1 else 0
                    H[b*ids[y]+(r-winding) % b, b*col+r] -= 1
    exact = eigh(H, eigvals_only=True)
    predicted = np.sort(np.concatenate([eigh(charge(L, M, 2*pi*r/b), eigvals_only=True) for r in range(b)]))
    err = float(np.max(np.abs(exact-predicted)))
    assert err < 1e-10
    assert exact[1]-exact[0] > 1e-10
    return dict(L=L, M=M, b=b, dimension=C*b, error=err, gap=float(exact[1]-exact[0]))


def projection_test(L, M, seed):
    q = 2
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(q*q,q*q)) + 1j*rng.normal(size=(q*q,q*q))
    h = (X+X.conj().T)/2
    words = list(itertools.product(range(q), repeat=M))
    wi = {w:i for i,w in enumerate(words)}
    S = len(words)
    R = np.zeros((S,S),complex)
    hs = [np.zeros((S,S),complex) for _ in range(M)]
    for col, w in enumerate(words):
        R[wi[(w[-1],)+w[:-1]],col] = 1
        for j in range(M):
            k = (j+1) % M
            for a in range(q):
                for c in range(q):
                    new = list(w); new[j] = a; new[k] = c
                    hs[j][wi[tuple(new)],col] += h[a*q+c,w[j]*q+w[k]]
    Hs = sum(hs)
    xs, ids = positions(L,M)
    C = len(xs)
    V = np.zeros((C*S,C*S),complex)
    for ix,x in enumerate(xs):
        block = sum(hs[j] for j in range(M-1) if x[j+1] == x[j]+1)
        if x[0] == 0 and x[-1] == L-1:
            block = block + hs[-1]
        V[ix*S:(ix+1)*S,ix*S:(ix+1)*S] = block
    cols, expected, branch_count = [], [], 0
    alpha = M/L - sin(pi*M/L)**2/(M*L*sin(pi/L)**2)
    for r in range(M):
        theta = 2*pi*r/M
        Pr = sum(np.exp(-1j*theta*k)*np.linalg.matrix_power(R,k) for k in range(M))/M
        vals, vecs = eigh((Pr+Pr.conj().T)/2)
        U = vecs[:, vals > .5]
        if not U.shape[1]:
            continue
        # Resolve degenerate seas by their analytic determinant wavefunctions.
        seas = [0] if 2*r < M else [1] if 2*r > M else [0,1]
        for ell in seas:
            ks = (2*pi*np.arange(ell,ell+M)-theta-pi*(M-1))/L
            f = np.array([np.linalg.det(np.exp(1j*np.outer(x,ks))/np.sqrt(L)) for x in xs])
            f /= np.linalg.norm(f)
            E = -2*sum(np.cos(ks))
            assert np.linalg.norm(charge(L,M,theta)@f-E*f) < 1e-10
            cols.append(np.kron(f[:,None],U))
            expected.append(alpha*U.conj().T@Hs@U)
            branch_count += 1
    P = np.concatenate(cols,axis=1)
    target = np.zeros((P.shape[1],P.shape[1]),complex)
    offset = 0
    for B in expected:
        n = len(B); target[offset:offset+n,offset:offset+n] = B; offset += n
    err = float(np.max(np.abs(P.conj().T@V@P-target)))
    assert np.max(np.abs(P.conj().T@P-np.eye(P.shape[1]))) < 1e-10
    assert err < 1e-10
    return dict(L=L,M=M,q=q,dimension=C*S,retained_dimension=len(target),branches=branch_count,error=err)


def flux_asymptotics():
    out=[]
    b=24; Q=.71; r=1
    theta=(Q+2*pi*r)/b
    for L in (96,192,384,768):
        M=L//2
        A=2*sin(pi*M/L)/sin(pi/L)
        exact=A*(1-cos(theta/L))
        prediction=2*pi*(2*sin(pi*M/L))/L*((r+Q/(2*pi))/b)**2
        out.append(dict(L=L,exact=exact,prediction=prediction,relative_error=abs(exact/prediction-1)))
    return out


def number_rounding():
    # A soluble reference with a fixed positive penalty in M != 0 mod b.
    # It is a diagnostic of the comparison, not the original h24 model.
    rho=np.sqrt(2)/3; b=24; gap=.01
    mu=-2*cos(pi*rho)
    out=[]
    for L in (4000,8000,16000,32000,64000):
        center=int(rho*L)
        candidates=range(center-2*b,center+2*b+1)
        energies=[-2*sin(pi*M/L)/sin(pi/L)-mu*M+(gap if M % b else 0) for M in candidates]
        winner=list(candidates)[int(np.argmin(energies))]
        M=b*round(rho*L/b)
        f=lambda x:-2*sin(pi*x)/pi-mu*x
        rounding=L*(f(M/L)-f(rho))
        bound=2*pi*sin(pi*rho)*b*b/(8*L)
        out.append(dict(L=L,winner=winner,compatible=winner%b==0,rounding_cost=rounding,quadratic_bound=bound))
    return out


if __name__ == '__main__':
    report={'offset_spectra':[offset_test(*p) for p in ((5,2,2),(7,4,2),(7,3,3),(8,4,4),(9,6,3))],
            'first_order':[projection_test(*p) for p in ((5,3,17),(6,4,29),(6,3,41))],
            'flux_asymptotics':flux_asymptotics(),'number_rounding_toy':number_rounding()}
    # Exact rational version of the safe density/scaling window.
    from fractions import Fraction as F
    assert F(4,5)*3 > 2
    assert F(4,5)/F(3,5)**2 == F(20,9) > 2
    report['rational_stability_checks']={'higher_payload_lower_bound':'12/5','nonzero_umklapp_lower_bound':'20/9'}
    Path('checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
