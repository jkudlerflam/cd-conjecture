"""Small independent controls for C48, not a microscopic CFT calculation."""
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import expm, eigvalsh
from scipy.integrate import quad


def logsumexp(xs):
    a = max(xs)
    return a + math.log(sum(math.exp(x-a) for x in xs))


def parameter_checks():
    alpha = 1.0
    L = 8*(alpha+math.log(32)+10)
    logq = math.log(16)-7*L/8+alpha+6
    eta = math.exp(1/8)/8
    rows = []
    for n in (1, 2, 3, 23, 4463):
        b = n+3
        ld = math.log(4464)
        logC0 = logsumexp([ld, b*ld, (b+1)*ld,
                           math.log(b)+(b-1)*ld,
                           math.log(2)+(b+2)*ld])
        T = max(2*L, math.log(2*4464), logC0+2)
        ratio = L/(16*T)
        logA = 22*math.log(2)+L+math.log(2+16/ratio**2)
        assert T >= 2*L and logC0-T <= -2+1e-10
        rows.append(dict(n=n, L=L, T=T, kappa_over_gap=ratio,
                         beta_min_times_gap=2*T,
                         log10_A_dec=logA/math.log(10),
                         input_support_exponent=L+(n+1)*math.log(128)+1))
    q = math.exp(logq)
    assert eta < 1 and q < .01 and 5*q/(1-q) < 1
    return dict(eta=eta, animal_ratio=q, KP_ratio=5*q/(1-q), rows=rows)


def embedding_check():
    Delta, u = 1., 2.
    E = np.array([Delta]*3+[Delta+u])
    a = np.diag(np.r_[0., E])
    Aold = np.diag(np.r_[np.zeros(4), E])
    flip = np.block([[np.zeros((4,4)), np.eye(4)],
                     [np.eye(4), np.zeros((4,4))]])
    J = np.zeros((20, 8))
    for s in range(4):
        J[s*5, s] = 1
        J[s+1, s+4] = 1
    Anew = np.kron(np.eye(4), a)
    assert np.linalg.norm(J.T@J-np.eye(8)) == 0
    assert np.linalg.norm(Anew@J-J@Aold) == 0
    valid = np.diag(J@J.T).astype(bool)
    invalid_E = np.diag(Anew)[~valid]
    rows=[]
    for lam in (.1, .05, .025):
        old = Aold+lam*flip
        new = Anew+lam*J@flip@J.T
        spec_error = np.max(np.abs(eigvalsh(new)-np.sort(np.r_[eigvalsh(old),invalid_E])))
        assert spec_error < 1e-12
        for beta in (2., 4., 10.):
            Zo = np.exp(-beta*eigvalsh(old)).sum()
            Zn = np.exp(-beta*eigvalsh(new)).sum()
            actual = math.log(Zn/Zo)/beta
            zin = 3*np.exp(-beta*E).sum()
            zval = 4+np.exp(-beta*E).sum()
            bound = math.log1p(math.exp(beta*lam)*zin/zval)/beta
            assert -1e-13 <= actual <= bound+1e-13
            rows.append(dict(lam=lam,beta=beta,pressure_difference=actual,bound=bound))
        low = (E-np.sqrt(E**2+4*lam**2))/2
        err = np.max(np.abs(low+lam**2/E))
        assert err <= 1.01*lam**4/Delta**3
    # Integrate the exact resummed operator kernels, not the divergent raw series.
    kernel_integrals = np.array([quad(lambda t: math.exp(-energy*t),0,np.inf)[0]
                                  for energy in E])
    expected = np.array([1,1,1,1/3])
    assert np.max(np.abs(kernel_integrals-expected)) < 1e-12
    return dict(max_kernel_error=float(np.max(np.abs(kernel_integrals-expected))),
                enlarged_dimension=20, original_dimension=8,
                raw_order_two_ratio=u/Delta, thermal_checks=rows)


def circle_components(labels):
    occupied = {i for i,x in enumerate(labels) if x != 0}
    result=[]
    while occupied:
        seed=min(occupied)
        comp={seed}
        stack=[seed]
        occupied.remove(seed)
        while stack:
            i=stack.pop()
            for j in ((i-1)%len(labels),(i+1)%len(labels)):
                if j in occupied:
                    occupied.remove(j)
                    comp.add(j)
                    stack.append(j)
        result.append(comp)
    return result


def polymer_factorization_check():
    # A one-site time circle with noncommuting active block operators.
    # Enumerate every active/passive/vacuum history and test rank-one reset
    # factorization, including fully winding histories, before taking log.
    A=np.diag([0., 1.3])
    P=np.diag([1.,0.])
    Q=np.eye(2)-P
    V=np.array([[.03, .08j],[-.08j,-.02]])
    widths=[.7,.9,.8,.75,.85,.95]
    free=[expm(-dt*A) for dt in widths]
    full=[expm(-dt*(A+V)) for dt in widths]
    active=[T-F for T,F in zip(full,free)]
    mats=[[P, F@Q, B] for F,B in zip(free,active)]
    def trace_word(labels):
        ans=np.eye(2,dtype=complex)
        for t,label in enumerate(labels):
            ans=mats[t][label]@ans
        return np.trace(ans)
    total=0j
    maxerr=0.
    for labels in itertools.product(range(3),repeat=len(widths)):
        weight=trace_word(labels)
        fact=1+0j
        for comp in circle_components(labels):
            isolated=[labels[t] if t in comp else 0 for t in range(len(widths))]
            fact*=trace_word(isolated)
        maxerr=max(maxerr,abs(weight-fact))
        total+=fact
    direct=np.eye(2,dtype=complex)
    for T in full:
        direct=T@direct
    direct=np.trace(direct)
    assert maxerr < 1e-13 and abs(total-direct) < 1e-12
    return dict(histories=3**len(widths), max_factorization_error=maxerr,
                partition_error=float(abs(total-direct)), partition=float(direct.real))


def holder_check():
    rng=np.random.default_rng(48)
    energies=np.array([0., .3, .9, 2., 5., 10.])
    beta=7.
    durations=np.array([.2,1.1,2.3,3.4])
    mats=[]
    prodnorm=1.
    for _ in durations:
        O=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
        O/=np.linalg.norm(O,2)
        mats.append(O)
        prodnorm*=np.linalg.norm(O,2)
    ans=np.eye(6,dtype=complex)
    for dt,O in zip(durations,mats):
        ans=ans@np.diag(np.exp(-dt*energies/2))@O
    lhs=abs(np.trace(ans))
    rhs=np.exp(-beta*energies/2).sum()*prodnorm
    assert lhs <= rhs
    return dict(lhs=float(lhs),rhs=float(rhs))


if __name__=='__main__':
    results=dict(parameters=parameter_checks(),embedding=embedding_check(),
                 polymer=polymer_factorization_check(),holder=holder_check())
    Path('checks.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(dict(status='all checks passed',
          kernel_error=results['embedding']['max_kernel_error'],
          polymer_error=results['polymer']['max_factorization_error'],
          partition_error=results['polymer']['partition_error'],
          KP_ratio=results['parameters']['KP_ratio']),indent=2))
