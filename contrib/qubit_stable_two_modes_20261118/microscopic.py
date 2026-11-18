"""Validate finite matching against literal microscopic qubit rings.

Scalar avoided-crossing calibration only, at L=12,14 with one BLAS worker.
The L-dependent crossing fields here define finite matrix elements; they
are not proposed thermodynamic Hamiltonians or central-charge evidence.
"""
from __future__ import annotations

from quartet import Sector, match, np, la, sp, threadpool_limits
from scipy.optimize import brentq
from pathlib import Path
import json
import time


def matrices(L, blockade=False):
    sectors = [Sector(L, n, 0, blockade=blockade) for n in range(0, L+1, 2)]
    bases = [s.sigma_basis(1) for s in sectors]
    sizes = [b.shape[1] for b in bases]
    ends = np.cumsum([0]+sizes)
    dim = int(ends[-1])
    pair = np.zeros((dim, dim), complex)
    number = np.zeros(dim)
    for i in range(len(sectors)):
        number[ends[i]:ends[i+1]] = sectors[i].N
        if i+1 < len(sectors):
            pair[ends[i+1]:ends[i+2], ends[i]:ends[i+1]] = (
                bases[i+1].conj().T @ sectors[i].pair(sectors[i+1], 3 if blockade else 1)
                @ bases[i]).toarray()
    def reduced(**kwargs):
        return sp.block_diag([(U.conj().T @ s.hamiltonian(**kwargs) @ U)
                              for s,U in zip(sectors,bases)]).toarray()
    if blockade:
        kinetic = reduced(kappa=0)
        hop = reduced(kappa=1)-kinetic
        interaction = np.zeros_like(kinetic)
    else:
        interaction = reduced(kappa=0, t2=0, V1=1, V2=0, V3=0)
        hop = reduced(kappa=1, t2=0, V1=0, V2=0, V3=0)
        kinetic = np.zeros_like(hop)
    return kinetic, hop, interaction, pair+pair.conj().T, number


def avoided(L, kappa, epsilons, blockade=False):
    kin,hop,inter,pair,number = matrices(L,blockade)
    reference = match(L,0,kappa)
    mu0 = reference['matching_mu']
    predicted = (16/9)*reference['amplitude_abs']
    rows=[]
    for eps in epsilons:
        t = eps**2/2
        if blockade:
            # The same units as the leading effective Hamiltonian.
            base = kin+kappa*hop+(4*eps/3)*pair
        else:
            base = (inter+kappa*t*hop+eps*pair)/t
            base -= np.diag(6*number)
        def gap(mu, vectors=False):
            H=base-np.diag(mu*number)
            if vectors:
                vals,vecs=la.eigh(H,subset_by_index=(0,2))
                return vals, vecs
            vals=la.eigvalsh(H,subset_by_index=(0,1))
            return float(vals[1]-vals[0])
        def derivative(mu):
            _, vecs = gap(mu,True)
            return float(np.vdot(vecs[:,0],number*vecs[:,0]).real
                         -np.vdot(vecs[:,1],number*vecs[:,1]).real)
        crossing,opt=brentq(derivative,mu0-12*eps**2,mu0+12*eps**2,
                            xtol=2e-14,rtol=1e-14,full_output=True)
        vals,vecs=gap(crossing,True)
        min_gap=float(vals[1]-vals[0])
        rows.append(dict(epsilon=eps, mu_crossing=float(crossing),
                         gap_over_t=min_gap,
                         coefficient_gap_over_2_t_eps2=min_gap/(2*eps**2),
                         coefficient_gap_over_2_t_eps4=min_gap/(2*eps**4),
                         ratio_to_prediction=(min_gap/(2*eps**2)/predicted
                                              if predicted>1e-10 else None),
                         third_level_gap_over_t=float(vals[2]-vals[1]),
                         mean_numbers=[float(np.vdot(vecs[:,i],number*vecs[:,i]).real)
                                       for i in range(2)],
                         optimization_success=bool(opt.converged),
                         evaluations=int(opt.function_calls)))
    return dict(L=L,kappa=kappa,blockade_only=blockade,dimension=len(number),
                reference_amplitude=reference['amplitude_abs'],
                predicted_coefficient=predicted,rows=rows)


def fourth_order_dark():
    """Fourth-order degenerate matching with its chemical/energy counterterms."""
    L=12
    kin,hop,_,B,N=matrices(L,True)
    H=kin+hop+2*np.diag(N)
    # Fixed-number eigenstates define the external dressed basis.
    i0=np.where(N==0)[0]
    i4=np.where(N==4)[0]
    _,v4,_=__import__('quartet').ground(H[np.ix_(i4,i4)])
    P=np.zeros((len(N),2),complex)
    P[i0[0],0]=1
    P[i4,1]=v4
    Q=la.null_space(P.conj().T)
    R=la.inv(Q.conj().T@H@Q)
    W=Q.conj().T@B@P
    Bq=Q.conj().T@B@Q
    Nq=Q.conj().T@np.diag(N)@Q
    M2=-W.conj().T@R@W
    alpha=float(((M2[1,1]-M2[0,0])/4).real)
    e2=float(M2[0,0].real)
    M3=W.conj().T@R@Bq@R@W
    M4=-W.conj().T@(R@Bq@R@Bq@R+alpha*R@Nq@R+e2*R@R)@W
    return dict(dimension=len(N), Q_gap=float(la.eigvalsh(Q.conj().T@H@Q)[0]),
                dark_pair_source_norm=float(la.norm(B@P[:,1])),
                second_order_diagonal=np.diag(M2).real.tolist(),
                second_order_off_diagonal_abs=float(abs(M2[1,0])),
                third_order_max_abs=float(np.max(abs(M3))),
                chemical_counterterm_per_gamma2=alpha,
                energy_counterterm_per_gamma2=e2,
                fourth_order_off_diagonal_abs=float(abs(M4[1,0])),
                predicted_gap_over_2_t_eps4=float(abs(M4[1,0]))*(4/3)**4)


def microscopic_dark_series():
    """Energy-dependent two-stage elimination through physical epsilon^6.

    H=D+epsilon A+epsilon^2 B+epsilon^4 C. High-energy D costs are
    retained exactly. The second elimination retains every other blockaded
    state. Counterterms follow from degeneracy, not an adjusted amplitude.
    """
    L=12
    _,hop,D,A,N=matrices(L,False)
    B=hop/2-2*np.diag(N)  # kappa=1, mu_eff=-2
    p=np.flatnonzero(np.diag(D).real<1e-10)
    q=np.flatnonzero(np.diag(D).real>1e-10)
    def pp(X): return X[np.ix_(p,p)]
    def qq(X): return X[np.ix_(q,q)]
    Apq=A[np.ix_(p,q)]; Bpq=B[np.ix_(p,q)]
    Aq=qq(A); Bq=qq(B)
    Di=la.inv(qq(D))
    e0=-float(L)
    def coefficients(mu2=0.0,e2=0.0):
        C=-mu2*np.diag(N)/2
        X={1:Aq,2:Bq-e0*np.eye(len(q)),4:qq(C)-e2*np.eye(len(q))}
        R=[-Di]
        for order in range(1,5):
            R.append(-Di@sum((X[j]@R[order-j] for j in X if j<=order),
                             np.zeros_like(Di)))
        F={}
        for order in range(2,7):
            F[order]=Apq@R[order-2]@Apq.conj().T
            if order>=3:
                F[order]+=Apq@R[order-3]@Bpq.conj().T+Bpq@R[order-3]@Apq.conj().T
            if order>=4:
                F[order]+=Bpq@R[order-4]@Bpq.conj().T
        F[2]+=pp(B)
        F[4]+=pp(C)
        return F
    F=coefficients()
    pN=N[p]
    i0=np.flatnonzero(pN==0); i4=np.flatnonzero(pN==4)
    _,v4,_=__import__('quartet').ground(F[2][np.ix_(i4,i4)])
    T=np.zeros((len(p),2),complex); T[i0[0],0]=1; T[i4,1]=v4
    U=la.null_space(T.conj().T)
    R=U@la.inv(e0*np.eye(U.shape[1])-U.conj().T@F[2]@U)@U.conj().T
    M2=T.conj().T@(F[4]+F[3]@R@F[3])@T
    mu2=float(((M2[1,1]-M2[0,0])/2).real)
    e2=float(M2[0,0].real)
    F=coefficients(mu2,e2)
    X,Y,Z,W=F[3],F[4],F[5],F[6]
    M4=T.conj().T@(W+Z@R@X+X@R@Z+Y@R@Y
                   +Y@R@X@R@X+X@R@Y@R@X+X@R@X@R@Y
                   +X@R@X@R@X@R@X-e2*X@R@R@X)@T
    return dict(blockaded_dimension=len(p), high_energy_dimension=len(q),
                leading_energy=e0, energy_epsilon4=e2,
                mu_eff_epsilon2=mu2,
                leading_off_diagonal_abs=float(abs(M2[1,0])),
                epsilon6_off_diagonal_abs=float(abs(M4[1,0])),
                predicted_gap_over_2_t_eps4=float(2*abs(M4[1,0])),
                leading_null_residual=float(la.norm((F[2]-e0*np.eye(len(p)))@T)),
                coefficient_hermiticity_error=float(max(la.norm(x-x.conj().T) for x in F.values())))


def main():
    start=time.perf_counter()
    with threadpool_limits(limits=1):
        out=dict(results=[avoided(14,2,[1/10,1/16,1/24]),
                          avoided(14,2,[1/10,1/16,1/24],True),
                          avoided(12,1,[1/12,1/20,1/32,1/48]),
                          avoided(12,1,[1/12,1/20,1/32,1/48],True)],
                 dark_fourth_order=fourth_order_dark(),
                 microscopic_dark_series=microscopic_dark_series())
    out['seconds']=time.perf_counter()-start
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
