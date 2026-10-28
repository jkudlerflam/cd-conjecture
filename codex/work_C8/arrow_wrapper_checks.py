"""Complete open-cell wrapper with typed vacancies and endpoint rewards."""
import itertools
import json
import numpy as np
from scipy import sparse as sp
from length_filter_checks import filter_data, plus, tensors
from vacancy_parent import charge_hamiltonian


def arrow_data(q):
    starts = [0,1,2]+[0]*q+[1]*q+[1]*q
    ends = [0,1,2]+[1]*q+[1]*q+[2]*q
    return np.array(starts), np.array(ends)


def legal_words(q, n):
    starts, ends = arrow_data(q)
    words = [(a,) for a in range(3+3*q)]
    for _ in range(n-1):
        words = [s+(a,) for s in words for a in range(3+3*q)
                 if ends[s[-1]] == starts[a]]
    return words


def local_bond(q, h0, a, b, U, w, mu, K):
    d = 3+3*q
    start, end = arrow_data(q)
    bad = np.array([end[x] != start[y] for x in range(d) for y in range(d)])
    charge = np.diag([U*((x<3) != (y<3)) for x in range(d) for y in range(d)])
    for x in range(q):
        for pair1, pair2 in [((0,3+x),(3+x,1)),
                             ((1,3+q+x),(3+q+x,1)),
                             ((1,3+2*q+x),(3+2*q+x,2))]:
            i, j = pair1[0]*d+pair1[1], pair2[0]*d+pair2[1]
            charge[i,j] = charge[j,i] = -w
    h = charge.astype(complex)+K*np.diag(bad)
    for ta, tb in ((0,1),(0,2),(1,1),(1,2)):
        ids = [(3+ta*q+x)*d+3+tb*q+y for x in range(q) for y in range(q)]
        h[np.ix_(ids, ids)] += h0
    aa, bb = np.zeros((d,d),complex), np.zeros((d,d),complex)
    aa[3:3+q,3:3+q] = a-mu*np.eye(q)
    bb[3+2*q:,3+2*q:] = b-mu*np.eye(q)
    h += np.kron(aa,np.eye(d))+np.kron(np.eye(d),bb)
    return h, charge, np.diag(bad.astype(float))


def restricted_h(words, local, d):
    index = {s:i for i,s in enumerate(words)}
    out = np.zeros((len(words),len(words)),complex)
    transitions = [[(a//d,a%d,local[a,b]) for a in np.flatnonzero(local[:,b])]
                   for b in range(d*d)]
    for col, s in enumerate(words):
        for j in range(len(s)-1):
            for x,y,z in transitions[s[j]*d+s[j+1]]:
                target=s[:j]+(x,y)+s[j+2:]
                assert target in index  # legal subspace must be invariant.
                out[index[target],col] += z
    return out


def target_frame(k, p, phase, words, U, w):
    q=2
    holes,hc,chi,ec=charge_hamiltonian(k,p,U,w)
    index={s:i for i,s in enumerate(words)}
    out=np.zeros(len(words),complex)
    payload=tensors([plus(j*phase) for j in range(k)])
    for hs, amp in zip(holes, chi):
        occ=[j for j in range(k+p) if j not in hs]
        for bs in itertools.product(range(q),repeat=k):
            occupied=iter([3+tag*q+x for tag,x in zip([0]+[1]*(k-2)+[2],bs)])
            s=tuple((0 if j<occ[0] else 2 if j>occ[-1] else 1)
                    if j in hs else next(occupied) for j in range(k+p))
            out[index[s]]=amp*payload[np.ravel_multi_index(bs,(q,)*k)]
    return out,ec[0],holes,chi


def main():
    k,q=3,2
    U,w,mu=1.,.1,2.
    K=2*mu+2*U
    phase,hf,af,bf,epsilon=filter_data(k)
    scale=2*U/epsilon
    h,charge,bad=local_bond(q,scale*hf,scale*af,scale*bf,U,w,mu,K)
    d=3+3*q
    h12=sp.kron(sp.csr_matrix(h),sp.eye(d),format='csr')
    h23=sp.kron(sp.eye(d),sp.csr_matrix(h),format='csr')
    b12=sp.kron(sp.csr_matrix(bad),sp.eye(d),format='csr')
    b23=sp.kron(sp.eye(d),sp.csr_matrix(bad),format='csr')
    comm=max(float(sp.linalg.norm(b12@h23-h23@b12)),
             float(sp.linalg.norm(b23@h12-h12@b23)),
             float(np.linalg.norm(bad@h-h@bad)))
    assert comm<1e-12 and np.linalg.eigvalsh(charge)[0]>=-1e-12
    # Check the combinatorial bound on every word at length 5.
    start,end=arrow_data(q)
    excess=-10
    for s in itertools.product(range(d),repeat=5):
        violations=sum(end[x]!=start[y] for x,y in zip(s[:-1],s[1:]))
        endpoints=sum((3<=x<3+q) or (3+2*q<=x<d) for x in s)
        excess=max(excess,endpoints-2*(violations+1))
    assert excess<=0
    rows=[]; saved=[]
    for p in (1,2):
        n=k+p
        words=legal_words(q,n)
        H=restricted_h(words,h,d)
        frame,ec,holes,chi=target_frame(k,p,phase,words,U,w)
        ev=np.linalg.eigvalsh(H)
        energy=ec-2*mu
        residual=np.linalg.norm(H@frame-energy*frame)
        assert residual<1e-10 and abs(ev[0]-energy)<1e-10
        assert ev[1]-ev[0]>1e-8
        rows.append(dict(length=n, legal_dimension=len(words),
                         full_dimension=d**n, ground_energy=float(ev[0]),
                         predicted_ground_energy=float(energy),
                         code_residual=float(residual), measured_gap=float(ev[1]-ev[0]),
                         guaranteed_illegal_sector_separation=K-2*mu-U))
        saved.append((words,frame,holes,chi))
    before={s+t:a*b for s,a in zip(saved[0][0],saved[0][1]) if a!=0
            for t,b in zip(saved[1][0],saved[1][1]) if b!=0}
    after={s+t:a*b for s,a in zip(saved[1][0],saved[1][1]) if a!=0
           for t,b in zip(saved[0][0],saved[0][1]) if b!=0}
    plain=sum(np.conj(after.get(s,0))*v for s,v in before.items())
    dressed=0j
    for s,v in before.items():
        if s[k+1]==0:
            changed=s[:k+1]+(2,)+s[k+2:]
            dressed+=np.conj(after.get(changed,0))*v
    lookup={s:i for i,s in enumerate(saved[1][2])}
    factor=sum(saved[0][3][i]*saved[1][3][lookup[(s[0],k+1)]]
               for i,s in enumerate(saved[0][2]))
    beta=factor**2
    assert abs(plain)<1e-12 and abs(dressed-beta)<1e-12
    result=dict(k=k,q=q,local_dimension=d,U=U,w=w,mu=mu,K=K,
                filter_scale=float(scale),max_violation_commutator=comm,
                endpoint_count_excess=int(excess),cells=rows,
                plain_transport_magnitude=float(abs(plain)),beta=float(beta),
                converted_transport_error=float(abs(dressed-beta)))
    with open('arrow_wrapper_checks.json','w') as f:
        json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
