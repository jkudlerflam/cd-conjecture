"""C24 independent checks. No imports from the candidate implementation."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from itertools import product, combinations, combinations_with_replacement
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.linalg import block_diag


def stage(state, edges, upper, right=False):
    n=len(edges)
    for a in range(upper):
        yield (a,('g',a),state), 1/np.sqrt(2*upper)
    for e,(a,b) in enumerate(edges):
        occupied=(state>>e)&1
        phase=(-1)**((state & ((1<<e)-1)).bit_count())
        if right: phase*=(-1)**(state.bit_count()-occupied)
        yield (b,('h',e) if occupied else ('g',a),state^(1<<e)),phase/np.sqrt(2*n)


def ports(upper,lower,edges):
    n=len(edges);d=1<<n;q=upper+lower
    rows={}; entries=[]
    for x in range(d):
        for (l,fl,y),vl in stage(x,edges,upper):
            for (r,fr,z),vr in stage(y,edges,upper,True):
                key=(l,fl,r,fr,z)
                if key not in rows:rows[key]=len(rows)
                entries.append((rows[key],x,vl*vr))
    rr,cc,vv=zip(*entries); W=coo_matrix((vv,(rr,cc)),shape=(len(rows),d)).tocsr()
    eye=np.eye(q,dtype=int);const=4*np.ones(q,dtype=int)
    weights=np.array([eye[a]-eye[b] for a,b in edges])
    def flag(f):
        kind,i=f
        if kind=='g':return 2*np.ones(q,dtype=int)-eye[i]
        a,b=edges[i];return 2*np.ones(q,dtype=int)+eye[a]-2*eye[b]
    def logical(s):return sum((weights[e] for e in range(n) if (s>>e)&1),np.zeros(q,dtype=int))
    for row,x,_ in entries:
        pass
    keys=list(rows)
    for row,x,_ in entries:
        l,fl,r,fr,z=keys[row]
        assert np.array_equal(eye[l]+flag(fl)+eye[r]+flag(fr)+logical(z),const+logical(x))
    errors=[np.max(np.abs((W.T@W).toarray()-np.eye(d)))];eps=1/(2*np.sqrt(upper*n))
    for e,(a,b) in enumerate(edges):
        f=np.zeros((d,d))
        for x in range(d):
            if not ((x>>e)&1): f[x^(1<<e),x]=(-1)**((x&((1<<e)-1)).bit_count())
        for side in [0,2]:
            ar=[];ac=[]
            for key,j in rows.items():
                if key[side]!=b:continue
                k=list(key);k[side]=a;k=tuple(k)
                if k in rows:ar.append(rows[k]);ac.append(j)
            F=coo_matrix((np.ones(len(ar)),(ar,ac)),shape=(len(rows),len(rows))).tocsr()
            target=eps*f
            if side==2:target=target@(np.diag([(-1)**x.bit_count() for x in range(d)]))*(1-1/(2*n))
            errors.append(np.max(np.abs((W.T@F@W).toarray()-target)))
    return dict(N=n,q=q,joint_nonzero_rows=len(rows),fock_dimension=d,maximum_entry_error=float(max(errors)),charge_checks=len(entries))


def plane(p):
    # Representatives with first nonzero coordinate 1, generated without normalization.
    v=[(1,a,b) for a in range(p) for b in range(p)]+[(0,1,a) for a in range(p)]+[(0,0,1)]
    v=np.array(v,dtype=int);A=(v@v.T)%p==0;s=len(v)
    edges=list(zip(*np.nonzero(A))); W=np.zeros((2*s,len(edges)),dtype=int)
    for e,(a,b) in enumerate(edges):W[a,e]=1;W[s+b,e]=-1
    assert len(edges)==(p+1)*s and np.all(A.sum(axis=0)==p+1)
    assert np.array_equal(A.astype(int)@A.astype(int).T,p*np.eye(s,dtype=int)+np.ones((s,s),dtype=int))
    # Any <=4-field neutral conversion is equality of two multisets of <=2 edges.
    sums={}
    for k in [1,2]:
        sums={}
        for ids in combinations_with_replacement(range(len(edges)),k):
            charge=tuple(W[:,ids].sum(axis=1))
            assert charge not in sums or sums[charge]==ids
            sums[charge]=ids
    # Enumerate six-cycles by noncollinear point triples and all chiral choices.
    edge_id={tuple(e):i for i,e in enumerate(edges)};cycles=0;zero=0;example=None
    if p<=3:
        n=len(edges);den=12*(n+1);k_num=3*(n+1)+np.arange(1,n+1)
        for a,b,c in combinations(range(s),3):
            ab=int(np.flatnonzero(A[a]&A[b])[0]);bc=int(np.flatnonzero(A[b]&A[c])[0]);ca=int(np.flatnonzero(A[c]&A[a])[0])
            if len({ab,bc,ca})<3:continue
            ids=[edge_id[e] for e in [(a,ab),(b,ab),(b,bc),(c,bc),(c,ca),(a,ca)]]
            z=np.array([1,-1,1,-1,1,-1]);assert np.all(W[:,ids]@z==0);cycles+=1
            for m in product([-1,1],repeat=6):
                if sum(int(k_num[i])*mi for i,mi in zip(ids,m))%(2*den)==0:
                    zero+=1
                    if example is None:example=dict(edge_labels=[i+1 for i in ids],z=z.tolist(),m=list(m))
        assert cycles==s*(s-1)*(s-2)//6-s*(p+1)*p*(p-1)//6
    return dict(p=p,q=2*s,N=len(edges),one_two_edge_multisets=len(sums),six_cycles=cycles if p<=3 else None,resonant_cycle_chiral_assignments=zero if p<=3 else None,example=example)


def acceptance_tests():
    rng=np.random.default_rng(24027);out=[]
    for d,r,g in [(8,4,2),(15,7,3),(25,10,4)]:
        for rep in range(20):
            gamma=.1+rng.random();H=np.diag([0.]*r+[gamma]*(d-r))
            X=rng.normal(size=(d,d-g))+1j*rng.normal(size=(d,d-g))
            A=block_diag(np.zeros((g,g)),X.conj().T@X)
            a=np.linalg.eigvalsh(A[g:r,g:r])[0];B=np.linalg.eigvalsh(A)[-1]
            bound=a*gamma/(a+B+gamma);ev=np.linalg.eigvalsh(H+A)
            assert np.count_nonzero(np.abs(ev)<1e-10)==g and ev[g]>=bound-1e-10
            out.append(float(ev[g]/bound))
    return dict(cases=len(out),dimensions=[8,15,25],minimum_actual_gap_over_bound=min(out))


def gaussian():
    rng=np.random.default_rng(24);out=[]
    for n,b in [(1,4),(6,20),(21,100)]:
        nu=.5;a=b+nu;f=.25+np.arange(1,n+1)/(12*(n+1));v=np.r_[1,np.tan(np.pi*f)]
        T=np.eye(n+1);T[0,0]=a;T[1:,0]=-f;Ti=np.linalg.inv(T)
        A=T.T@np.diag(v)@T;B=Ti@np.diag(v)@Ti.T;K=Ti@Ti.T
        err=np.linalg.norm(K@A@K-B);assert err<1e-10
        full=block_diag(A,B);eig=np.linalg.eigvalsh(full)
        assert eig[0]>0
        # Exact diagonal reference normalized perturbation, generally includes density-current blocks.
        R=rng.normal(size=(2*n+2,2*n+2));R=(R+R.T)/2;R*=.01/np.linalg.norm(R,2)
        vals,U=np.linalg.eigh(full);root=(U*np.sqrt(vals))@U.T
        pert=root@(np.eye(len(full))+R)@root
        e,Q=np.linalg.eigh(pert);sqrt=(Q*np.sqrt(e))@Q.T;inv=(Q*(1/np.sqrt(e)))@Q.T
        J=np.block([[np.zeros_like(A),np.eye(n+1)],[np.eye(n+1),np.zeros_like(A)]])
        c,O=np.linalg.eigh(sqrt@J@sqrt);metric=inv@(O*np.abs(c))@O.T@inv
        ref=block_diag(K,np.linalg.inv(K));e,Q=np.linalg.eigh(ref);refinv=(Q*(1/np.sqrt(e)))@Q.T
        kappa=np.linalg.eigvalsh(refinv@metric@refinv)[0]
        assert kappa>2/3
        out.append(dict(N=n,b=b,reference_stiffness_min=float(eig[0]),KAK_error=float(err),normalized_scaling_metric_min=float(kappa),six_cycle_dimension_lower_bound=float(3*kappa)))
    return out


def history_toy():
    # Literal history parent with a noncommuting acceptance pin, good dimension four.
    # Input bits: two logical bits, one invalid-format bit, one trusted scratch bit.
    d=16;times=9;Pinit=np.diag([int((x>>3)&1) for x in range(d)])
    U=np.eye(d);states=[U.copy()]
    gates=[]
    cnot=np.zeros((d,d))
    for x in range(d):cnot[x^(((x>>2)&1)<<3),x]=1
    # Compute FAIL into trusted bit, keep it, then erase, repeated cycle.
    gates=[cnot,np.eye(d),np.eye(d),cnot]*2
    for g in gates:U=g@U;states.append(U.copy())
    H=np.zeros((times*d,times*d))
    for t,g in enumerate(gates):
        l=slice(t*d,(t+1)*d);r=slice((t+1)*d,(t+2)*d)
        H[l,l]+=.5*np.eye(d);H[r,r]+=.5*np.eye(d);H[r,l]-=.5*g;H[l,r]-=.5*g.T
    H[:d,:d]+=Pinit
    # All-time FAIL; unlike terminal only, catches bad inputs that reset.
    fail=np.diag([int((x>>3)&1) for x in range(d)])
    alltime=H+np.kron(np.eye(times),fail)
    terminal=H.copy();terminal[-d:,-d:]+=fail
    ea=np.linalg.eigvalsh(alltime);et=np.linalg.eigvalsh(terminal)
    assert sum(abs(ea)<1e-10)==4 and sum(abs(et)<1e-10)==8
    return dict(total_dimension=times*d,history_vertices=times,all_time_kernel=4,terminal_only_kernel=8,all_time_gap=float(ea[4]))


def main():
    graphs=[(1,1,[(0,1)]),(2,2,[(0,2),(1,2),(1,3)]),(2,2,[(0,2),(0,3),(1,2),(1,3)]),(3,3,[(a,3+b) for a in range(3) for b in range(3) if a!=b])]
    out={'ports':[ports(*g) for g in graphs], 'planes':[plane(p) for p in [2,3,5,7]], 'acceptance':acceptance_tests(),'history_toy':history_toy(),'gaussian':gaussian()}
    Path('independent_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
