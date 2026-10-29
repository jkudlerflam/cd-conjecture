"""Independent integer verifier, plus small floating consistency diagnostics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh


def ints(a):
    return np.array([int(x) for x in np.asarray(a).ravel()],dtype=object).reshape(np.asarray(a).shape)


def verify(path):
    c=json.loads(path.read_text()); d=c['d']; k=c['k']; D=d**k; n=d**(k-1)
    h=np.array(c['bond']); assert np.array_equal(h*16,np.rint(h*16))
    h16=ints(h*16); eye=lambda n:np.eye(n,dtype=object)
    K=np.kron(h16,eye(d**(k-2)))*2**28
    coeff=iter(c['X_coeff_integer'])
    for a in range(n):
        for b in range(a,n):
            if a==b==n-1: continue
            X=np.zeros((n,n),dtype=object); X[a,b]=X[b,a]=next(coeff)
            K+=(np.kron(X,eye(d))-np.kron(eye(d),X))*2**16
    assert next(coeff,None) is None
    if c['ground_constraints']:
        H=np.kron(h16,eye(d))+np.kron(eye(d),h16)
        ops=[]
        for a in range(d):
            for b in range(d):
                o=np.zeros((d,d),dtype=object); o[a,b]=1
                ops.append(np.kron(eye(d),np.kron(o,eye(d))))
        yi=iter(c['stationarity_integer'])
        for a in range(d):
            for b in range(a+1,d):
                o=ops[d*a+b]-ops[d*b+a]
                K+=next(yi)*np.kron(H@o-o@H,eye(d**(k-3)))*2**12
        B=ints(c['B_integer']); W=B@B.T
        for a,oa in enumerate(ops):
            for b,ob in enumerate(ops):
                Gnum=oa.T@(H@ob-ob@H)+(oa.T@H-H@oa.T)@ob
                K-=W[a,b]*np.kron(Gnum,eye(d**(k-3)))*8
    assert np.array_equal(K,ints(c['K_integer']))
    assert np.array_equal(K,K.T)
    Q=ints(c['Q_integer']); ei=c['e_trial_integer']
    R=K*2**16-ei*2**16*eye(D)-Q@Q.T
    delta=max(sum(abs(x) for x in row) for row in R)
    assert delta==c['residual_row_bound_integer']
    assert ei*2**16-delta==c['bound_numerator']
    assert c['bound_denominator']==2**48
    # Symmetric R has operator norm at most its largest absolute row sum.
    return c['bound_numerator']/c['bound_denominator']


def reduced(psi,d,L,k):
    a=psi.reshape(d**k,d**(L-k))
    return a@a.T


def diagnostics():
    from bootstrap import bonds,ring,setup
    worst_stationarity=0.; worst_marginal=0.; min_g=0.; min_slack=1e9
    for name,d,h in bonds():
        L=4
        H=ring(h,d,L); en,U=eigh(H); p=U[:,0]
        # Average translates to handle degenerate ground spaces without symmetry assumptions.
        rho=np.zeros((d**3,d**3))
        tensor=p.reshape((d,)*L)
        for t in range(L):
            pp=np.transpose(tensor,tuple(range(t,L))+tuple(range(t))).reshape(-1)
            rho+=reduced(pp,d,L,3)/L
        four=rho.reshape(d,d,d,d,d,d)
        left=np.trace(four,axis1=2,axis2=5).reshape(d*d,d*d)
        right=np.trace(four,axis1=0,axis2=3).reshape(d*d,d*d)
        worst_marginal=max(worst_marginal,float(np.max(abs(left-right))))
        C,DX,S,G=setup(h,d,3,True)
        worst_stationarity=max(worst_stationarity,float(np.max(abs(np.einsum('ij,aji->a',rho,S)))))
        gram=np.einsum('ij,abji->ab',rho,G)
        min_g=min(min_g,float(np.linalg.eigvalsh(gram)[0]))
        for path in Path('certificates').glob(name+'_*.json'):
            lb=verify(path)
            if json.loads(path.read_text())['k'] <= L:
                min_slack=min(min_slack,en[0]/L-lb)
    assert worst_marginal<1e-10 and worst_stationarity<1e-10 and min_g>-1e-10 and min_slack>-1e-10
    return dict(max_marginal_residual=worst_marginal,max_stationarity_residual=worst_stationarity,
                min_ground_gram_eigenvalue=min_g,min_ring_minus_bound_per_site=min_slack)

if __name__=='__main__':
    paths=sorted(Path('certificates').glob('*.json'))
    for path in paths: print('EXACT CERTIFICATE PASS',path,verify(path))
    print('DIAGNOSTICS',json.dumps(diagnostics(),sort_keys=True))
