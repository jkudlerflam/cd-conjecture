"""Inspect optimal numerical duals and independently verify the analytic XX dual."""
from core import *
def embed(h,j,n):return np.kron(np.kron(np.eye(2**j),h),np.eye(2**(n-j-2)))
def inspect(x,h,n):
    m=n-1;D=2**m;x=x-np.trace(x)*np.eye(D)/D
    basis=[embed(h,j,m) for j in range(m-1)]
    A=np.stack([b.flatten() for b in basis],axis=1)
    w=la.lstsq(A,x.flatten())[0];res=x-sum(wj*b for wj,b in zip(w,basis))
    c=pauli_transform(x,m).flatten();weights={};spans={};largest=[]
    for i,z in enumerate(c):
        digs=[];v=i
        for _ in range(m):digs.append(v%4);v//=4
        digs=digs[::-1];nz=[j for j,a in enumerate(digs) if a];weight=len(nz);span=nz[-1]-nz[0]+1 if nz else 0
        weights[weight]=weights.get(weight,0)+float(abs(z)**2);spans[span]=spans.get(span,0)+float(abs(z)**2)
        if abs(z)>1e-6:largest.append((''.join('IXYZ'[a] for a in digs),float(z.real)))
    total=sum(weights.values())
    return dict(bond_coefficients=w.real.tolist(),relative_bond_fit_residual=float(la.norm(res)/la.norm(x)),pauli_weight_fractions={k:v/total for k,v in weights.items()},pauli_span_fractions={k:v/total for k,v in spans.items()},largest=sorted(largest,key=lambda v:-abs(v[1]))[:20])
def verify():
    out=[]
    for n in range(2,33):
        v=np.sqrt(2/(n+1))*np.sin(np.pi*np.arange(1,n+1)/(n+1));s=np.cos(np.pi/(n+1));a=v[:-1]*v[1:]/s;x=np.cumsum(a)[:-1]-1
        one=np.diag(a/2,1)+np.diag(a/2,-1);ev=la.eigvalsh(one);E=float(ev[ev<0].sum());t=-1/(np.pi*s)
        residual=None
        if n<=8:
            X=sum((x[j]*embed(xxz(0),j,n-1) for j in range(n-2)),np.zeros((2**(n-1),)*2,complex))
            K=np.kron(xxz(0),np.eye(2**(n-2)))+np.kron(X,I)-np.kron(I,X)
            expected=sum((a[j]*embed(xxz(0),j,n) for j in range(n-1)),np.zeros((2**n,)*2,complex))
            residual=float(la.norm(K-expected));assert residual<1e-12
            assert abs(la.eigvalsh(K,subset_by_index=[0,0])[0]-E)<1e-12
        assert E>=t-1e-12
        assert n*n*(1/(np.pi*s)-1/np.pi)/.5<=np.pi+1e-12
        out.append(dict(n=n,weighted_ground=E,analytic_lower=t,identity_residual=residual))
    write_json('analytic_checks.json',out)
    structures=[]
    for name in ['XX','Heisenberg','TFIM_1']:
        for n in range(4,9):
            path=Path(f'matrices/{name}_{n}.npz')
            if path.exists():
                f=np.load(path);structures.append(dict(name=name,n=n,**inspect(f['X'],f['h'],n)))
    write_json('dual_structure.json',structures)
    print('Analytic XX certificate checks passed for n=2..32; telescoping matrix checks n=2..8.')
if __name__=='__main__':verify()
