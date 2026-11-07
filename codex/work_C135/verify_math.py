"""Independent small spin-matrix checks of the analytic arguments."""
from core import *
from certificates import embed
def jw(n):
    sm=np.array([[0,0],[1,0]],complex);out=[]
    for j in range(n):
        a=np.ones((1,1),complex)
        for k in range(n):a=np.kron(a,-P[3] if k<j else (sm if k==j else I))
        out.append(a)
    return out
out={}
h=xxz(0,.5);X=-xxz(0)/2-np.kron(P[3],I)/12-np.kron(I,P[3])/24
K=np.kron(h,I)+np.kron(X,I)-np.kron(I,X)
out['exact_n3_eigenvalues']=la.eigvalsh(K).tolist()
assert abs(out['exact_n3_eigenvalues'][0]+5/12)<1e-13
# Computational bit 0 is up. One up spin over the all-down reference.
w=np.zeros(8);w[[3,5,6]]=[1,-1,1];w/=np.sqrt(3);rho=np.outer(w,w)
out['w_lti_residual']=float(la.norm(left(rho)-right(rho)))
rot=w.reshape(2,2,2).transpose(2,0,1).reshape(8);rrot=np.outer(rot,rot)
out['cyclic_window_lti_residual']=float(la.norm(left(rrot)-right(rrot)))
seam=ring_matrix(h,3).toarray()-embed(h,0,3)-embed(h,1,3)
g1=np.kron(P[3],np.eye(4));tw=embed(h,0,3)+embed(h,1,3)+g1@seam@g1
out['twisted_n3_eigenvalues']=la.eigvalsh(tw).tolist();out['w_twisted_residual']=float(la.norm(tw@w+1.25*w))
assert out['w_twisted_residual']<1e-13
ints=[3187,-2870,2007,-756];C=np.array([[ints[abs(i-j)] for j in range(4)] for i in range(4)])/10000
ev,U=la.eigh(C);B0=(U*np.log(ev/(1-ev)))@U.T;c=jw(4)
Q=sum(B0[i,j]*c[i].conj().T@c[j] for i in range(4) for j in range(4));r=la.expm(Q);r/=np.trace(r)
out['rational_n4']=dict(lti_residual=float(la.norm(left(r)-right(r))),energy=float(np.trace(np.kron(h,np.eye(4))@r).real),covariance_residual=float(la.norm(np.array([[np.trace(c[i].conj().T@c[j]@r) for j in range(4)] for i in range(4)])-C)))
assert abs(out['rational_n4']['energy']+7553/20000)<1e-12
assert out['rational_n4']['lti_residual']<1e-12
out['tail_lower']=float(-np.sqrt(3)/(2*np.pi)-1/12-np.pi**2/768)
assert out['tail_lower']>-.375
out['twisted_rings']={}
for L in range(4,8):
    vals=[-np.sin(np.pi*m/L)/(L*np.sin(np.pi/L))+.5*(m/L-.5) for m in range(L+1)]
    out['twisted_rings'][L]=dict(minimum=min(vals),occupation=int(np.argmin(vals)))
    assert min(vals)>=-.375-1e-13
# Paired generalized covariance, including particle-hole compatibility.
checks=[]
for name in ['XX','TFIM_1','Heisenberg']:
    f=np.load(f'matrices/{name}_4.npz');r=f['rho'];n=4;c=jw(n);W=[a for cj in c for a in [cj,cj.conj().T]]
    G0=np.array([[np.trace(a@b.conj().T@r) for b in W] for a in W]);J=np.kron(np.eye(n),P[1]);v=np.sqrt(2/(n+1))*np.sin(np.pi*np.arange(1,n+1)/(n+1));s=np.cos(np.pi/(n+1))
    def symbol(k):
        V=np.kron((v*np.exp(1j*k*np.arange(n)))[:,None],I);return V.conj().T@G0@V
    grid=2*np.pi*np.arange(64)/64;gs=np.array([symbol(k) for k in grid]);mineig=min(la.eigvalsh(g)[0] for g in gs);maxeig=max(la.eigvalsh(g)[-1] for g in gs)
    ph=max(la.norm(symbol(k)+P[1]@symbol(-k).T@P[1]-I) for k in grid)
    lag0=gs.mean(0);lag1=np.einsum('k,kij->ij',np.exp(-1j*grid),gs)/len(grid)
    rec=dict(name=name,min_symbol=float(mineig),max_symbol=float(maxeig),particle_hole_residual=float(ph),onsite_residual=float(la.norm(lag0-G0[:2,:2])),neighbor_residual=float(la.norm(lag1-s*G0[:2,2:4])))
    assert mineig>-1e-10 and maxeig<1+1e-10 and ph<1e-10 and rec['neighbor_residual']<1e-10
    checks.append(rec)
out['covariance_smoothing']=checks;write_json('math_checks.json',out);print(json.dumps(out,indent=2))
