"""Polynomial-size exact LTI optimization for even quadratic qubit bonds."""
from core import *
def coefficient(h):
    sm=np.array([[0,0],[1,0]],complex);c=[np.kron(sm,I),np.kron(-P[3],sm)];W=[a for cj in c for a in [cj,cj.conj().T]]
    ops=[a@b.conj().T for a in W for b in W]
    k=la.lstsq(np.stack([o.flatten() for o in ops],axis=1),h.flatten())[0].reshape(4,4)
    assert la.norm(sum(z*o for z,o in zip(k.flatten(),ops))-h)<1e-10
    return herm(k.T)
def positive(a):
    w,v=la.eigh(herm(a));return (v*np.maximum(w,0))@v.conj().T
def solve(h,n):
    start=time.time();D=2*n;A=np.zeros((D,D));A[:4,:4]=coefficient(h).real
    assert np.max(abs(h.imag))<1e-13
    Gv=cp.Variable((D,D),symmetric=True);J=np.kron(np.eye(n),P[1].real);idn=np.eye(D)
    c1=Gv>>0;c2=idn-Gv>>0;eq=Gv[:-2,:-2]-Gv[2:,2:]==0;ph=Gv+J@Gv.T@J==idn
    prob=cp.Problem(cp.Minimize(cp.trace(A@Gv)),[c1,c2,eq,ph])
    prob.solve(solver='CLARABEL',tol_feas=1e-9,tol_gap_abs=1e-9,tol_gap_rel=1e-9)
    solver='CLARABEL covariance'
    if prob.status!='optimal':
        prob.solve(solver='SCS',eps=1e-9,max_iters=20000)
        solver='SCS covariance'
    X=herm(eq.dual_value);Y=herm(ph.dual_value);Zp=positive(c1.dual_value);Zn=positive(c2.dual_value)
    R=A+Y+J@Y.T@J-Zp+Zn;R[:-2,:-2]+=X;R[2:,2:]-=X
    re=la.eigvalsh(herm(R));lower=float(-np.trace(Y).real-np.trace(Zn).real+np.minimum(re,0).sum())-1e-10
    raw=herm(Gv.value);Q=np.zeros_like(raw)
    for r in range(-n+1,n):
        pairs=[(i,i+r) for i in range(n) if 0<=i+r<n]
        block=sum(raw[2*i:2*i+2,2*j:2*j+2] for i,j in pairs)/len(pairs)
        for i,j in pairs:Q[2*i:2*i+2,2*j:2*j+2]=block
    Q=(Q+idn-J@Q.T@J)/2;Q=herm(Q)
    w=la.eigvalsh(Q);delta=max(0,-w[0],w[-1]-1)+1e-12;Q=(Q+delta*idn)/(1+2*delta)
    upper=float(np.trace(A@Q).real)
    rec=dict(n=n,value=float(prob.value),lower=lower,upper=upper,status=prob.status,solver=solver,seconds=time.time()-start,dual_gap=upper-lower,stationarity_residual=float(la.norm(R)),lti_residual=float(la.norm(raw[:-2,:-2]-raw[2:,2:])),particle_hole_residual=float(la.norm(raw+J@raw.T@J-idn)),min_covariance=float(la.eigvalsh(raw)[0]),repair_mix=float(2*delta/(1+2*delta)))
    return rec,Q
if __name__=='__main__':
    out=[]
    for path in sorted(Path('results').glob('*.json')):
        rec=json.loads(path.read_text());name=rec['name']
        if name not in ['XX','XX_field'] and not name.startswith('TFIM'):continue
        h=unpack(rec['h']);rec.setdefault('covariance_lti',{})
        for n in range(2,11):
            if str(n) in rec['covariance_lti'] and rec['covariance_lti'][str(n)].get('version',0)==2:continue
            row,Q=solve(h,n);row['real_formulation']=True;row['version']=2;rec['covariance_lti'][str(n)]=row
            if n<=8 and not rec['lti'][str(n)]['solver'].endswith('covariance'):
                rec.setdefault('full_spin_lti',{})[str(n)]=rec['lti'][str(n)]
                row['full_spin_value_difference']=abs(row['value']-rec['lti'][str(n)]['value'])
            rec['lti'][str(n)]=row
            write_json(path,rec);print(name,n,row,flush=True)
