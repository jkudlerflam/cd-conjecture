"""C132 support-reduced CPTP recovery SDP. All indices use input x output order."""
import numpy as np
import scipy.linalg as la
import cvxpy as cp
import fr_lib as F

def support(r,cut=1e-11):
    w,u=la.eigh((r+r.conj().T)/2)
    mask=w>cut*max(w.max(),1e-300)
    return w[mask],u[:,mask]

def solve(h,r,l,d,eps=2e-7,max_iters=30000,return_matrix=False,solver="SCS"):
    rb=F.ptrace_first(r,d,l)
    wb,u=support(rb); wo,v=support(r)
    b,q=len(wb),len(wo)
    inp=np.kron(np.eye(d),u)
    x=(inp.conj().T@r@inp).reshape(d,b,d,b)
    hh=h.reshape(d,d,d,d)
    C=np.zeros((b*q,b*q),complex)
    for a in range(d):
        for ap in range(d):
            z=v.conj().T@np.kron(hh[ap,:,a,:],np.eye(d**(l-1)))@v
            C+=np.kron(x[a,:,ap,:].T,z)
    C=(C+C.conj().T)/2
    real=np.max(abs(C.imag))<1e-12
    J=cp.Variable((b*q,b*q),symmetric=True) if real else cp.Variable((b*q,b*q),hermitian=True)
    tp=cp.bmat([[cp.reshape(cp.trace(J[i*q:(i+1)*q,j*q:(j+1)*q]),(1,1),order='C') for j in range(b)] for i in range(b)])
    rec=sum(wb[i]*J[i*q:(i+1)*q,i*q:(i+1)*q] for i in range(b))
    cons=[J>>0,tp==np.eye(b),rec==np.diag(wo)]
    prob=cp.Problem(cp.Minimize(cp.real(cp.trace(C@J))),cons)
    if solver == 'SCS':
        val=prob.solve(solver='SCS',eps=eps,max_iters=max_iters,verbose=False)
    else:
        val=prob.solve(solver=solver,tol_gap_abs=eps,tol_feas=eps,tol_gap_rel=eps,verbose=False)
    j=J.value
    if j is None: raise RuntimeError(prob.status)
    tpn=np.array([[np.trace(j[i*q:(i+1)*q,k*q:(k+1)*q]) for k in range(b)] for i in range(b)])
    rn=sum(wb[i]*j[i*q:(i+1)*q,i*q:(i+1)*q] for i in range(b))
    eig=la.eigvalsh(j)
    # CVXPY equality multipliers use the opposite sign to our maximizing dual.
    y=-cons[1].dual_value; z=-cons[2].dual_value
    slack=C-np.kron(y,np.eye(q))-np.kron(np.diag(wb),z)
    slack=(slack+slack.conj().T)/2
    mine=la.eigvalsh(slack)[0]
    dual=float((np.trace(y)+np.trace(z@np.diag(wo))).real)
    # Shifting Y down by the negative slack eigenvalue makes a numerical dual bound.
    corrected_dual=dual+b*min(0,float(mine))
    # Reconstruct Petz Choi in the same coordinates, a feasible objective check.
    jp=np.zeros_like(j,dtype=complex)
    root=F.psd_power(r,.5); inv=F.psd_power(rb,-.5)
    for c in range(d):
        K=v.conj().T@root@np.kron(inv@u,np.eye(d)[:,c:c+1])
        kv=K.T.reshape(-1)
        jp+=np.outer(kv,kv.conj())
    out=dict(energy=float(val),status=prob.status,input_rank=b,output_rank=q,choi_size=b*q,input_discarded_mass=float(np.trace(rb).real-wb.sum()),output_discarded_mass=float(np.trace(r).real-wo.sum()),
        dual=dual,dual_min_eig=float(mine),corrected_dual=corrected_dual,gap=float(val-corrected_dual),tp_res=float(np.max(abs(tpn-np.eye(b)))),rec_res=float(np.max(abs(rn-np.diag(wo)))),
        min_eig=float(eig[0]),choi_rank=int(np.sum(eig>1e-5)),petz_energy=float(np.trace(C@jp).real),
        distance_petz=float(la.norm(j-jp)),iterations=prob.solver_stats.num_iters)
    if return_matrix: out.update(J=j,C=C,U=u,V=v,wb=wb,wo=wo)
    return out

def gauge_scale(h):
    d=F.bond_dim(h); o=cp.Variable((d,d),hermitian=True); a=cp.Variable(); m=cp.Variable()
    z=h+cp.kron(o,np.eye(d))-cp.kron(np.eye(d),o)+a*np.eye(d*d)
    p=cp.Problem(cp.Minimize(m),[z<<m*np.eye(d*d),z>>-m*np.eye(d*d),cp.trace(o)==0])
    p.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
    return float(m.value)

if __name__=='__main__':
    h=F.heis_bond(); L=6; _,v,e=F.ground_states(h,L)
    for l in range(2,L):
        r=F.ti_marginal(v,2,L,l); z=solve(h,r,l,2)
        pet=F.petz_defect(h,v,L,l)
        assert abs(pet-z['petz_energy'])<1e-9
        print(l,z,'scaled',L*L*(z['energy']-e/L)/.5,flush=True)
