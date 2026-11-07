import numpy as np
import scipy.linalg as la
import fr_lib as F
from optimal_recovery import solve

def symstate(l):
    inds=np.arange(2**l);n=np.array([int(i).bit_count() for i in inds]);v=np.zeros((2**l,l+1))
    for j in range(l+1):
        ids=n==j;v[ids,j]=1/np.sqrt(np.sum(ids))
    return v@v.T/(l+1)

for l in [2,3,4,5,6]:
    r=symstate(l); h=-F.heis_bond()
    pet=np.trace(h@F.petz_periodize(r,2,l,steps=1)[0]).real
    theory=-.25+1/(2*l*(l+1))
    assert abs(pet-theory)<1e-12,(l,pet,theory)
    if l<=4:
        z=solve(h,r,l,2,eps=1e-9,solver='CLARABEL')
        assert abs(z['energy']-theory)<1e-7,z
    print('PASS FM analytic optimum l=',l,flush=True)

# A non-Petz channel, a convex mixture with the replacement channel.
h=F.heis_bond();L=6;l=3;d=2
_,v,e=F.ground_states(h,L);r=F.ti_marginal(v,d,L,l);rb=F.ptrace_first(r,d,l)
root=F.psd_power(r,.5);inv=F.psd_power(rb,-.5)
def R(x):return .37*root@np.kron(inv@x@inv,np.eye(d))@root+.63*np.trace(x)*r
assert la.norm(R(rb)-r)<1e-12
outs=[];tau=r.copy()
for _ in range(3):
    blocks=tau.reshape(d,4,d,4)
    ext=np.block([[R(blocks[a,:,ap,:]) for ap in range(d)] for a in range(d)])
    rr=ext.reshape(4,4,4,4);outs.append(np.einsum('aibi->ab',rr))
    tau=F.ptrace_first(ext,d,l+1)
assert max(la.norm(tau-r),la.norm(outs[0]-outs[1]),la.norm(outs[0]-outs[2]))<1e-12
assert la.norm(F.ptrace_first(outs[0],d,2)-F.ptrace_last(outs[0],d,2))<1e-12
print('PASS arbitrary CPTP stabilization and equal one-site marginals',flush=True)

# Direct Choi contraction against a genuinely complex bond.
h=F.c106_bond(1,.03,2);L=4;l=2;d=4
_,v,e=F.ground_states(h,L,k=8);r=F.ti_marginal(v,d,L,l)
z=solve(h,r,l,d,eps=1e-7,return_matrix=True)
b,q=z['input_rank'],z['output_rank'];u=z['U'];vv=z['V'];j=z['J']
t=np.kron(np.eye(d),u);x=(t.conj().T@r@t).reshape(d,b,d,b)
ext=np.zeros((d*d**l,d*d**l),complex)
for a in range(d):
    for ap in range(d):
        block=sum(x[a,i,ap,k]*j[i*q:(i+1)*q,k*q:(k+1)*q] for i in range(b) for k in range(b))
        ext[a*d**l:(a+1)*d**l,ap*d**l:(ap+1)*d**l]=vv@block@vv.conj().T
val=np.trace(np.kron(h,np.eye(d**(l-1)))@ext).real
assert abs(val-z['energy'])<1e-9
pet=F.petz_defect(h,v,L,l)
assert abs(pet-z['petz_energy'])<1e-9
print('PASS complex Choi objective and Petz compression',flush=True)
