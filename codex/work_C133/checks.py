"""Deterministic controls and the C106 long-crossover Casimir test."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import eigh

checks={}
charges=np.repeat([-1.,0.,1.,2.],[8,19,15,4])
qsite=np.r_[np.tile(charges,96),np.zeros(48)]
assert len(qsite)==4464
charge_index=float(np.sum((qsite-qsite.mean())**2))
assert abs(charge_index-101664/31)<1e-9
checks['route_A_charge_index']=charge_index
vis=[]
for n in [1,2,4,8,16,32]:
    for q in [.2,.8,1.,1.25,5.]:
        x=2*np.arange(n+2)*np.log(q); shift=max(x)
        logS=shift+np.log(np.exp(x-shift).sum());V=np.exp((2*n+2)*np.log(q)-2*logS)
        assert V <= 1/(n+2)**2+1e-14
        vis.append([n,q,V*(n+2)**2])
checks['visibility_max_ratio']=max(x[2] for x in vis)
gram=[]
for n in [1,2,4,8,12]:
    a=4.;w=.25;delta=a-2*w
    A=a*np.eye(n)-w*(np.eye(n,k=1)+np.eye(n,k=-1));ev,U=eigh(A)
    B=U[0,:][None,:]*U
    def functions(t):
        diff=ev[:,None]-ev[None,:]
        F=np.empty_like(diff);same=np.abs(diff)<1e-12
        ep=np.broadcast_to(ev[:,None],diff.shape);eq=np.broadcast_to(ev[None,:],diff.shape)
        F[same]=t*np.exp(-ep[same]*t)
        F[~same]=(np.exp(-eq[~same]*t)-np.exp(-ep[~same]*t))/diff[~same]
        return np.einsum('jp,pq,jq->j',B,F,B)
    G,err=quad_vec(lambda t:np.outer(functions(t),functions(t)),0,np.inf,epsabs=1e-13,epsrel=1e-11)
    for t in [.1,.5,1.]:
        f=functions(t); assert np.max(np.abs(f.sum()-t*np.sum(U[0,:]**2*np.exp(-t*ev))))<1e-12
    ratio=float(np.trace(G)*4*delta**3);assert ratio<=1+1e-10
    gram.append(dict(n=n,trace=float(np.trace(G)),bound=1/(4*delta**3),ratio=ratio,integration_error=float(err)))
checks['gram']=gram
rows=[json.loads(s) for s in Path('ed_results.jsonl').read_text().splitlines()]
xx=[];ising=[];at=[]
for r in rows:
    if r['family']=='xxz' and r['parameter']==0:xx.append(abs(r['E']+1/np.sin(np.pi/r['L'])))
    if r['family']=='fz_parafermion' and r['d']==2:ising.append(abs(r['E']+2/np.sin(np.pi/(2*r['L']))))
    if r['family']=='ashkin_teller' and r['parameter']==0:at.append(abs(r['E']+4/np.sin(np.pi/(2*r['L']))))
checks['xx_energy_error']=max(xx);checks['ising_energy_error']=max(ising);checks['ashkin_teller_decoupled_error']=max(at)
assert max(max(xx),max(ising),max(at))<1e-8

def fermion_ring(L,eps,a=2.):
    values=[]
    for flux,parity in [(np.pi,0),(0.,1)]:
        k=(2*np.pi*np.arange(L)+flux)/L;root=np.sqrt(4*np.sin(k)**2+eps**2)
        energies=np.concatenate([-a*eps-root,-a*eps+root]);negative=energies<0
        E=float(energies[negative].sum());N=int(negative.sum())
        if N%2!=parity:E+=float(np.min(np.abs(energies)))
        values.append(E)
    return min(values), ('AP' if values[0]<=values[1]+1e-11 else 'P')

cross=[]
for ell in [2,4,8,16,32,64]:
    eps=2*np.sin(np.pi/(2*ell))/np.sqrt(3);M=np.sqrt(4+eps*eps);v=np.sqrt(3)*np.cos(np.pi/(2*ell))
    for factor in [4,8,16,32,64]:
        L=ell*factor;e,s=fermion_ring(L,eps);e2,s2=fermion_ring(2*L,eps)
        Q=L*(e2-2*e)/M
        cross.append(dict(ell=ell,epsilon=float(eps),L=L,M=float(M),v_over_M=float(v/M),A_exact_over_M=float(2*v/M),Q=float(Q),A_eff_over_M=float(4*Q/np.pi),sector=s,sector2=s2))
assert all(r['sector']=='AP' and r['sector2']=='AP' for r in cross)
checks['c106_max_relative_error_last']=max(abs(r['A_eff_over_M']/r['A_exact_over_M']-1) for r in cross if r['L']==64*r['ell'])
from scan import bonds,ground
mf=[]
for family,par,d,h,ls in bonds():
    if (family,par) in [('xxz',.5),('ashkin_teller',.5),('spin1_blbq',-np.pi/4)]:
        L=6 if d==2 else 4
        e,r,_=ground(h,d,L);em,rm,_=ground(h,d,L,matrixfree=True)
        assert abs(e-em)<1e-8
        mf.append(dict(family=family,L=L,energy_difference=abs(e-em),matrixfree_residual=rm))
checks['matrixfree_controls']=mf
F=np.diag(np.linspace(-.7,.9,7))+.12*(np.eye(7,k=1)+np.eye(7,k=-1))
a=np.arange(1.,8.);b=np.cos(np.arange(7.));T=np.outer(a,b)/20
ks=2*np.pi*np.arange(64)/64
values=np.array([np.linalg.det(F+np.exp(1j*k)*T+np.exp(-1j*k)*T.T) for k in ks])
coeff=np.fft.fft(values)/64
checks['rank_one_determinant_forbidden_coefficients']=float(np.max(np.abs(coeff[2:-1])))
assert checks['rank_one_determinant_forbidden_coefficients']<1e-12
Path('checks.json').write_text(json.dumps(checks,indent=2)+'\n');Path('c106_results.json').write_text(json.dumps(cross,indent=2)+'\n')
text=['# C133 checks','', '[numerical] Independent exact-energy and Gram integral controls passed.','',json.dumps(checks,indent=2),'','[numerical] C106 two-band free-chain test, t=1, a=2, k_F=pi/(2 ell), epsilon=2 sin(k_F)/sqrt(3). Both fermionic boundary sectors and their physical parity constraints were compared. At L a multiple of 4 ell the antiperiodic even-parity sector wins. This fixes shell contamination without tuning the Hamiltonian between L and 2L. The table uses L=64 ell.','', '| ell | epsilon | v/M | Q_L | A_eff/M | A_exact/M |','| --: | --: | --: | --: | --: | --: |']
for r in cross:
    if r['L']==64*r['ell']:text.append('| '+ ' | '.join(f'{r[k]:.9g}' for k in ['ell','epsilon','v_over_M','Q','A_eff_over_M','A_exact_over_M'])+' |')
Path('CHECKS.md').write_text('\n'.join(text)+'\n');print('\n'.join(text))

fixed=[]
r=.5
for ell in [2,4,8,16,32,64,128]:
    k=np.pi/(2*ell);ss=np.sin(k)**2;cc=1-ss;rr=r*r
    x=(rr*(1+ss)+np.sqrt(rr*rr*(1+ss)**2+4*(cc-rr)*rr*ss))/(2*(cc-rr))
    a=np.sqrt(1+x);eps=2*np.sin(k)/np.sqrt(x);M=np.sqrt(4+eps**2);v=2*np.sqrt(x/(1+x))*np.cos(k)
    assert abs(v/M-r)<1e-12
    L=64*ell;e,sector=fermion_ring(L,eps,a);e2,sector2=fermion_ring(2*L,eps,a)
    assert sector==sector2=='AP'
    Q=L*(e2-2*e)/M
    fixed.append(dict(ell=ell,a=float(a),epsilon=float(eps),L=L,v_over_M=float(v/M),A_exact_over_M=1.,Q=float(Q),A_eff_over_M=float(4*Q/np.pi)))
Path('c106_fixed_ratio.json').write_text(json.dumps(fixed,indent=2)+'\n')
with Path('CHECKS.md').open('a') as f:
    f.write('\n[numerical] A second C106 test adjusts a with epsilon to keep v/M=1/2 at every parameter, while keeping the same commensurate Fermi points. Hence A/M=1 throughout. At L=64 ell the computed values are\n\n| ell | epsilon | Q_L | A_eff/M |\n| --: | --: | --: | --: |\n')
    for z in fixed:f.write(f"| {z['ell']} | {z['epsilon']:.9g} | {z['Q']:.9g} | {z['A_eff_over_M']:.9g} |\n")
print('fixed ratio Casimir estimates',[z['A_eff_over_M'] for z in fixed])
