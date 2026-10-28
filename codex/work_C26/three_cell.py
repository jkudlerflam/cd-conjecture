"""Connected fourth-order adjacent-bond coefficient. Streaming middle-energy rows."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time
from pathlib import Path
import numpy as np
import scipy.linalg as la
import sidon_cell as sc
from checks_sidon import spectral,norm,tl,residual,pair_engine

def pair_data(e1,e2,ps1,ps2,N):
    d1,d2=len(e1),len(e2)
    b=sum(np.einsum('ia,jb->ijab',r[:,:N],l[:,:N]) for r,l in zip(ps1,ps2))
    energy=e1[:,None]+e2[None,:];rr=np.zeros_like(energy);np.divide(1,energy,out=rr,where=energy>0)
    w=rr[:,:,None,None]*b
    y=np.zeros_like(w)
    for r,l in zip(ps1,ps2):
        first=(r@w.reshape(d1,-1)).reshape(d1,d2,N,N)
        y+=(l@first.transpose(1,0,2,3).reshape(d2,-1)).reshape(d2,d1,N,N).transpose(1,0,2,3)
    bb=b.reshape(-1,N*N);ww=w.reshape(-1,N*N);yy=y.reshape(-1,N*N)
    A=bb.conj().T@ww;T=ww.conj().T@ww
    k4=-yy.conj().T@(rr.ravel()[:,None]*yy)+(A@T+T@A)/2
    return dict(w=w,y=y,A=A,T=T,K4=k4)

def run(n,scale=1.):
    start=time.time();sc.RHO_SCALE=scale;c=sc.reflection_cell(n,True);N=c['N'];rho2=c['rho']**2
    e1,f1,_=spectral(c,'R');e2,f2,_=spectral(c);e3,f3,_=spectral(c,'L')
    A=pair_data(e1,e2,f1['R'],f2['L'],N);B=pair_data(e2,e3,f2['R'],f3['L'],N)
    independent=pair_engine(c)
    pair_error=max(norm(A['K4']-independent['K4']),norm(B['K4']-independent['K4']))
    d1,d2,d3=len(e1),len(e2),len(e3)
    g1=np.eye(d1,N);g3=np.eye(d3,N)
    # V_B R V_A and V_A R V_B, respectively, before the remaining outer port action.
    wa=[np.einsum('jh,ihab->ijab',r,A['w'],optimize=True) for r in f2['R']]
    wb=[np.einsum('jh,hkbc->jkbc',l,B['w'],optimize=True) for l in f2['L']]
    ir=np.zeros((N**3,N**3),complex)
    for first in range(0,d2,8):
        ss=slice(first,min(first+8,d2))
        y=np.einsum('ijab,kc->ijkabc',A['y'][:,ss],g3)
        y+=np.einsum('ia,jkbc->ijkabc',g1,B['y'][ss])
        for s in range(3):
            y+=np.einsum('ijab,kc->ijkabc',wa[s][:,ss],f3['L'][s][:,:N])
            y+=np.einsum('ia,jkbc->ijkabc',f1['R'][s][:,:N],wb[s][ss])
        energy=e1[:,None,None]+e2[ss][None,:,None]+e3[None,None,:]
        rr=np.zeros_like(energy);np.divide(1,energy,out=rr,where=energy>0)
        yy=y.reshape(-1,N**3)
        ir-=yy.conj().T@(rr.ravel()[:,None]*yy)
    AA=np.kron(A['A'],np.eye(N))+np.kron(np.eye(N),B['A'])
    TT=np.kron(A['T'],np.eye(N))+np.kron(np.eye(N),B['T'])
    folded=(AA@TT+TT@AA)/2
    k4=ir+folded
    connected=k4-np.kron(A['K4'],np.eye(N))-np.kron(np.eye(N),B['K4'])
    par=sc.kron_all([c['logical_parity']]*3)
    rec=dict(n=n,rho_multiplier=scale,rho2=rho2,dimensions=[d1,d2,d3],physical_dimension=d1*d2*d3,
        code_dimension=N**3,pair_comparison_error=pair_error,
        full_K4_traceless_over_rho2=norm(tl(k4))/rho2,
        connected_scalar=float(np.trace(connected).real/N**3),connected_traceless_over_rho2=norm(tl(connected))/rho2,
        connected_nonGaussian_over_rho2=residual(connected,n,3)/rho2,
        hermiticity=norm(connected-connected.conj().T),parity_error=norm(connected@par-par@connected),seconds=time.time()-start)
    assert pair_error<5e-12
    assert rec['hermiticity']<5e-12 and rec['parity_error']<5e-12
    if scale==1:np.savez('three_cell_coefficients_n'+str(n)+'.npz',K4=k4,connected=connected,rho2=rho2)
    print(json.dumps(rec),flush=True)
    return rec

if __name__=='__main__':
    out=[]
    for n in [1,2]:
        for scale in [1.,.5,.25]:
            out.append(run(n,scale))
            Path('three_cell_results.json').write_text(json.dumps(out,indent=2))
