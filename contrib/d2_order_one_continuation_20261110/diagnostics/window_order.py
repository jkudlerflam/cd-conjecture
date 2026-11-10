#!/usr/bin/env python3
"""Whole-window order measurements on frozen checkpoints, without optimization."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import json
import time
import numpy as np
import analyze_saved as A
Z=A.Z


def window_square(psi,ops,weights):
    # Transfer environments for I, Q and Q†Q. Pauli operators obey O_j²=I.
    e0=np.ones((1,1),dtype=complex)
    e1=np.zeros((1,1),dtype=complex)
    e2=np.zeros((1,1),dtype=complex)
    for m,op,f in zip(psi.M,ops,weights):
        n0=Z._apply_site(e0,m,None)
        n1=Z._apply_site(e1,m,None)
        n2=Z._apply_site(e2,m,None)
        if abs(f)>0:
            cross=Z._apply_site(e1,m,op)
            n2+=np.conj(f)*cross+f*cross.conj().T+abs(f)**2*n0
            n1+=f*Z._apply_site(e0,m,op)
        e0,e1,e2=n0,n1,n2
    norm=e0.item()
    assert abs(norm.imag)<1e-10 and abs(norm.real-1)<1e-9
    assert abs(e2.item().imag)<1e-8
    return float((e2/norm).item().real),complex((e1/norm).item())


def check_tiny():
    p=Z.random_mps(6,1,D0=4,seed=129,dtype=complex)
    v=p.M[0][0]
    for m in p.M[1:]:v=np.tensordot(v,m,axes=([-1],[0]))
    v=v.reshape(-1)
    ops=[np.cos(.31*j)*Z.PX+np.sin(.31*j)*Z.PY for j in range(6)]
    f=np.exp(.37j*np.arange(6));f[[0,5]]=0
    Q=np.zeros((64,64),dtype=complex)
    for j in range(6):
        term=np.array([[1.]])
        for k in range(6):term=np.kron(term,ops[k] if k==j else Z.I2)
        Q+=f[j]*term
    val,mean=window_square(p,ops,f)
    exact=np.vdot(Q@v,Q@v).real
    assert abs(val-exact)<1e-11 and abs(mean-np.vdot(v,Q@v))<1e-11
    return dict(transfer=val,exact_state_vector=float(exact),absolute_error=abs(val-exact))


def main():
    validation=check_tiny()
    rows=json.loads((A.HERE/'saved_diagnostics.json').read_text())
    selected=[]
    for key in sorted(set((r['name'],r['L']) for r in rows)):
        selected.append(min((r for r in rows if (r['name'],r['L'])==key),key=lambda r:r['E']))
    output=dict(tiny_validation=validation,rows=[])
    for r in selected:
        start=time.time();raw=json.loads((A.ROOT/r['source_record']).read_text());L=r['L'];W=L//2
        p=Z.MPS.load(r['checkpoint']);_,g,_=A.D.build_model(L,raw['couplings'],raw.get('distance',1))
        ops={'X':[gj*Z.SP+gj.conjugate()*Z.SM for gj in g],
             'Y':[-1j*gj*Z.SP+1j*gj.conjugate()*Z.SM for gj in g]}
        dist,y=np.array(raw['measurements']['long_connected_correlations']['Y']).T
        qs=np.linspace(0,np.pi,2001)
        q=float(qs[np.argmax(np.abs(np.exp(-1j*np.outer(qs,dist))@y))])
        result=dict(name=r['name'],L=L,parity=r['parity'],source_record=r['source_record'],
                    checkpoint_sha256=r['checkpoint_sha256'],window=[L//4,3*L//4],channels={})
        for label,axis,momentum in [('X_pi','X',np.pi),('Y_zero','Y',0.),('Y_adaptive','Y',q)]:
            weights=np.exp(1j*momentum*np.arange(L));weights[:L//4]=0;weights[3*L//4:]=0
            sq,mean=window_square(p,ops[axis],weights)
            result['channels'][label]=dict(axis=axis,q=momentum,square=sq,structure_factor=sq/W,
                order_square=sq/W**2,mean_real=mean.real,mean_imag=mean.imag)
        if abs(raw['couplings'].get('Jx',0)+raw['couplings'].get('Jy',0))<1e-14:
            a=result['channels']['X_pi']['square'];b=result['channels']['Y_zero']['square']
            assert abs(a-b)<1e-7
        result['seconds']=time.time()-start;output['rows'].append(result)
        (A.HERE/'window_order.json').write_text(json.dumps(output,indent=2)+'\n')
        print(r['name'],L,r['parity'],result['channels'],flush=True)


if __name__=='__main__':main()
