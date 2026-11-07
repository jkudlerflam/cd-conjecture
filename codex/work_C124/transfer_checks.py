"""Independent full-matrix evaluation of the thermal-transfer contractions."""
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh
from checks import ROOT, embed_bond, histories_dense


def main():
    out=[]
    for model,d,delta,L,M in [('XX',2,0.,4,.5),('Heisenberg',2,1.,4,.5),
            ('XXZ_m09',2,-.9,4,.5),('XXZ_p05',2,.5,4,.5),('SU3',3,1.,3,1.)]:
        h=np.zeros((d*d,d*d))
        for a,b in product(range(d),repeat=2):
            if d==3: h[d*b+a,d*a+b]=1.
            else:
                h[d*a+b,d*a+b]=delta*(.5-a)*(.5-b)
                if a!=b: h[d*b+a,d*a+b]=.5
        H={n:sum((embed_bond(h,n,d,i,(i+1)%n) for i in range(n)),
                 np.zeros((d**n,d**n),complex)) for n in [L,2*L]}
        ev,v=eigh(H[L]); evd,vd=eigh(H[2*L])
        tau=L/M
        e=(v*np.exp(-tau*(ev-ev[0])/2))@v.conj().T
        bD=(vd*np.exp(-tau*(evd-evd[0])/2))@vd.conj().T
        bS=np.kron(e,e)
        ix=np.arange(d**(2*L))
        a=(ix//(d**L))%d; b=ix%d
        labD=b+d*a
        maskS=(a[:,None]==a[None,:])&(a[:,None]==b[:,None])&(a[None,:]==b[None,:])
        maskD=labD[:,None]==labD[None,:]
        us=bS*maskS; ud=bD*maskD
        xs=us.copy(); xd=ud.copy()
        previous=None
        # Independent history enumeration only on the small numerator ring.
        for N in range(1,5):
            sn=float(np.vdot(us,xs).real); dn=float(np.vdot(ud,xd).real)
            q=histories_dense(H[L]-ev[0]*np.eye(d**L),d,L,N*tau,N,[0])
            err=float(abs(sn-q@q)/sn)
            assert err<2e-10,(model,N,err)
            G=N*tau*(evd[0]-2*ev[0])+np.log(sn/dn)
            row=dict(model=model,L=L,N=N,G=float(G),numerator_error=err)
            if previous is not None: row['increment']=float(G-previous)
            out.append(row)
            previous=G
            xs=(bS@xs@bS)*maskS
            xd=(bD@xd@bD)*maskD
        print(model,flush=True)
    # Independent diagonal matrix-unit calculation of the counterexample.
    for L in [2,3,4]:
        eps=L**-2; tau=2*L
        levels=np.array([0.,eps,1.])
        n=2*L
        states=np.arange(3**n)
        digits=(states[:,None]//(3**np.arange(n))[None,:])%3
        energy=levels[digits].sum(axis=1)
        keep=digits[:,0]==digits[:,L]
        norm=float(np.exp(-tau*energy[keep]).sum())
        z=lambda t:float(np.exp(-t*levels).sum())
        expected=z(2*tau)*z(tau)**(2*L-2)
        assert abs(norm/expected-1)<1e-12
        # Only the all-zero matrix unit has eigenvalue one; its u entry is one.
        assert np.count_nonzero(energy==0)==1
        out.append(dict(model='onsite',L=L,boundary=float(np.log(norm)),error=abs(norm/expected-1)))
    (ROOT/'transfer_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('complete',len(out))


if __name__=='__main__': main()
