"""All-sector sixth-order checks after the dilute fourth constraint.

The complex-generic family is left as a Schur expression; four slices have
explicit local formulas here. Independent epsilon spectra check q=1 and 3.
"""
import json
import time
from pathlib import Path
import numpy as np

from algebra import dilute_higher,maxabs,occupations


def local_sixth(L,P,q):
    ns=occupations(L)[P];pos={int(s):i for i,s in enumerate(P)}
    N=np.sum(ns,axis=1)
    D={r:np.sum(ns*np.roll(ns,-r,axis=1),axis=1) for r in range(3,9)}
    if q==0:
        diag=-16*N+25/3*D[5]+25/12*D[6]+7/12*D[7]+D[8]/8
        hops={1:8.}
    elif q==-1:
        diag=2*L-4*N+256/3*D[5]+128/3*D[6]+256/11*D[7]+8*D[8]
        hops={}
    elif q==1:
        diag=2*L-50*N-2*D[4]+10/3*D[5]+14/3*D[6]
        hops={1:24.,4:-1.}
    elif q==3:
        gap8=np.sum(ns*np.roll(ns,-8,axis=1)*(1-np.roll(ns,-3,axis=1))
                    *(1-np.roll(ns,-5,axis=1)),axis=1)
        diag=18*L-332/3*N-256/3*D[3]+8*gap8
        hops={1:80.,3:-128/3}
    else:raise ValueError(q)
    M=np.diag(diag).astype(complex)
    for col,state in enumerate(P):
        for r,t in hops.items():
            for j in range(L):
                k=(j+r)%L
                if ns[col,j] and not ns[col,k]:
                    target=state^(1<<j)^(1<<k)
                    if target in pos:
                        M[pos[target],col]+=t;M[col,pos[target]]+=np.conj(t)
    return M


def main():
    start=time.perf_counter();rows=[];spectra=[]
    for L in (10,11,12,13,15,16):
        for q in (0,-1,1,3):
            P,K,H=dilute_higher(L,q)
            err=maxabs(K[2]-local_sixth(L,P,q))
            assert maxabs(K[1])<1e-9 and err<1e-8
            rows.append({'L':L,'q':q,'dimension':len(P),'fifth_order_norm':maxabs(K[1]),
                         'sixth_full_matrix_error':err})
            print(rows[-1],flush=True)
            if L==10 and q in (1,3):
                epsrows=[];lead=K[2].toarray()
                for eps in (.075,.05,.035,.025):
                    full=sum(eps**n*h for n,h in enumerate(H)).toarray()
                    ev,U=np.linalg.eigh(full);m=len(P)
                    a,sv,b=np.linalg.svd(U[P,:m]);polar=a@b
                    effective=(polar*ev[:m])@polar.conj().T
                    epsrows.append({'epsilon':eps,'error_over_epsilon6':maxabs(effective/eps**6-lead),
                                    'min_overlap':float(min(sv)),
                                    'band_gap_over_epsilon4':float((ev[m]-ev[m-1])/eps**4)})
                    # A separated finite-ring band is required. Its prefactor
                    # need not exceed one; q=3 at .075 has a small positive gap.
                    assert min(sv)>.7 and ev[m]-ev[m-1]>0
                assert epsrows[-1]['error_over_epsilon6']<.35*epsrows[0]['error_over_epsilon6']
                spectra.append({'L':L,'q':q,'rows':epsrows})
    generic=[]
    for q in (1j,.5+1j):
        for L in (12,13):
            P,K,_=dilute_higher(L,q);pos={int(s):i for i,s in enumerate(P)}
            numbers=np.array([int(s).bit_count() for s in P]);M=K[2].toarray()
            generic.append({'L':L,'q':[q.real,q.imag],'dimension':len(P),
                            'fifth_order_norm':maxabs(K[1]),
                            'number_changing_sixth_norm':maxabs(M*(numbers[:,None]!=numbers[None,:])),
                            'one_particle_hops':{str(r):[float(M[pos[1<<r],pos[1]].real),
                                                          float(M[pos[1<<r],pos[1]].imag)]
                                                 for r in range(1,L//2+1)}})
    out={'full_matrix_checks':rows,'direct_spectra':spectra,'generic_diagnostics':generic,
         'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('dilute_sixth_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
