#!/usr/bin/env python3
"""Independent ED checks for all new observables, including overlapping bonds."""
import json,time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from phase_common import HERE,Z,build,write_json
from observables import environments,operator_matrix,pair_correlations


def embed(O,i,width,L):return np.kron(np.eye(2**i),np.kron(O,np.eye(2**(L-i-width))))


def main():
    start=time.time();rows=[];L=8
    with threadpool_limits(limits=1):
        mpo,gph,shift,_=build(dict(L=L));n=np.diag([0.,1.]);zz=np.kron(Z.ZZ_,Z.ZZ_)
        H=sum(embed(b,j,2,L) for j,b in enumerate(mpo.h2))+sum(embed(h,j,1,L) for j,h in enumerate(mpo.h1_in))
        for P in (0,1):
            ix=[i for i in range(2**L) if i.bit_count()%2==P];e,v=eigh(H[np.ix_(ix,ix)],subset_by_index=[0,0])
            ve=np.zeros(2**L,dtype=complex);ve[ix]=v[:,0]
            psi=Z.random_mps(L,P,D0=8,seed=17,dtype=float)
            psi,hist=Z.dmrg(psi,mpo,[dict(chi=32,alpha=a,res_tol=1e-12,maxit=300) for a in (1e-3,0,0,0)])
            left=environments(psi);errors={}
            pol=np.kron(Z.PX,Z.PX)+np.kron(Z.PY,Z.PY)
            pops=[Z.rotate_bond(pol,gph[j],gph[j+1]) for j in range(L-1)]
            for name,ops,w in [('density',[n]*L,1),('bond',mpo.h2,2),('kinetic',[b-zz for b in mpo.h2],2),('valley_bond',pops,2)]:
                ids=list(range(L-w+1));got=operator_matrix(psi,left,ops,w,ids)
                vecs=np.array([embed(O,j,w,L)@ve for j,O in enumerate(ops)])
                means=np.real(vecs.conj()@ve);raw=np.real(vecs.conj()@vecs.T)
                conn=raw-np.outer(means,means)
                err=max(np.max(abs(np.array(got['profile'])-means)),np.max(abs(np.array(got['raw'])-raw)),
                        np.max(abs(np.array(got['connected'])-conn)))
                assert err<1e-8;errors[name]=float(err)
            pairs=pair_correlations(psi,left,gph,0,L)
            A0=np.kron(np.kron(Z.SP,Z.I2),Z.SP);errp=0.
            for row in pairs:
                i=row['i'];Ai=gph[i]*gph[i+2]*embed(A0,i,3,L);mi=np.vdot(ve,Ai@ve)
                for x in row['correlations']:
                    j=x['j'];Aj=gph[j]*gph[j+2]*embed(A0,j,3,L);mj=np.vdot(ve,Aj@ve)
                    raw=np.vdot(Ai@ve,Aj@ve);conn=raw-np.conj(mi)*mj
                    errp=max(errp,abs(complex(*x['raw'])-raw),abs(complex(*x['connected'])-conn))
            assert errp<1e-8;errors['pair']=float(errp)
            energy_error=abs(hist[-1]['E']+mpo.const-e[0]);assert energy_error<1e-8
            rows.append(dict(P=P,L=L,errors=errors,energy_error=energy_error))
    write_json(HERE/'validation.json',dict(rows=rows,runtime_seconds=time.time()-start))
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
