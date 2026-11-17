#!/usr/bin/env python3
"""Independent literal Pauli matrices, gauge and ED/DMRG observable gates."""
import json,time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from common import HERE,Z,build,product_expectation


def literal(L,b,half_ends=False):
    X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);z=np.diag([1.,-1.]);I=np.eye(2)
    h=.5*(np.kron(X,X)-np.kron(Y,Y))+b['D']*(np.kron(X,Y)-np.kron(Y,X))+b['Jz']*np.kron(z,z)
    H=np.zeros((2**L,2**L),complex)
    for j in range(L-1):H+=np.kron(np.eye(2**j),np.kron(h,np.eye(2**(L-j-2))))
    for j in range(L):H+=(.5 if half_ends and j in (0,L-1) else 1)*b['hz']*np.kron(np.eye(2**j),np.kron(z,np.eye(2**(L-j-1))))
    return H


def dense_mpo(mpo):
    terms={0:np.ones((1,1))}
    for W in mpo.W:
        nxt={}
        for a,T in terms.items():
            for b in range(mpo.w):
                if np.any(W[a,b]):nxt[b]=nxt.get(b,0)+np.kron(T,W[a,b])
        terms=nxt
    return terms[2]+mpo.const*np.eye(2**mpo.L)


def main():
    start=time.time();rows=[];L=8
    with threadpool_limits(limits=1):
        for cfg in [dict(L=L,V=4,a=1),dict(L=L,V=8,a=1),dict(L=L,V=12,a=.96),dict(L=L,V=40,a=1,hz=-19.976047250613814)]:
            mpo,gph,shift,b=build(cfg);H=literal(L,b);ph=np.ones(1,complex)
            for g in gph:
                th=-np.angle(g);ph=np.kron(ph,[np.exp(-.5j*th),np.exp(.5j*th)])
            Hrot=ph[:,None]*H*ph.conj()[None,:]
            me=float(np.max(abs(dense_mpo(mpo)-L*shift*np.eye(2**L)-Hrot)))
            enderr=float(np.linalg.norm(H-literal(L,b,True),2))
            assert me<1e-11 and abs(enderr-abs(b['hz']))<1e-11
            for P in (0,1):
                ix=[i for i in range(2**L) if i.bit_count()%2==P]
                e,v=eigh(H[np.ix_(ix,ix)],subset_by_index=[0,0]);ve=np.zeros(2**L,complex);ve[ix]=v[:,0]
                Se=[Z.entropy(np.linalg.svd(ve.reshape(2**l,-1),compute_uv=False)**2) for l in range(1,L)]
                Ze=[sum(abs(ve[s])**2*(1-2*((s>>(L-j-1))&1)) for s in range(2**L)) for j in range(L)]
                psi=Z.random_mps(L,P,D0=8,seed=3,dtype=float)
                psi,hist=Z.dmrg(psi,mpo,[dict(chi=32,alpha=a,res_tol=1e-11,maxit=300) for a in (1e-3,0,0,0)])
                S=Z.entropies(psi);m=Z.measure(psi,mpo,gph=gph,gwin=4,cp_rmax=2)
                errE=abs(hist[-1]['E']+mpo.const-L*shift-e[0]);errS=float(np.max(abs(S-Se)));errZ=float(np.max(abs(np.array(m['Z'])-Ze)))
                # Original-frame pair expectation at separation two, independent ED product.
                ops={2:np.conj(gph[2])*Z.SM,4:np.conj(gph[4])*Z.SM};pv=product_expectation(psi,ops)
                O=np.ones((1,1))
                for j in range(L):O=np.kron(O,Z.SM if j in (2,4) else Z.I2)
                pe=np.vdot(ve,O@ve);errp=abs(pv-pe)
                assert max(errE,errS,errZ,errp)<1e-8
                rows.append(dict(V=b['V'],a=b['a'],P=P,L=L,matrix_error=me,endpoint_half_error=enderr,
                                 energy_error=errE,entropy_error=errS,density_Z_error=errZ,pair_error=errp,residual=hist[-1]['res_max']))
                print(rows[-1],flush=True)
    (HERE/'validation.json').write_text(json.dumps(dict(rows=rows,runtime_seconds=time.time()-start),indent=2)+'\n')


if __name__=='__main__':main()
