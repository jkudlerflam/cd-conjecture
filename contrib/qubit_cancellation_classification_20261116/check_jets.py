"""Check that unused microscopic Taylor coefficients do not remove witnesses."""
import json
import time
from pathlib import Path
import numpy as np
from scipy import sparse
from algebra import bond,chain,projectors,schur_series,maxabs


def modsq(z):
    return np.array([sum(z[k]*np.conj(z[n-k]) for k in range(n+1)).real for n in range(3)])


def main():
    start=time.perf_counter();rng=np.random.default_rng(20261117);rows=[];L=9;dim=1<<L
    def random():return rng.normal()+1j*rng.normal()
    for dense in (True,False):
        for q0,f0 in ([(0,np.sqrt(3)),(1,0),(3,2j*np.sqrt(3))] if dense else [(1,0),(3,0)]):
            for _ in range(4):
                qs=np.array([q0,random(),random()]);fs=np.array([f0,random(),random()])
                if dense and f0!=0:
                    # Tangency to the fourth-scalar cone, so order five is scalar.
                    fs[1]=(3*np.real(np.conj(q0-1)*qs[1])/abs(f0)**2+1j*rng.normal())*f0
                Q=modsq(qs);F=modsq(fs);shift=qs.copy();shift[0]-=1;R=modsq(shift)
                unit=np.array([1.,0.,0.])
                if dense:
                    C=.5*unit-qs.real+(Q+F)/6
                    B=-3.5*unit+3*qs.real-5*Q/6-2*F/3
                    tau=(R-F)/4+.5j*np.array([np.imag(3*fs[n]+sum(np.conj(qs[k])*fs[n-k]
                                                                             for k in range(n+1))) for n in range(3)])
                else:
                    C=-Q;B=4*unit+4*qs.real+2*Q;tau=-2*unit-2*qs.real
                eye=sparse.eye(dim,format='csr',dtype=complex)
                H=[chain(L,bond(v=1)),chain(L,bond(a=1,b=-1))]
                for order in range(2,7):
                    params={}
                    if order==2:params.update(u=-4 if dense else 2,t=-1)
                    if 2<=order<=4:params['p']=qs[order-2]
                    if dense and 3<=order<=5:params['f']=fs[order-3]
                    if order==6:params.update(f=random())
                    if order>=4:params.update(u=-B[order-4],t=-tau[order-4])
                    if order==6:
                        params['u']+=rng.normal();params['t']+=random()
                    h=chain(L,bond(**params))
                    if order>=4:h-=L*C[order-4]*eye
                    if dense and order==2:h+=2*L*eye
                    H.append(h)
                P=projectors(L);S=schur_series(H,P)
                keep=np.flatnonzero(abs(S[2].diagonal())<1e-9);K=schur_series(S[2:],keep)
                words=P[keep]
                if not dense:
                    keep=np.flatnonzero(abs(K[2].diagonal())<1e-9);K=schur_series(K[2:],keep);words=words[keep]
                    fifth=K[1];sixth=K[2];source=1;dest=1<<(4 if q0==1 else 3)
                    expected=-1 if q0==1 else -128/3
                else:
                    fifth=K[3];sixth=K[4];source=sum(1<<j for j in (0,2,4,7));dest=sum(1<<j for j in (0,3,5,7))
                    expected=-.5 if f0==0 else (15*q0*q0+4*q0-3)/16+.25j*(5*q0-1)*complex(f0).imag
                pos={int(s):i for i,s in enumerate(words)};actual=complex(sixth[pos[dest],pos[source]])
                assert maxabs(fifth)<1e-8 and abs(actual-expected)<1e-8
                rows.append({'dense':dense,'q':q0,'fifth_norm':maxabs(fifth),
                             'witness':[actual.real,actual.imag],'error':abs(actual-expected)})
    out={'seed':20261117,'L':L,'cases':rows,'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('jet_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Passed',len(rows),'higher-jet checks; maximum error',max(r['error'] for r in rows))


if __name__=='__main__':main()
