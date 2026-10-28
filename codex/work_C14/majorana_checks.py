"""Finite-dimensional checks of the Majorana proof's norm and doubling factors."""
import json
import numpy as np

def gammas(n):
    I=np.eye(2); Z=np.diag([1.,-1.]); X=np.array([[0.,1.],[1.,0.]]); Y=np.array([[0.,-1j],[1j,0.]])
    out=[]
    for j in range(n):
      for W in [X,Y]:
        t=np.array([[1.]])
        for a in range(n): t=np.kron(t,Z if a<j else W if a==j else I)
        out.append(t)
    return out

def main():
    rng=np.random.default_rng(1401); out={'cross_norm':[]}
    for m in [1,2,3]:
      g=gammas(2*m); B=rng.normal(size=(2*m,2*m)); dim=2**(2*m)
      cross=sum(.5j*B[a,b]*(g[a]@g[2*m+b]) for a in range(2*m) for b in range(2*m))
      eig=np.linalg.eigvalsh(cross); pred=np.linalg.svd(B,compute_uv=False).sum()/2
      d=2**m; o=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));o=(o+o.conj().T)/2
      gauged=cross+np.kron(o,np.eye(d))-np.kron(np.eye(d),o)
      P=np.diag([(-1.)**int(j).bit_count() for j in range(d)])
      PT=np.kron(P,P);P1=np.kron(P,np.eye(d))
      even=(gauged+PT@gauged@PT)/2; extracted=(even-P1@even@P1)/2
      out['cross_norm'].append(dict(m=m,norm=np.ptp(eig)/2,nuclear_half=pred,
        gauge_spread=np.ptp(np.linalg.eigvalsh(gauged))/2,
        extraction_error=float(np.max(abs(extracted-cross)))))
    g=gammas(4); left=[g[0],g[1]];el=[g[2],g[3]];right=[g[4],g[5]];er=[g[6],g[7]]
    fl=[(a+1j*b)/2 for a,b in zip(left,el)];fr=[(a+1j*b)/2 for a,b in zip(right,er)]
    B=rng.normal(size=(2,2))
    hd=sum(.5j*B[a,b]*(left[a]@right[b]+el[a]@er[b]) for a in range(2) for b in range(2))
    hf=sum(1j*B[a,b]*(fl[a].conj().T@fr[b]-fr[b].conj().T@fl[a]) for a in range(2) for b in range(2))
    Q=sum(f.conj().T@f for f in fl+fr)
    out['doubling']=dict(matrix_error=float(np.max(abs(hd-hf))),
       number_commutator=float(np.max(abs(hd@Q-Q@hd))))
    assert all(abs(r['norm']-r['nuclear_half'])<1e-10 and r['extraction_error']<1e-10 for r in out['cross_norm'])
    assert out['doubling']['matrix_error']<1e-10 and out['doubling']['number_commutator']<1e-10
    open('majorana_checks.json','w').write(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__': main()
