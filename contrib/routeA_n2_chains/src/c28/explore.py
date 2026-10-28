"""C28 algebra, based on C22's bus, with unequal charge conversions."""
import json
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from bus_base import one_bus, data_ops, h2_g2, h2, quadratic_fit

PHASES=(.83,.31,-.07)
RATIO=.45

def coefficients(n, phases=PHASES, ratio=RATIO, order=4):
    A,(L,R)=one_bus(n,g=0,phases=phases,minus_ratio=ratio)
    C=one_bus(n,g=1,phases=phases,minus_ratio=ratio)[0]-A
    e,U=eigh(A); C=U.conj().T@C@U
    L=[U.conj().T@s for s in L];R=[U.conj().T@s for s in R]
    W=np.kron(R[0],L[0])+1j*np.kron(R[1],L[2])-1j*np.kron(R[2],L[1])
    den=(e[:,None]+e[None,:]).reshape(-1,1)
    X=W/den;out=[]
    for r in range(order+1):
        out.append(-(-1)**r*W.conj().T@X)
        Z=X.reshape(len(e),len(e),-1)
        X=(np.einsum('ij,jka->ika',C,Z)+np.einsum('ij,kja->kia',C,Z)).reshape(W.shape)/den
    return out

def normal_order(K, modes):
    cs,_,_,_=data_ops(modes);D=len(K)
    vacuum=np.eye(D)[:,0]
    one=np.column_stack([c.conj().T@vacuum for c in cs])
    const=K[0,0].real
    h=one.conj().T@K@one-const*np.eye(modes)
    Hq=const*np.eye(D,dtype=complex)
    for i in range(modes):
        for j in range(modes):Hq+=h[i,j]*cs[i].conj().T@cs[j]
    pairs=[(i,j) for i in range(modes) for j in range(i+1,modes)]
    two=np.column_stack([cs[i].conj().T@cs[j].conj().T@vacuum for i,j in pairs])
    vertex=two.conj().T@(K-Hq)@two
    H4=np.zeros_like(K)
    for a,(i,j) in enumerate(pairs):
        for b,(k,l) in enumerate(pairs):
            H4+=vertex[a,b]*cs[i].conj().T@cs[j].conj().T@cs[l]@cs[k]
    return const,h,vertex,float(np.linalg.norm(K-Hq-H4,2))

def symbol(h, fields,k):
    n=len(h)//2;T=h[:n,n:];M=h[:n,:n]+h[n:,n:]+np.diag(fields)
    return M+np.exp(1j*k)*T+np.exp(-1j*k)*T.conj().T

def roots(h,fields,grid=2001):
    n=len(fields);ks=np.linspace(-np.pi,np.pi,grid)
    es=np.array([np.linalg.eigvalsh(symbol(h,fields,k)) for k in ks])
    out=[]
    for b in range(n):
        for j in np.flatnonzero(es[:-1,b]*es[1:,b]<0):
            k=brentq(lambda k:np.linalg.eigvalsh(symbol(h,fields,k))[b],ks[j],ks[j+1],xtol=1e-14)
            e,U=eigh(symbol(h,fields,k));u=U[:,b]
            T=h[:n,n:]
            v=np.vdot(u,(1j*np.exp(1j*k)*T-1j*np.exp(-1j*k)*T.conj().T)@u).real
            out.append((k,v,u))
    return sorted(out,key=lambda x:x[0])

def wedge(x,y):return np.array([x[i]*y[j]-x[j]*y[i] for i in range(len(x)) for j in range(i+1,len(x))])

def scattering(U, out1,out2,in1,in2):
    def wave(r):
        k,v,u=r
        return np.concatenate([u,np.exp(1j*k)*u])
    return np.vdot(wedge(wave(out1),wave(out2)),U@wedge(wave(in1),wave(in2)))

def main():
    results=[]
    for n in (1,2):
        K=coefficients(n);c,h,V,err=normal_order(K[2],2*n)
        c4,h4,V4,err4=normal_order(K[4],2*n)
        print('n',n,'T',h[:n,n:],'M',h[:n,:n]+h[n:,n:],'quartic norm',np.linalg.norm(V4),flush=True)
        if n==2:
            # Cancel diagonal self energies, then scan two endpoint fields.
            M=h[:n,:n]+h[n:,n:];scale=np.linalg.norm(h[:n,n:],2)
            rng=np.random.default_rng(28);best=None
            for it in range(1000):
                offset=rng.uniform(-1.8,1.8,2)*scale
                fields=-M.diagonal().real+offset
                rr=roots(h,fields,101)
                if len(rr)==4:
                    score=min(abs(r[1]) for r in rr)
                    if best is None or score>best[0]:best=(score,fields,rr)
            if best:
                score,fields,rr=best
                print('BEST',score,fields,[(r[0],r[1]) for r in rr],flush=True)
            else:
                fields=-M.diagonal().real;rr=roots(h,fields)
                print('NO FOUR ROOTS',flush=True)
        else:
            fields=-(h[:1,:1]+h[1:,1:]).diagonal().real+np.array([abs(h[0,1])*.37]);rr=roots(h,fields)
        row=dict(n=n,phases=PHASES,ratio=RATIO,fields=fields.tolist(),h_real=h.real.tolist(),h_imag=h.imag.tolist(),
                 quartic_g2=float(np.linalg.norm(V)),degree_above_four_g4=err4,
                 roots=[dict(k=r[0],v=r[1]) for r in rr],
                 quartic_g4_real=V4.real.tolist(),quartic_g4_imag=V4.imag.tolist())
        if n==2 and len(rr)==4:
            right=[r for r in rr if r[1]>0];left=[r for r in rr if r[1]<0]
            row['pair_transfer_abs']=float(abs(scattering(V4,right[1],left[1],right[0],left[0])))
            row['pair_mismatch']=float(np.angle(np.exp(1j*(right[1][0]+left[1][0]-right[0][0]-left[0][0]))))
            row['forward_RL']=[[float(scattering(V4,a,b,a,b).real) for b in left] for a in right]
        results.append(row)
    with open('algebra_results.json','w') as f:json.dump(results,f,indent=2)

if __name__=='__main__':main()
