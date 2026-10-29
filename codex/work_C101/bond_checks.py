"""Small exact-spectrum consistency checks, no large ring diagonalizations."""
import json
import numpy as np
from itertools import product

def width(h):
    e=np.linalg.eigvalsh(h); return (e[-1]-e[0])/2

def hubbard(t,U,mu):
    cs=[]
    for a in range(4):
        c=np.zeros((16,16))
        for s in range(16):
            if (s>>a)&1:c[s^(1<<a),s]=(-1)**((s&((1<<a)-1)).bit_count())
        cs.append(c)
    ns=[c.T@c for c in cs]
    h=-t*sum(cs[a].T@cs[a+2]+cs[a+2].T@cs[a] for a in range(2))
    h+=U/2*(ns[0]@ns[1]+ns[2]@ns[3])-mu/2*sum(ns)
    return h

def hubbard_spectrum(t,U,mu):
    R=np.sqrt(4*t*t+U*U/16)
    return sorted([0,U-2*mu]+[-mu/2-t,-mu/2+t]*2+[-mu]*3+[U/2-mu,U/4-mu-R,U/4-mu+R]+[U/2-3*mu/2-t,U/2-3*mu/2+t]*2)

def tj(q,t,mu):
    d=q+1;h=np.zeros((d*d,d*d))
    for a,b in product(range(d),repeat=2):
        h[a*d+b,a*d+b]-=mu/2*((a>0)+(b>0))
        if a and b:
            h[b*d+a,a*d+b]+=t;h[a*d+b,a*d+b]-=t
        elif a!=b:h[b*d+a,a*d+b]-=t
    return h

def tj_M(t,mu):
    e=[0,-mu/2-t,-mu/2+t,-mu,-mu-2*t]
    return (max(e)-min(e))/2

def graded(par,fields,J):
    d=len(par);h=np.zeros((d*d,d*d));es=[]
    for a,b in product(range(d),repeat=2):
        h[b*d+a,a*d+b]+=J*(-1)**(par[a]*par[b])
        h[a*d+b,a*d+b]+=(fields[a]+fields[b])/2
    for a in range(d):
        es.append(J*(-1)**par[a]+fields[a])
        for b in range(a):es.extend([(fields[a]+fields[b])/2-J,(fields[a]+fields[b])/2+J])
    return h,sorted(es)

if __name__=='__main__':
    rng=np.random.default_rng(101);out={};err=0
    for _ in range(30):
        t,U,mu=rng.uniform(.1,2),rng.uniform(-10,10),rng.uniform(-6,6)
        h=hubbard(t,U,mu); es=hubbard_spectrum(t,U,mu)
        err=max(err,float(np.max(np.abs(np.linalg.eigvalsh(h)-es))))
        assert width(h)>=2*t-1e-10
    out['hubbard_spectrum_max_error']=err
    F=np.zeros((16,16));P=np.diag([1,-1,-1,1]);I=np.eye(4)
    for a,b in product(range(4),repeat=2):
        F[b+4*a,a+4*b]=(-1)**(a.bit_count()*b.bit_count())
    gauge_error=0
    for _ in range(10):
        h=hubbard(1,float(rng.normal()),float(rng.normal()))
        z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));o=z+z.conj().T
        oe=(o+P@o@P)/2
        D=np.kron(oe,I)-np.kron(I,oe)
        gauge_error=max(gauge_error,float(np.max(np.abs(F@h@F-h))),float(np.max(np.abs(F@D@F+D))))
        assert width(h+np.kron(o,I)-np.kron(I,o))>=width(h)-1e-10
    out['fermionic_reflection_error']=gauge_error
    err=0
    for q,mu in product(range(2,7),[-5,-2,-1,0,1,3]):
        err=max(err,abs(width(tj(q,1,mu))-tj_M(1,mu)))
    out['tj_width_max_error']=err
    err=0
    for d in range(2,8):
        par=rng.integers(0,2,d);fields=rng.normal(size=d);J=.7
        h,es=graded(par,fields,J)
        err=max(err,float(np.max(np.abs(np.linalg.eigvalsh(h)-es))))
    out['graded_spectrum_max_error']=err
    X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1,-1]);err=0
    for _ in range(30):
        jj=rng.normal(size=3);h=sum(j*np.kron(s,s)/4 for j,s in zip(jj,[X,Y,Z]))
        err=max(err,abs(width(h)-sum(sorted(abs(jj))[-2:])/4))
    out['xyz_width_max_error']=err
    assert max(out.values())<1e-10
    with open('bond_results.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
