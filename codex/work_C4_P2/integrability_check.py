import numpy as np

def comm(a,b):return a@b-b@a

def data(d,epsilon):
    p=np.zeros((d*d,d*d),dtype=np.int64)
    for a in range(d):
        for b in range(d):p[b*d+a,a*d+b]=1
    if epsilon==1:
        omega=np.eye(d,dtype=np.int64).reshape(-1)
    else:
        assert d%2==0
        j=np.zeros((d,d),dtype=np.int64)
        for a in range(d//2):j[a,a+d//2]=1;j[a+d//2,a]=-1
        omega=j.reshape(-1)
    k=np.outer(omega,omega)
    unit=np.eye(d,dtype=np.int64)
    s,t=np.kron(p,unit),np.kron(unit,p)
    e,f=np.kron(k,unit),np.kron(unit,k)
    a,b=t-s,f-e
    c0=comm(s+t,comm(s,t))
    c1=comm(e+f,comm(s,t))+comm(s+t,comm(e,t)+comm(s,f))
    c2=comm(s+t,comm(e,f))+comm(e+f,comm(e,t)+comm(s,f))
    c3=comm(e+f,comm(e,f))
    assert np.array_equal(c0,2*a)
    assert np.array_equal(c3,2*b)
    # Matrix element used in the all-dimension proof in P2_REPORT.md.
    if epsilon==1:
        bra,ket=(0,0,1),(2,1,2)
    else:
        bra,ket=(0,0,d//2),(1,0,d//2+1)
    index=lambda z:z[0]*d*d+z[1]*d+z[2]
    ij=(index(bra),index(ket))
    assert a[ij]==b[ij]==c0[ij]==c3[ij]==0
    assert c1[ij]==2 and c2[ij]==epsilon*d-2
    # Projection onto a,b using entries where only one is nonzero.
    ia=tuple(np.argwhere((a!=0)&(b==0))[0])
    ib=tuple(np.argwhere((a==0)&(b!=0))[0])
    def project(c):
        ca=int(c[ia])//int(a[ia]);cb=int(c[ib])//int(b[ib])
        assert ca*a[ia]==c[ia] and cb*b[ib]==c[ib]
        return ca,cb,c-ca*a-cb*b
    z1,w1,r1=project(c1);z2,w2,r2=project(c2)
    factor=d-2*epsilon
    assert np.array_equal(2*epsilon*r2,factor*r1)
    assert np.any(r1)
    print(f'd={d}, epsilon={epsilon}: C1 divergence=({z1},{w1}), C2 divergence=({z2},{w2}); 2 epsilon R2=(d-2 epsilon)R1; residual nonzero')
    # Stronger identity, retaining a canonical residual C1.
    # C2 - [(d-2epsilon)/(2epsilon)]C1 is a pure divergence.
    zz=2*epsilon*c2-factor*c1
    z,w,r=project(zz)
    assert not np.any(r)
    print('  2 epsilon C2-(d-2 epsilon)C1 =',z,'(P23-P12) +',w,'(K23-K12)')

if __name__=='__main__':
    for d in [3,4,5,6]:data(d,1)
    for d in [4,6]:data(d,-1)
