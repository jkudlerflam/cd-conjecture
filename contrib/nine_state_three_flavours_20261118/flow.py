"""Complete generic-momentum three-channel spinless one-loop equations."""
import itertools
import numpy as np
from physical import onsite
from mechanism import scalar,fermi_modes,vertex,bloch

EDGES=((0,1),(0,2),(1,2))


def determinant_coefficients(T,M):
    nodes=np.array([-1.,-.3,.3,1.])
    values=[np.linalg.det(bloch(T,M,np.arccos(x))).real for x in nodes]
    return np.linalg.solve(np.polynomial.polynomial.polyvander(nodes,3),values)


def match(T,M,quartic,require_generic=True):
    roots=np.polynomial.polynomial.polyroots(determinant_coefficients(T,M))
    assert np.max(abs(roots.imag))<1e-8 and max(abs(roots.real))<1
    modes=fermi_modes(T,M,roots.real)
    channels=[(i,5-i) if modes[i]['v']>0 else (5-i,i) for i in range(3)]
    v=np.array([modes[R]['v'] for R,L in channels]);assert min(v)>0
    def vertex_at(a,b,c,d):return vertex(quartic['W'],quartic['pairs'],modes,a,b,c,d)
    F=np.array([[vertex_at(Ri,Lj,Ri,Lj).real for Rj,Lj in channels] for Ri,Li in channels])
    J=np.array([[vertex_at(Ri,Li,Rj,Lj) for Rj,Lj in channels] for Ri,Li in channels])
    assert np.max(abs(J.imag))<1e-10 and np.max(abs(F-F.T))<1e-10
    J=J.real
    b=np.diag(F)/(2*np.pi*v);f=F/(np.pi*(v[:,None]+v[None,:]));np.fill_diagonal(f,0)
    C=J/(2*np.pi*np.sqrt(v[:,None]*v[None,:]));np.fill_diagonal(C,0)
    D=b[:,None]+b[None,:]-2*f
    a=np.array([4*v[i]*v[j]/(v[i]+v[j])**2 for i,j in EDGES])
    d=np.array([D[i,j] for i,j in EDGES]);c=np.array([C[i,j] for i,j in EDGES])
    pairs=list(itertools.combinations_with_replacement(range(6),2));collisions=[];sep=[]
    for ab,cd in itertools.combinations(pairs,2):
        delta=abs(np.exp(1j*sum(modes[i]['k'] for i in ab))-np.exp(1j*sum(modes[i]['k'] for i in cd)))
        if delta<1e-9:collisions.append([ab,cd])
        else:sep.append(delta)
    generic=len(collisions)==3 and all(sum(ab)==sum(cd)==5 for ab,cd in collisions)
    if require_generic:assert generic
    allowed=[]
    for ab in itertools.combinations(range(6),2):
        for cd in itertools.combinations(range(6),2):
            if abs(np.exp(1j*sum(modes[i]['k'] for i in ab))-np.exp(1j*sum(modes[i]['k'] for i in cd)))<1e-9:
                allowed.append(dict(outgoing=ab,incoming=cd,value=vertex_at(*ab,*cd)))
    return dict(cos_roots=roots.real,modes=modes,channels=channels,velocities=v,forward_RL=F,pair_matrix=J,
                b=b,f=f,D=d,c=c,a=a,allowed_quartics=allowed,pair_sum_collisions=collisions,
                other_pair_sum_chord_margin=min(sep),generic_momentum_channels=generic)


def beta(t,y,a):
    """y=(D12,D13,D23,c12,c13,c23), time tau=g² ell."""
    d,c=y[:3],y[3:]
    dd=-(1+2*a)*c*c-np.sum(c*c)
    dc=-d*c-np.array([c[1]*c[2],c[0]*c[2],c[0]*c[1]])
    return np.r_[dd,dc]


def complete_beta(b,f,c,a):
    """Independent 3x3 matrix implementation, also retaining forward terms."""
    db=-np.sum(c*c,axis=1)
    df=a*c*c
    D=b[:,None]+b[None,:]-2*f
    dc=-D*c-c@c;np.fill_diagonal(dc,0)
    return db,df,dc
