"""Finite Gram matrices and n=1,2 bands of the positive mobile core."""
import itertools
import json
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import brentq
import mobile as m


def gram(n,mass=1.):
    A=np.diag(np.full(n,mass+2.))
    if n>1:
        A[0,0]+=1; A[-1,-1]-=1
        A+=np.diag(-np.ones(n-1),1)+np.diag(-np.ones(n-1),-1)
    E=np.zeros((n,n)); E[0,0]=1
    B=np.block([[np.kron(np.eye(n),A),-np.kron(np.eye(n),E)],
                [np.zeros((n*n,n*n)),np.kron(np.eye(n),A)]])
    v=np.r_[np.zeros(n*n),np.eye(n).ravel(order='F')]
    P=solve_continuous_lyapunov(B,np.outer(v,v))
    ids=[j*(n+1) for j in range(n)]
    G=P[np.ix_(ids,ids)]
    return G,float(np.linalg.norm(B@P+P@B.T-np.outer(v,v)))


def fit(H,n):
    mm=2*n; ops=[]
    for j in range(mm):
        O=np.ones((1,1))
        for a in range(mm):
            fac=np.diag([1.,-1.]) if a<j else np.array([[0,1],[0,0]]) if a==j else np.eye(2)
            O=np.kron(O,fac)
        ops.append(O)
    X=np.array([np.eye(2**mm).ravel()]+[(a.T@b).ravel() for a in ops for b in ops]).T
    coef=np.linalg.lstsq(X,H.ravel(),rcond=None)[0]
    h=coef[1:].reshape(mm,mm); h=(h+h.conj().T)/2
    return h[:n,n:],h[:n,:n]+h[n:,n:],float(np.linalg.norm(H-(X@coef).reshape(H.shape),2))


def bands(n):
    c=m.core(n,g=.05,mass=.2,express=1.,loop=.001)
    me=m.measure(c); H=m.bond(me,me); T,mu,ng=fit(H,n)
    phase=np.angle(T[0,0]); R=T*np.exp(-1j*phase)
    gamma=np.linalg.eigvalsh((R+R.conj().T)/2)
    field=gamma[0]/2; mu+=field*np.eye(n)
    scale=np.max(abs(T)); T/=scale; mu/=scale
    def vals(k):return np.linalg.eigvalsh(mu+T*np.exp(1j*k)+T.conj().T*np.exp(-1j*k))
    grid=np.linspace(-np.pi,np.pi,6001); energies=np.array([vals(k) for k in grid])
    roots=[]; centres=[]; velocities=[]
    for a in range(n):
        rr=[brentq(lambda k:vals(k)[a],grid[j],grid[j+1],xtol=1e-14)
            for j in range(len(grid)-1) if energies[j,a]*energies[j+1,a]<0]
        roots.extend(rr); centres.append(float(sum(rr)/2))
        for k in rr:
            e,U=np.linalg.eigh(mu+T*np.exp(1j*k)+T.conj().T*np.exp(-1j*k))
            v=U[:,a]; der=1j*(T*np.exp(1j*k)-T.conj().T*np.exp(-1j*k))
            velocities.append(float(abs(np.vdot(v,der@v))))
    dist=lambda x:abs((x+np.pi)%(2*np.pi)-np.pi)
    d2=min(dist(a-b) for a,b in itertools.combinations(roots,2))
    pairs=list(itertools.combinations_with_replacement(range(len(roots)),2))
    d4=min(dist(roots[a]+roots[b]-roots[c]-roots[d])
           for (a,b),(c,d) in itertools.combinations(pairs,2))
    return {'n':n,'core_dimension':c['d'],'T_scale':float(scale),
            'phase_rotated_eigenvalues_over_T':(gamma/scale).tolist(),
            'field':float(field),'roots':roots,'centres_mod_pi':centres,
            'delta2':float(d2),'delta4':float(d4),
            'minimum_velocity_over_T':min(velocities),'nonGaussian_over_T':ng/scale}


if __name__=='__main__':
    out={'status':'numerical','gram':[],'bands':[]}
    for n in range(1,6):
        G,res=gram(n); ev=np.linalg.eigvalsh(G)
        rho=float(np.linalg.det(G)/np.trace(G)**(n-1))
        out['gram'].append({'n':n,'lyapunov_dimension':2*n*n,'residual':res,
                            'eigenvalues':ev.tolist(),'determinant_lower_bound':rho,
                            'smallest_over_largest':float(ev[0]/ev[-1])})
    for n in (1,2):out['bands'].append(bands(n))
    print(json.dumps(out,indent=2))
    open('band_checks.json','w').write(json.dumps(out,indent=2)+'\n')
