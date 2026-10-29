"""Zero-field repulsive Hubbard Lieb-Wu quadrature, t=1.
Finite numerical check, not an all-parameter proof.
"""
import json
import numpy as np
from scipy.special import digamma,iv
from scipy.linalg import solve
from numpy.polynomial.legendre import leggauss
from bond_checks import hubbard_spectrum

def kernel(x,u):return np.real(digamma(1+1j*x/(4*u))-digamma(.5+1j*x/(4*u)))/(4*np.pi*u)
def point(U,Q,N):
    u=U/4;x,w=leggauss(N);k=Q*x;w=Q*w;s=np.sin(k);c=np.cos(k)
    R=kernel(s[:,None]-s[None,:],u)
    rho=solve(np.eye(N)-c[:,None]*R*w[None,:],np.ones(N)/(2*np.pi))
    mat=np.eye(N)-R*(w*c)[None,:]
    f,g=solve(mat,np.column_stack([-2*c,np.ones(N)])).T
    edge=kernel(np.sin(Q)-s,u)*w*c
    mu=(-2*np.cos(Q)+edge@f)/(1+edge@g)
    eps=f-mu*g
    def ep(z):return -2*np.cos(z)-mu+(kernel(np.sin(z)-s,u)*w*c)@eps
    delta=1e-5
    derivative=(ep(Q-2*delta)-8*ep(Q-delta)+8*ep(Q+delta)-ep(Q+2*delta))/(12*delta)
    rhoQ=1/(2*np.pi)+np.cos(Q)*(kernel(np.sin(Q)-s,u)*w)@rho
    vc=derivative/(2*np.pi*rhoQ)
    a=np.pi/(2*u);ex=np.exp(a*(s-max(s)))
    vs=-(w*c*ex)@eps/(4*u*((w*ex)@rho))
    es=hubbard_spectrum(1,U,mu);M=(es[-1]-es[0])/2
    return dict(U=U,Q=Q,N=N,n=float(w@rho),mu=float(mu),vc=float(vc),vs=float(vs),A=float(vc+vs),M=float(M),ratio=float((vc+vs)/(2*np.pi*M)),edge_residual=float(abs(ep(Q))))
if __name__=='__main__':
    results=[];maxerr=0
    for U in [.5,1,2,4,8,16]:
        for frac in [.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99]:
            p=point(U,np.pi*frac,128);r=point(U,np.pi*frac,192)
            err=max(abs(r[k]-p[k]) for k in ['vc','vs','n','mu'])
            if err>2e-6:
                p=r;r=point(U,np.pi*frac,288)
                err=max(abs(r[k]-p[k]) for k in ['vc','vs','n','mu'])
            r['grid_change']=err;maxerr=max(maxerr,err)
            assert 0<r['n']<1+1e-9 and r['vc']>=0 and r['vs']>=0
            assert r['A']<=6*np.sin(np.pi*r['n']/2)+1e-6
            results.append(r)
    half=[]
    for U in [1,2,4,8]:
        p=point(U,np.pi,192);exact=2*iv(1,2*np.pi/U)/iv(0,2*np.pi/U)
        half.append(dict(U=U,vs=p['vs'],exact=exact,error=abs(p['vs']-exact)))
    out=dict(points=results,half_filling_checks=half,max_grid_change=maxerr)
    with open('hubbard_quadrature.json','w') as f:json.dump(out,f,indent=2)
    print('points',len(results),'max grid change',maxerr)
    print('maximum sampled ratio',max(results,key=lambda x:x['ratio']))
    print('half filling check',half)
