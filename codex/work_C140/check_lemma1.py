"""Checks of the dual identity sigma=r0+J*(sigma dx-mu) and of Lemma 1 (two-neighbour bound). Floating point. J=1."""
import json
import numpy as np
from scipy.integrate import quad
from bethe_scan import solve
out=[]
for delta in [0.1,0.5,0.9,0.99]:
    g=np.arccos(delta); b=np.pi/g
    a1=lambda x: np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
    K=lambda x: np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
    r0=lambda x: 1/(2*g*np.cosh(b*x))
    Kh=lambda w: np.exp(-g*w)*(-np.expm1(-(np.pi-2*g)*w))/(-np.expm1(-np.pi*w))
    Jh=lambda w: Kh(w)/(1+Kh(w))
    # J on a grid by inverse cosine transform
    W=np.linspace(1e-9,400/g,400001); dw=W[1]-W[0]; JW=Jh(W)
    xg=np.linspace(0,60*g,6001)
    Jg=np.array([np.trapezoid(JW*np.cos(W*xx),W) for xx in xg])/np.pi
    Jf=lambda x: np.interp(np.abs(x),xg,Jg,right=0.)
    J0=Jg[0]; J0b=(np.pi-2*g)/((np.pi-g)*2*g)
    for L,n in [(16,8),(17,8),(64,32),(64,20),(256,128),(256,127),(256,64)]:
        E,d,x=solve(L,n,delta,return_roots=True)
        sig=lambda y: a1(y)-K(np.asarray(y)[:,None]-x[None,:]).sum(1)/L
        yy=np.linspace(-x[-1]-40*g,x[-1]+40*g,40001); dy=yy[1]-yy[0]; sy=sig(yy)
        pts=np.linspace(0,x[-1]+3*g,301)
        idres=0.; worst=np.inf; worst_bulk=np.inf
        for p in pts:
            conv=np.sum(Jf(p-yy)*sy)*dy
            rhs=r0(p)+conv-Jf(p-x).sum()/L
            sp=sig(np.array([p]))[0]
            idres=max(idres,abs(sp-rhs)/r0(p))
            k=np.searchsorted(x,p)
            nb=[x[i] for i in (k-1,k) if 0<=i<n]
            outm=np.abs(yy)>x[-1]
            lb=r0(p)+np.sum(Jf(p-yy[outm])*sy[outm])*dy-max(Jf(p-np.array(nb)))/L
            worst=min(worst,sp-lb)
            worst_bulk=min(worst_bulk,sp-(r0(p)-J0/L))
        row=dict(delta=delta,L=L,n=n,J0=float(J0),J0_bound=float(J0b),identity_rel_err=float(idres),
                 min_sig_minus_lemma1=float(worst),min_sig_minus_bulk=float(worst_bulk))
        out.append(row); print(row,flush=True)
json.dump(out,open('check_lemma1.json','w'),indent=0)
