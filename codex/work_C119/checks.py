"""C119 small exact-trace and analytic-formula checks. NumPy/SciPy only.
Larger scans are guarded by SLURM_JOB_ID. All output stays beside this file.
"""
from pathlib import Path
import argparse
import json
import os
import numpy as np
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parent
PI = np.longdouble(str(np.pi))

def projected(eps, beta, parity):
    eps = np.asarray(eps, dtype=np.longdouble).ravel()
    x = beta * np.abs(eps)
    base = -beta * np.minimum(eps, 0).sum() + np.log1p(np.exp(-x)).sum()
    sign = (-1)**(int(np.count_nonzero(eps < 0)) + parity)
    with np.errstate(divide='ignore'):
        lt = np.log(np.tanh(x/2)).sum()
    correction = np.log1p(np.exp(lt)) if sign > 0 else np.log(-np.expm1(lt))
    return base + correction - np.log(np.longdouble(2))

def logz(n, beta, bands):
    k0 = 2*PI*np.arange(n, dtype=np.longdouble)/n
    return np.logaddexp(projected(bands(k0+PI/n), beta, 0),
                        projected(bands(k0), beta, 1))

def ratio(n, tau, bands):
    beta = np.longdouble(tau)*n
    return float(2*logz(n, beta, bands)-logz(2*n, beta, bands))

def xx(k):
    return np.cos(k)

def two_band(eps, a, t=1):
    def bands(k):
        e = np.sqrt((2*t*np.sin(k))**2+eps**2)
        return np.array([-a*eps+e, -a*eps-e])
    return bands

def theta_log(s):
    if s < 1:
        return -0.5*np.log(s)+theta_log(1/s)
    n = np.arange(1, max(3, int(np.sqrt(45/(np.pi*s)))+2))
    return np.log1p(2*np.exp(-np.pi*s*n*n).sum())

def eta_log(y):
    if y < 1:
        return eta_log(1/y)-0.5*np.log(y)
    n = np.arange(1, max(3, int(45/(2*np.pi*y))+2))
    return -np.pi*y/12 + np.log1p(-np.exp(-2*np.pi*y*n)).sum()

def boson_f(y, K=1):
    return theta_log(y/(2*K))+theta_log(2*K*y)-2*eta_log(y)

def boson_r(y, K=1):
    # Modular channel avoids subtracting the extensive thermal term.
    if y < 0.3:
        def excitation(z):
            return theta_log(z/(2*K))+theta_log(2*K*z)-2*(eta_log(z)+np.pi*z/12)
        return 2*excitation(1/y)-excitation(2/y)
    return 2*boson_f(y,K)-boson_f(y/2,K)

def xx_f2(y):
    nmax = max(20, int(60/(2*np.pi*y))+2)
    r = np.arange(nmax)+0.5
    q = np.exp(-2*np.pi*y*r)
    logs, ds = [], []
    for sign in [1, -1]:
        logs.append(np.pi*y/6+4*np.log1p(sign*q).sum()-np.log(2))
        ds.append(7*np.pi**3*y/360+16*np.pi**3*y/3*np.sum(r**3*sign*q/(1+sign*q)))
    r = np.arange(1,nmax+1)
    q = np.exp(-2*np.pi*y*r)
    logs.append(-np.pi*y/3+4*np.log1p(q).sum()+np.log(2))
    ds.append(-np.pi**3*y/45+16*np.pi**3*y/3*np.sum(r**3*q/(1+q)))
    weights = np.exp(np.array(logs)-logsumexp(logs))
    assert abs(logsumexp(logs)-boson_f(y)) < 2e-12
    return float(np.dot(weights, ds))

def common_parity_f(y):
    nmax=max(20,int(60/(2*np.pi*y))+2)
    r=np.arange(nmax)+.5
    q=np.exp(-2*np.pi*y*r)
    logs=[np.pi*y/3+8*np.log1p(s*q).sum()-np.log(2) for s in [1,-1]]
    r=np.arange(1,nmax+1)
    logs.append(-2*np.pi*y/3+8*np.log1p(np.exp(-2*np.pi*y*r)).sum()+np.log(8))
    return float(logsumexp(logs))

def c106_f2(tau,eps,a):
    k=np.arcsin(eps*np.sqrt(a*a-1)/2)
    e=a*eps
    d1=4*np.sin(2*k)
    d2=8*np.cos(2*k)
    d3=-16*np.sin(2*k)
    v=d1/(2*e)
    kap=d2/(2*e)-d1*d1/(4*e**3)
    lam=d3/(2*e)-3*d1*d2/(4*e**3)+3*d1**3/(8*e**5)
    y=tau*v
    nmax=max(30,int(60/(2*np.pi*y))+2)
    logs,ds=[],[]
    for shift,sign,B4,weight in [(.5,1,7/240,.5),(.5,-1,7/240,.5),(1,1,-1/30,8)]:
        r=np.arange(nmax)+shift
        q=np.exp(-2*np.pi*y*r)
        w=sign*q/(1+sign*q)
        vacuum=np.pi*y/3 if shift==.5 else -2*np.pi*y/3
        logs.append(vacuum+8*np.log1p(sign*q).sum()+np.log(weight))
        ds.append(-4*tau*np.pi**3*lam*B4/3
                  -4*tau*lam*(2*np.pi)**3/3*np.sum(r**3*w)
                  +tau*tau*kap*kap*(2*np.pi)**4*np.sum(r**4*w*(1-w)))
    return float(np.dot(np.exp(np.array(logs)-logsumexp(logs)),ds)),v

def fermions(m):
    out=[]
    for j in range(m):
        c=np.zeros((2**m,2**m),complex)
        for state in range(2**m):
            if state & (1<<j):
                c[state^(1<<j),state]=(-1)**((state & ((1<<j)-1)).bit_count())
        out.append(c)
    return out

def c106_bond(eps,a,t=1):
    c=fermions(4)
    B=np.diag([-1j*t,1j*t])
    X=np.array([[-a*eps,eps],[eps,-a*eps]])
    h=np.zeros((16,16),complex)
    for u in range(2):
        for v in range(2):
            term=B[u,v]*c[u].conj().T@c[2+v]
            h+=term+term.conj().T
            h+=X[u,v]*(c[u].conj().T@c[v]+c[2+u].conj().T@c[2+v])/2
    return h

def spin_ring(h,d,n):
    dim=d**n
    mat=np.zeros((dim,dim),complex)
    for state in range(dim):
        for i in range(n):
            j=(i+1)%n
            ai,aj=(state//d**i)%d,(state//d**j)%d
            for bi in range(d):
                for bj in range(d):
                    # Local mode ordering uses the first site as low-order digit.
                    mat[state+(bi-ai)*d**i+(bj-aj)*d**j,state]+=h[bi+d*bj,ai+d*aj]
    return mat

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cluster',action='store_true')
    args=ap.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Large scan requires a SLURM allocation.')
    out={}
    ed=[]
    for n in [3,4]:
        ev=np.linalg.eigvalsh(spin_ring(c106_bond(.17,1.2),4,n))
        for beta in [0, .7, 5]:
            err=abs(logsumexp(-beta*ev)-logz(n,beta,two_band(.17,1.2)))
            assert err<2e-11, (n,beta,err)
            ed.append(float(err))
    out['C106_ED_max_logZ_error']=max(ed)
    rows=[]
    for y in [.5,1,2,4]:
        coeff=2*xx_f2(y)-xx_f2(y/2)/4
        for n in [32,64,128,256]:
            r=ratio(n,y,xx)
            residual=r-boson_r(y)-coeff/n**2
            rows.append(dict(L=n,y=y,R=r,r=boson_r(y),coefficient=coeff,
                             residual_times_L4=residual*n**4))
        assert abs(rows[-1]['residual_times_L4']) < max(300,abs(rows[-4]['residual_times_L4'])*2)
    out['XX_subleading']=rows
    out['slow_sectors']=[dict(y=y,r=boson_r(y),leading=4*np.exp(-np.pi/(2*y)))
                         for y in [.5,.25,.125,.0625,.03125]]
    out['own_scale']=[dict(v=v,tau=1/v,r=boson_r(1)) for v in [1,.1,.001]]
    out['stacks']=[dict(epsilon=e,R=boson_r(1/(1+e))+boson_r(e/(1+e)))
                   for e in [.5,.1,.02,.005]]
    out['C106_crossover']=[]
    for n in [32,128,512,2048]:
        for eps in [0,1e-5,.001,.01,.1]:
            a=1.2
            M=np.sqrt(4+eps**2)
            r=ratio(n,1/M,two_band(eps,a))
            out['C106_crossover'].append(dict(L=n,epsilon=eps,u=n*eps/2,R=r,
                bound=3*np.log(2)+10*np.pi/3))
            assert r <= 3*np.log(2)+10*np.pi/3+1e-8
    out['C106_subleading']=[]
    a=1.2
    eps=2*np.sin(np.pi/8)/np.sqrt(a*a-1)
    for tau in [1,2]:
        f2,v=c106_f2(tau,eps,a)
        f2half,_=c106_f2(tau/2,eps,a)
        coeff=2*f2-f2half/4
        rinf=2*common_parity_f(tau*v)-common_parity_f(tau*v/2)
        for n in [64,128,256,512]:
            r=ratio(n,tau,two_band(eps,a))
            out['C106_subleading'].append(dict(L=n,tau=tau,R=r,r=rinf,
                coefficient=coeff,residual_times_L4=(r-rinf-coeff/n**2)*n**4))
        # The measured L^-2 coefficient must converge to the derived one.
        measured=(out['C106_subleading'][-1]['R']-rinf)*512**2
        assert abs(measured-coeff)<.01*max(1,abs(coeff)), (measured,coeff)
    out['massive']=[]
    for n in [32,128,512]:
        for gap in [.001,.01,.1]:
            r=ratio(n,.5,two_band(gap,0))
            out['massive'].append(dict(L=n,gap=gap,z=n*gap/2,R=r))
    out['radius']=[]
    for K in [1,10,100,1000]:
        out['radius'].append(dict(K=K,r_fixed_y=boson_r(1,K)))
    out['XXZ_ferromagnetic']=[]
    for eta in [.1,.01,.001]:
        K=np.pi/(2*eta)
        v=np.pi*np.sin(eta)/(2*(np.pi-eta))
        for tau in [1,4,16]:
            limit=2*theta_log(2/(np.pi*tau))-theta_log(4/(np.pi*tau))
            out['XXZ_ferromagnetic'].append(dict(eta=eta,tau=tau,
                r=boson_r(v*tau,K),limit=limit))
    ys=np.geomspace(.02,100,200)
    out['XX_scaling_ceiling_max_excess']=max(boson_r(y)-np.pi*y/4 for y in ys)
    assert out['XX_scaling_ceiling_max_excess'] < 1e-10
    # Adversarial but deterministic free-family scan. Extended mode remains linear algebra free.
    ns=[4,8,16,32,64] if not args.cluster else [4,8,16,32,64,128,256,512,1024,2048]
    es=np.geomspace(1e-5,.5,8 if not args.cluster else 50)
    best=(-np.inf,None)
    for n in ns:
        for eps in es:
            for a in [1.001,1.2,2]:
                M=np.sqrt(4+eps**2)
                for s in [.1,1,4]:
                    r=ratio(n,s/M,two_band(eps,a))
                    excess=r-3*np.log(2)-10*np.pi*s/3
                    if excess>best[0]: best=(excess,[n,float(eps),a,s,r])
                    assert excess < 1e-9
    out['uniform_scan_largest_excess']=best
    name='checks_cluster.json' if args.cluster else 'checks_small.json'
    (ROOT/name).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(file=name,ED_error=out['C106_ED_max_logZ_error'],
        largest_bound_excess=best,XX_ceiling_excess=out['XX_scaling_ceiling_max_excess']),indent=2))

if __name__=='__main__':
    main()
