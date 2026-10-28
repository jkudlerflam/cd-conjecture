"""C35 finite algebra and formula checks. Numpy/scipy, no external imports."""
import json
import numpy as np
from scipy.special import gammaln
from scipy.optimize import minimize_scalar

out={}
H=np.cumsum(1/np.arange(1,2001,dtype=float))
rows=[]
best=(-1,None)
second=(-1,None)
for N in range(2,2001):
    k=np.arange(1,2001)
    ld=gammaln(N+k)-gammaln(N)-gammaln(k+1)
    q=2*k*(N*N-1)/(N*(N+k)*H)/(ld/np.log(2))
    i=int(np.argmax(q))
    if q[i]>best[0]: best=(float(q[i]),[N,int(k[i])])
    q2=q.copy()
    if N==2:q2[0]=-1
    i=int(np.argmax(q2))
    if q2[i]>second[0]:second=(float(q2[i]),[N,int(k[i])])
    assert np.max(q)<1+1e-11
out['symmetric_scan']={'N':[2,2000],'k':[1,2000], 'pairs':1999*2000,'max_q':best,'max_except_heisenberg':second}
for N,k in [(2,1),(3,1),(2,2),(7,2),(8,2),(10,10),(100,100),(1000,1000)]:
    ld=float(gammaln(N+k)-gammaln(N)-gammaln(k+1))
    R=k*(N*N-1)/(N*(N+k)*H[k-1])
    rows.append({'N':N,'k':k,'log_d':ld,'A_over_2piM':float(R),'q':float(2*R*np.log(2)/ld)})
out['examples']=rows
# Coefficient in R ~ F(alpha) log d / log log d, k=alpha N.
def coef(x):
    a=np.exp(x)
    f=np.log1p(a)+a*np.log1p(1/a)
    return a/((1+a)*f)
opt=minimize_scalar(lambda x:-coef(x),bounds=(-20,20),method='bounded')
out['asymptotic_ray_optimization']={'log_alpha':float(opt.x),'coefficient_R':float(-opt.fun),'at_alpha_1':float(coef(0))}
# SO/Sp local matrices, and the swap projection that removes arbitrary gauges.
rng=np.random.default_rng(3501)
loc=[]
for group,ds in [('SO',range(3,9)),('Sp',[2,4,6,8])]:
    for d in ds:
        P=np.eye(d*d).reshape(d,d,d,d).transpose(1,0,2,3).reshape(d*d,d*d)
        if group=='SO':
            vec=np.eye(d).ravel(); h=P-2*np.outer(vec,vec)/(d-2); expected=d/(d-2)
        else:
            n=d//2; J=np.zeros((d,d));J[:n,n:]=np.eye(n);J[n:,:n]=-np.eye(n)
            vec=J.ravel();h=P+np.outer(vec,vec)/(n+1);expected=1 if n>=2 else .5
        ev=np.linalg.eigvalsh(h); M=(ev[-1]-ev[0])/2
        o=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));o=(o+o.conj().T)/2
        hg=h+np.kron(o,np.eye(d))-np.kron(np.eye(d),o)
        err=np.max(np.abs((hg+P@hg@P)/2-h))
        assert abs(M-expected)<1e-12 and err<1e-12
        loc.append({'group':group,'d':d,'M':float(M),'distinct_eigenvalues':np.unique(np.round(ev,10)).tolist(),'swap_average_error':float(err)})
out['local_checks']=loc
# Rectangular Weyl dimensions: recurrence/unimodality and d>=d_sym.
worst=1e99
for N in range(2,61):
    for k in range(1,31):
        dims=[]
        for r in range(1,N):
            ld=sum(np.log1p(k/(j-i)) for i in range(1,r+1) for j in range(r+1,N+1))
            dims.append(ld)
        worst=min(worst,min(dims)-dims[0])
        assert min(dims)>dims[0]-1e-10
out['rectangular_dimension_check']={'N':[2,60],'k':[1,30],'all_heights':True,'minimum_log_dimension_difference':float(worst)}
with open('checks.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps({k:v for k,v in out.items() if k not in ['local_checks']},indent=2))
