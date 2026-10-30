"""C123 exact finite-grid XX determinants. All output remains in this folder."""
from pathlib import Path
from itertools import product
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp

ROOT=Path(__file__).resolve().parent

def xx_parts(n,beta):
    # Inherited C121 full physical trace, evaluated without bulk subtraction.
    pieces=[]
    for shift,wanted in [(0.5,0),(0.,1)]:
        eps=np.cos(2*np.pi*(np.arange(n)+shift)/n)
        eps[np.abs(eps)<1e-13]=0.
        parity=int(np.count_nonzero(eps<0)%2)
        even,odd=0.,-np.inf
        for energy in np.abs(eps):
            x=-beta*energy
            even,odd=np.logaddexp(even,odd+x),np.logaddexp(odd,even+x)
        lower=(shift==.5) if n%4==0 else (shift==0.)
        offset=0. if lower else np.tan(np.pi/(2*n))
        pieces.append((even,odd)[wanted^parity]-beta*offset)
    return float(np.logaddexp(*pieces))

def xx_R(n,beta):
    x=np.pi/(2*n)
    defect=2*np.sin(x/2)**2/(np.sin(x)*np.cos(x))
    return float(beta*defect+2*xx_parts(n,beta)-xx_parts(2*n,beta))

def contour(n,beta,count,sites):
    times=np.repeat(np.arange(2*count)*beta/(2*count),len(sites))
    positions=np.tile(sites,2*count)
    return times,positions

def xx_kernels(n,beta,count,sites):
    times,positions=contour(n,beta,count,sites)
    size=len(times)
    pieces=[]
    for shift,wanted in [(0.5,0),(0.,1)]:
        mom=2*np.pi*(np.arange(n)+shift)/n
        eps=np.cos(mom)
        eps[np.abs(eps)<1e-13]=0.
        neg=eps<0; pos=eps>0; zero=eps==0; nz=~zero
        lower=(shift==.5) if n%4==0 else (shift==0.)
        offset=0. if lower else np.tan(np.pi/(2*n))
        boltz=np.exp(-beta*np.abs(eps))
        weights=[float(np.sum(np.log1p(boltz))-beta*offset),
                 float(np.sum(np.log(-np.expm1(-beta*np.abs(eps[nz])))))-beta*offset]
        matrices=[np.zeros((size,size),complex),np.zeros((size,size),complex)]
        for i in range(size):
            for j in range(size):
                dt=times[i]-times[j]
                # Equal-time distinct sites anticommute, so either ordering agrees.
                forward=dt>=0
                decay=np.where(pos, (beta-dt if forward else -dt),
                               (dt if forward else beta+dt))
                mag=np.exp(-np.abs(eps)*decay)
                phase=np.exp(1j*mom*(positions[j]-positions[i]))/n
                val=mag/(1+boltz)
                if not forward: val=-val
                matrices[0][i,j]=np.dot(phase,val)
                signed=np.where(pos,-1.,1.)
                valp=np.zeros(n)
                valp[nz]=signed[nz]*mag[nz]/(-np.expm1(-beta*np.abs(eps[nz])))
                # Remove the divergent zero-mode occupation term, retain -theta.
                if not forward: valp[zero]=-1.
                matrices[1][i,j]=np.dot(phase,valp)
        u=np.exp(-1j*positions[:,None]*mom[zero][None,:])/np.sqrt(n)
        v=u.conj().T
        g=int(np.count_nonzero(zero))
        bordered=np.block([[matrices[1],u],[-v,np.zeros((g,g))]])
        sign=(-1.)**(wanted+int(np.count_nonzero(neg))+g)
        pieces.append((matrices[0],bordered,weights,sign))
    return pieces

def xx_histories(n,beta,count,sites):
    pieces=xx_kernels(n,beta,count,sites)
    size=2*count*len(sites)
    vals=[]
    for flat in product(range(2),repeat=count*len(sites)):
        labels=np.array(flat).reshape(count,len(sites))
        empty=1-np.concatenate([labels,labels[::-1]]).ravel()
        total=0j
        for ordinary,bordered,weights,sign in pieces:
            a=ordinary.copy(); b=bordered.copy()
            idx=np.arange(size)
            a[idx,idx]-=empty; b[idx,idx]-=empty
            # There are an even number of empty projectors on this contour.
            total+=.5*(np.exp(weights[0])*np.linalg.det(a)
                      +sign*np.exp(weights[1])*np.linalg.det(b))
        assert abs(total.imag)<1e-9*max(1,abs(total.real)),total
        assert total.real>-1e-12,total
        vals.append(float(total.real))
    return np.array(vals)

def xx_case(n,a,count):
    beta=2*a*n*count
    q=xx_histories(n,beta,count,[0])
    qd=xx_histories(2*n,beta,count,[0,n])
    x=np.pi/(2*n)
    defect=2*np.sin(x/2)**2/(np.sin(x)*np.cos(x))
    g=beta*defect+np.log(q@q)-np.log(qd.sum())
    r=xx_R(n,beta)
    loss=xx_parts(n,beta)-np.log(q.sum())
    lossd=xx_parts(2*n,beta)-np.log(qd.sum())
    assert -1e-9<=loss<=count*np.log(2)+1e-9
    assert -1e-9<=lossd<=2*count*np.log(2)+1e-9
    assert r-3*count*np.log(2)-1e-9<=g<=r+2*count*np.log(2)+1e-9
    return dict(L=n,a=a,N=count,G=float(g),G_per_N=float(g/count),R=r,
                loss=float(loss),loss_D=float(lossd))

def dense_xx(n):
    dim=2**n
    h=np.zeros((dim,dim))
    for s in range(dim):
        for i in range(n):
            j=(i+1)%n
            if ((s>>i)&1)!=((s>>j)&1): h[s^(1<<i)^(1<<j),s]+=.5
    return eigh(h)

def dense_histories(ev,vec,beta,count,sites):
    dim=len(ev)
    f=(vec*np.exp(-beta*(ev-ev[0])/(4*count)))@vec.T
    e=f@f
    vals=[]
    masks=[(np.arange(dim)>>site)&1 for site in sites]
    for flat in product(range(2),repeat=count*len(sites)):
        labels=np.array(flat).reshape(count,len(sites))
        w=f.copy()
        for j in range(count):
            keep=np.ones(dim,dtype=bool)
            for k in range(len(sites)): keep &= masks[k]==labels[j,k]
            w=keep[:,None]*w
            w=(f if j==count-1 else e)@w
        vals.append(float(np.sum(abs(w)**2)))
    return np.array(vals)

def calibrate():
    rows=[]
    for n in [4,6,8]:
        ev,vec=dense_xx(n)
        for beta in [.2,2.,8.]:
            for count in [1,2,3]:
                for sites in [[0],[0,n//2]]:
                    if n==8 and count==3 and len(sites)==2: continue
                    exact=dense_histories(ev,vec,beta,count,sites)
                    gauss=xx_histories(n,beta,count,sites)
                    err=float(np.max(abs(exact-gauss))/max(exact))
                    assert err<2e-9,(n,beta,count,sites,err,exact,gauss)
                    rows.append(dict(L=n,beta=beta,N=count,sites=sites,error=err))
    return rows

def main():
    start=time.time()
    cal=calibrate()
    print('calibration',len(cal),max(x['error'] for x in cal),flush=True)
    rows=[]
    for n in [4,8,32,128,512,2048,8192]:
        for a in [.125,.25,.5,1.,2.,4.,16.]:
            for count in [1,2,4]:
                rows.append(xx_case(n,a,count))
        print('completed L',n,flush=True)
    for n in [32,512,8192]:
        for a in [.25,1.,4.]: rows.append(xx_case(n,a,6))
    out=dict(calibration=cal,xx=rows,seconds=time.time()-start)
    (ROOT/'xx_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(seconds=out['seconds'],cases=len(rows),
                         large=[r for r in rows if r['L']==8192]),indent=2))

if __name__=='__main__': main()
