"""Exact finite-grid XX free-fermion formulas, copied without change from
codex/work_C123/checks.py (functions xx_parts, xx_R, contour, xx_kernels,
xx_histories, xx_case). See that folder for the derivation and calibration."""
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
