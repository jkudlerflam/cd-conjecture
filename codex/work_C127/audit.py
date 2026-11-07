"""Momentum-sector comparisons and a compact summary of the C127 checks."""
from pathlib import Path
import json
import numpy as np
from scipy.special import logsumexp
from checks import boson_log
ROOT=Path(__file__).resolve().parent

def coefficients(m,N=80):
    a=np.zeros((N+1,m));a[0,0]=1
    for k in range(1,N+1):
        for sign in [-1,1]:
            for species in range(3):
                for n in range(k,N+1):a[n]+=np.roll(a[n-k],sign*k)
    return a

def folded_boson(t,m,coef):
    if t<.5:
        logs=np.full_like(coef,-np.inf);np.log(coef,out=logs,where=coef>0)
        return logsumexp(logs-np.arange(len(coef))[:,None]/t,axis=0)
    cutoff=int(np.ceil(40*t));j=np.arange(1,cutoff+1)
    theta=2*np.pi*np.arange(m)/m
    logcharacters=-6*np.log(abs(1-np.exp(-j[None,:]/t+1j*theta[:,None]*j[None,:]))).sum(axis=1)
    largest=max(logcharacters)
    sectors=np.fft.fft(np.exp(logcharacters-largest)).real/m
    assert min(sectors)>0
    return np.log(sectors)+largest

def main():
    raw=json.loads((ROOT/'checks.json').read_text());out=[]
    coefs={m:coefficients(m) for m in range(1,5)}
    K=np.log(3)+8/(3*np.e)
    for row in raw['transfer']:
        V=4*np.e*row['M'];m=row['m']
        for c in row['comparisons']:
            L=c['L'];t=V*row['beta']/(2*np.pi*L)
            B=folded_boson(t,m,coefs[m])
            total=logsumexp(B);target,_=boson_log(t,6)
            assert abs(total-target)<1e-9
            for sector in row['sectors']:
                if not sector['gaps']:continue
                z=logsumexp(-L*np.array(sector['gaps']))
                out.append(dict(delta=row['delta'],beta=row['beta'],m=m,L=L,k=sector['k'],
                    logZ=float(z),reference_logZ=float(B[sector['k']]),Kstar_sector_margin=float(K+B[sector['k']]-z)))
    (ROOT/'momentum_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    worst=min(out,key=lambda r:r['Kstar_sector_margin'])
    print(json.dumps(dict(worst_folded_momentum_comparison=worst,
        sector_comparisons=len(out),failed_sector_comparisons=sum(r['Kstar_sector_margin']<-1e-8 for r in out)),indent=2))
if __name__=='__main__':main()
