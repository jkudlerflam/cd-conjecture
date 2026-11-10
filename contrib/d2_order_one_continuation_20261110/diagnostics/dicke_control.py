#!/usr/bin/env python3
"""Exact weighted-Dicke ground states of the D=0 rank-one parent."""
import json
import math
import analyze_saved as A
import numpy as np
from check_mappings import embed,product
Z=A.Z

def weights(L,jz):
    r=2*jz-math.sqrt(4*jz*jz-1)
    return np.array([-r if j%2==0 else 1 for j in range(L)])

def profiles(L,jz):
    w=weights(L,jz);N=L//2;left=[np.ones(1)];right=[None]*(L+1);right[L]=np.ones(1)
    for x in w:
        p=np.convolve(left[-1],[1,x*x]);left.append(p/p.max())
    for j in range(L-1,-1,-1):
        p=np.convolve(right[j+1],[1,w[j]*w[j]]);right[j]=p/p.max()
    cuts=[]
    for cut in range(1,L):
        k=np.arange(max(0,N-(L-cut)),min(cut,N)+1)
        p=left[cut][k]*right[cut][N-k];p/=p.sum()
        ww=[p[(k%2)==q].sum() for q in (0,1)]
        S=Z.entropy(p);binary=Z.entropy(ww)
        cuts.append(dict(cut=cut,S=S,S_binary=binary,S_conditional=S-binary))
    return cuts

def dense_check(L,jz):
    w=weights(L,jz);N=L//2;v=np.zeros(2**L)
    for s in range(2**L):
        bits=[(s>>(L-1-j))&1 for j in range(L)]
        if sum(bits)==N:v[s]=np.prod([w[j] for j,b in enumerate(bits) if b])
    v/=np.linalg.norm(v)
    rot=product([Z.PX if j%2 else Z.I2 for j in range(L)]);v=rot@v
    xx=[embed(Z.PX,j,L) for j in range(L)];yy=[embed(Z.PY,j,L) for j in range(L)]
    zz=[embed(Z.ZZ_,j,L) for j in range(L)];hf=-2*math.sqrt(jz*jz-.25)
    H=sum(.5*xx[j]@xx[(j+1)%L]-.5*yy[j]@yy[(j+1)%L]+jz*zz[j]@zz[(j+1)%L]
        for j in range(L))+hf*sum(zz)
    err=float(np.linalg.norm(H@v+L*jz*v));assert err<1e-12
    wrap=.5*xx[-1]@xx[0]-.5*yy[-1]@yy[0]+jz*zz[-1]@zz[0]+hf/2*(zz[-1]+zz[0])
    Hopen=H-wrap
    openerr=float(np.linalg.norm(Hopen@v+(L-1)*jz*v));assert openerr<1e-12
    ps=np.linalg.svd(v.reshape(2**(L//2),-1),compute_uv=False)**2
    entropy=Z.entropy(ps);poly=profiles(L,jz)[L//2-1]['S'];assert abs(entropy-poly)<1e-12
    return dict(L=L,Jz=jz,eigen_residual=err,open_half_edge_field_residual=openerr,
                dense_S=entropy,polynomial_S=poly)

def main():
    out=dict(checks=[],rows=[],slopes=[])
    for jz in [1.,5/8]:
        out['checks'] += [dense_check(L,jz) for L in [4,6,8]]
        prev=None
        for L in [24,48,96,192]:
            cuts=profiles(L,jz)
            stats={f:A.D.entropy_stats([c[f] for c in cuts],L) for f in ['S','S_binary','S_conditional']}
            row=dict(Jz=jz,L=L,D=0,hz=-2*math.sqrt(jz*jz-.25),N_transformed=L//2,cuts=cuts,stats=stats)
            out['rows'].append(row)
            if prev:
                slope=dict(Jz=jz,sizes=[prev['L'],L],c={})
                for f in stats:
                    slope['c'][f]={k:6*(stats[f][k]-prev['stats'][f][k])/math.log(L/prev['L']) for k in ['S_central3','S_fractional']}
                out['slopes'].append(slope);print(slope,flush=True)
            prev=row
    (A.HERE/'dicke_control.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
