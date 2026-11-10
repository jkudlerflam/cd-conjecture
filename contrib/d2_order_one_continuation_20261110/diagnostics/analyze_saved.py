#!/usr/bin/env python3
"""Postprocess frozen MPS checkpoints. No DMRG updates and no source writes."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=ROOT/'contrib/d2_order_one_20261110'
sys.path.insert(0,str(OLD))
import dmrg_screen as D
Z=D.Z
NAMES=['bridge_r4_u1_k4','bridge_r3_u1_k4','bridge_r1.5_u1.5_k3',
       'bridge_r2_u1.5_k4','control_xx','control_two_xx']


def schmidt_diagnostics(psi):
    psi=psi.copy();L=psi.L
    cuts=[]
    for j in range(L-1):
        Z._qr_move_right(psi,j)
        C=psi.M[j+1];Dl,_,Dr=C.shape
        n0l,n0r=psi.n0[j+1],psi.n0[j+2]
        pp=[]
        for q in (0,1):
            rows=Z.bsl(n0l,Dl,q)
            A=np.hstack([C[rows,0,Z.bsl(n0r,Dr,q)],C[rows,1,Z.bsl(n0r,Dr,1-q)]])
            p=np.linalg.svd(A,compute_uv=False)**2 if A.size else np.array([])
            pp.append(p)
        norm=sum(p.sum() for p in pp);pp=[p/norm for p in pp]
        weights=np.array([p.sum() for p in pp])
        total=Z.entropy(np.concatenate(pp))
        binary=Z.entropy(weights)
        conditional=sum(float(w)*Z.entropy(p) for w,p in zip(weights,pp) if w>1e-30)
        row=dict(cut=j+1,weights=weights.tolist(),S=total,S_binary=binary,
                 S_conditional=conditional,identity_error=total-binary-conditional)
        full=np.concatenate(pp)
        row['S2']=float(-np.log(np.sum(full**2)))
        row['S3']=float(-.5*np.log(np.sum(full**3)))
        if j+1==L//2:
            row['probabilities_by_parity']=[p.tolist() for p in pp]
            n=min(len(pp[0]),len(pp[1]),16)
            row['leading_pair_relative_differences']=(np.abs(pp[0][:n]-pp[1][:n])/
                np.maximum(pp[0][:n]+pp[1][:n],1e-30)).tolist()
        cuts.append(row)
    return cuts


def envelope(r):
    L=r['L'];corr=r['measurements']['long_connected_correlations']
    out={}
    for axis in ['X','Y','Z']:
        pairs=np.array(corr[axis]);dist,val=pairs.T
        axisout={}
        for label,lo,hi in [('near',L/8,L/4),('far',L/4,L/2)]:
            v=val[(dist>=lo)&(dist<=hi)]
            axisout[label]=dict(rms=float(np.sqrt(np.mean(v*v))),
                mean_abs=float(np.mean(np.abs(v))),max_abs=float(np.max(np.abs(v))))
        out[axis]=axisout
    return out


def density_diagnostics(r):
    n=np.array(r['measurements']['dens']);L=len(n)
    y=n[L//4:3*L//4];x=np.arange(len(y))
    trend=np.polyval(np.polyfit(x,y,1),x);osc=y-trend
    q=np.linspace(0,np.pi,2001)
    # Positive-frequency spectrum normalized by the central-window length.
    power=np.abs(np.exp(-1j*np.outer(q,x))@osc)**2/len(x)**2
    k=int(np.argmax(power[1:]))+1
    return dict(mean=float(y.mean()),range=float(np.ptp(y)),rms_detrended=float(np.sqrt(np.mean(osc**2))),
                linear_slope=float(np.polyfit(x,y,1)[0]),q_peak=float(q[k]),peak_amplitude=float(2*np.sqrt(power[k])))


def main():
    records=[]
    for name in NAMES:
        bykey={}
        for p in (OLD/'dmrg').glob(name+'_L*.json'):
            r=json.loads(p.read_text())
            if r['status']!='done':continue
            key=(r['L'],r['parity'])
            if key not in bykey or r['E_physical']<bykey[key][1]['E_physical']:bykey[key]=(p,r)
        for (L,P),(path,r) in sorted(bykey.items()):
            start=time.time();ckpt=Path(r['checkpoint'])
            with np.load(ckpt,allow_pickle=False) as z:
                meta=json.loads(str(z['meta']))
            assert meta['couplings']==r['couplings'] and len(meta['history'])==len(r['history'])
            psi=Z.MPS.load(ckpt)
            cuts=schmidt_diagnostics(psi)
            total=np.array([c['S'] for c in cuts])
            assert np.max(np.abs(total-np.array(r['Sprof'])))<1e-9
            summaries={}
            for field in ['S','S_binary','S_conditional']:
                summaries[field]=D.entropy_stats([c[field] for c in cuts],L)
            rec=dict(name=name,L=L,parity=P,seed=r['seed'],E=r['E_physical'],
                source_record=str(path.relative_to(ROOT)),record_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                checkpoint=str(ckpt),checkpoint_sha256=hashlib.sha256(ckpt.read_bytes()).hexdigest(),
                cuts=cuts,summaries=summaries,envelope=envelope(r),density=density_diagnostics(r),
                seconds=time.time()-start)
            records.append(rec)
            print(name,L,P,'binary',summaries['S_binary']['S_fractional'],
                  'conditional',summaries['S_conditional']['S_fractional'],flush=True)
            (HERE/'saved_diagnostics.json').write_text(json.dumps(records,indent=2)+'\n')
    pairs=[]
    for name in NAMES:
        rs=[r for r in records if r['name']==name]
        ls=sorted(set(r['L'] for r in rs))
        gs=[min([r for r in rs if r['L']==L],key=lambda x:x['E']) for L in ls]
        for a,b in zip(gs,gs[1:]):
            out=dict(name=name,sizes=[a['L'],b['L']],parities=[a['parity'],b['parity']],slopes={})
            for field in ['S','S_binary','S_conditional']:
                out['slopes'][field]={key:6*(b['summaries'][field][key]-a['summaries'][field][key])/math.log(b['L']/a['L'])
                    for key in ['S_central3','S_fractional']}
            pairs.append(out)
    (HERE/'entropy_decomposition_slopes.json').write_text(json.dumps(pairs,indent=2)+'\n')
    print(json.dumps(pairs,indent=2))


if __name__=='__main__':main()
