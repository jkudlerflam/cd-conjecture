#!/usr/bin/env python3
"""Incremental Renyi postprocessing of immutable completed snapshots."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'d2_order_one_20261110'
sys.path.insert(0,str(OLD))
import dmrg_screen as D
Z=D.Z

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def stats(values,L,alpha):
    out=D.entropy_stats(values,L)
    for fit in out['profiles'].values():fit['c']*=2/(1+1/alpha)
    return out

def spectra(psi):
    psi=psi.copy(); result={str(a):[] for a in (1,2,3)}
    for j in range(psi.L-1):
        Z._qr_move_right(psi,j)
        C=psi.M[j+1];Dl,_,Dr=C.shape
        n0l,n0r=psi.n0[j+1],psi.n0[j+2]
        blocks=[]
        for q in (0,1):
            ix=Z.bsl(n0l,Dl,q)
            A=np.hstack([C[ix,0,Z.bsl(n0r,Dr,q)],C[ix,1,Z.bsl(n0r,Dr,1-q)]])
            if A.size: blocks.append(np.linalg.svd(A,compute_uv=False)**2)
        p=np.concatenate(blocks);p=p/p.sum()
        result['1'].append(float(Z.entropy(p)))
        for a in (2,3):result[str(a)].append(float(np.log(np.sum(p**a))/(1-a)))
    return result

def main():
    output=HERE/'renyi.json'
    cache=json.loads(output.read_text()) if output.exists() else []
    paths=list((OLD/'dmrg').glob('*.json'))
    for folder in ('validation','attractive_scan/results','bridge_scan/dmrg','independent_seed','independent_sizes','simple_attractive/results'):
        paths+=list((HERE/folder).glob('*.json'))
    wanted={'bridge_r4_u1_k4','control_xx','control_two_xx','pref_ref','pref_rat','attr_ref','pref_zrat','b_r5_2_d1_10_q0_hm7_8','simple_attr_7_m1_5_m7'}
    for path in sorted(paths):
        raw=path.read_bytes();r=json.loads(raw)
        if r.get('status')!='done' or r.get('name') not in wanted:continue
        record_hash=hashlib.sha256(raw).hexdigest()
        if any(c['source_record']==str(path) and c['record_sha256']==record_hash for c in cache):continue
        checkpoint=Path(r['checkpoint']);checkpoint_hash=digest(checkpoint)
        with np.load(checkpoint,allow_pickle=False) as z:meta=json.loads(str(z['meta']))
        if len(meta['history'])!=len(r['history']) or meta['couplings']!=r['couplings']:continue
        ent=spectra(Z.MPS.load(checkpoint))
        if digest(path)!=record_hash or digest(checkpoint)!=checkpoint_hash:
            print('CHANGED; discarded',path.name,flush=True);continue
        error=float(np.max(np.abs(np.array(ent['1'])-r['Sprof'])))
        assert error<1e-9,(path,error)
        row=dict(name=r['name'],L=r['L'],parity=r['parity'],seed=r['seed'],chi=r['chi'],
                 E=r['E_physical'],source_record=str(path),record_sha256=record_hash,
                 sweep_converged=r['convergence']['sweep_converged'],status_at_measurement=r['status'],
                 checkpoint=str(checkpoint),checkpoint_sha256=checkpoint_hash,
                 von_neumann_reproduction_error=error,profiles=ent,
                 summaries={a:stats(v,r['L'],float(a)) for a,v in ent.items()})
        cache.append(row)
        D.atomic_json(output,cache)
        print(r['name'],r['L'],r['parity'],r['chi'],{a:s['S_fractional'] for a,s in row['summaries'].items()},flush=True)
    # Rebuild summary convention even for cached rows from an earlier script version.
    archive_sources={digest(p):p for p in HERE.glob('**/validation_snapshots/*.json')}
    for r in cache:
        r['summaries']={a:stats(v,r['L'],float(a)) for a,v in r['profiles'].items()}
        if r.get('sweep_converged') is None:
            path=Path(r['source_record'])
            if (not path.exists() or digest(path)!=r['record_sha256']) and r['record_sha256'] in archive_sources:
                path=archive_sources[r['record_sha256']]
                r['source_record_snapshot']=str(path)
            if path.exists() and digest(path)==r['record_sha256']:
                raw=json.loads(path.read_text())
                r['sweep_converged']=raw['convergence']['sweep_converged']
                r.pop('convergence_note',None)
            else:
                r['sweep_converged']=None
                r['convergence_note']='Historical cache retained, but excluded from selection without a verified convergence flag.'
    D.atomic_json(output,cache)
    slopes=[]
    for name in sorted({r['name'] for r in cache}):
        rs=[r for r in cache if r['name']==name and r.get('sweep_converged') is True]
        # Highest chi per (L,P,seed), then the lower variational energy.
        bykey={}
        for r in rs:
            key=(r['L'],r['parity'],r['seed'])
            if key not in bykey or (r['chi'],-r['E'])>(bykey[key]['chi'],-bykey[key]['E']):bykey[key]=r
        rs=list(bykey.values());sizes=sorted({r['L'] for r in rs})
        gs=[]
        for L in sizes:
            g=min((r for r in rs if r['L']==L),key=lambda r:r['E']).copy()
            g['completed_converged_parities']=sorted({r['parity'] for r in rs if r['L']==L})
            gs.append(g)
        for i,a in enumerate(gs):
            for b in gs[i+1:]:
                if b!=gs[i+1] and b['L']!=2*a['L']:continue
                slopes.append(dict(name=name,sizes=[a['L'],b['L']],parities=[a['parity'],b['parity']],
                    chis=[a['chi'],b['chi']],
                    selected_sweeps_converged=True,
                    both_parities_converged=all(r['completed_converged_parities']==[0,1] for r in (a,b)),
                    c_eff={str(alpha):{key:12*(b['summaries'][str(alpha)][key]-a['summaries'][str(alpha)][key])/((1+1/alpha)*math.log(b['L']/a['L']))
                        for key in ('S_central3','S_fractional')} for alpha in (1,2,3)}))
    D.atomic_json(HERE/'renyi_slopes.json',slopes)
    exact=[]
    for legs in (1,2):
        for L in (24,48,96,192):
            n=L//legs
            h=np.diag(-np.ones(n-1),1)+np.diag(-np.ones(n-1),-1)
            e,v=np.linalg.eigh(h);C=v[:,e<0]@v[:,e<0].T
            profiles={str(a):[0.] for a in (1,2,3)}
            for l in range(1,n):
                p=np.clip(np.linalg.eigvalsh(C[:l,:l]),1e-16,1-1e-16)
                profiles['1'].append(float(-np.sum(p*np.log(p)+(1-p)*np.log(1-p))))
                for a in (2,3):profiles[str(a)].append(float(np.sum(np.log(p**a+(1-p)**a))/(1-a)))
            for a in profiles:profiles[a].append(0.)
            if legs==1:ent={a:p[1:-1] for a,p in profiles.items()}
            else:ent={a:[p[(l+1)//2]+p[l//2] for l in range(1,L)] for a,p in profiles.items()}
            exact.append(dict(legs=legs,L=L,profiles=ent,summaries={a:stats(v,L,float(a)) for a,v in ent.items()}))
    pairs=[]
    for legs in (1,2):
        rs=[r for r in exact if r['legs']==legs]
        for a,b in zip(rs,rs[1:]):
            pairs.append(dict(legs=legs,sizes=[a['L'],b['L']],c_eff={str(alpha):{
                key:12*(b['summaries'][str(alpha)][key]-a['summaries'][str(alpha)][key])/((1+1/alpha)*math.log(b['L']/a['L']))
                for key in ('S_central3','S_fractional')} for alpha in (1,2,3)}))
    D.atomic_json(HERE/'renyi_exact_controls.json',dict(rows=exact,pairs=pairs))

if __name__=='__main__':main()
