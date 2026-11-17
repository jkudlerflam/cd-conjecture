#!/usr/bin/env python3
"""Exact open XX entropies and saved same-solver control audit.

For two disjoint legs the bipartition entropy is the sum of leg entropies.
No false identification of the spin next-neighbour hop with a single JW chain.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits
from audit import HERE, ROOT, sbar, sha, ctwo


def exact(L,mu,legs):
    n=L//legs
    k=np.arange(1,n+1)
    levels=-2*np.cos(np.pi*k/(n+1))-mu
    if min(abs(levels))<1e-10: raise ValueError('degenerate Fermi level')
    N=int(np.count_nonzero(levels<0))
    U=np.sqrt(2/(n+1))*np.sin(np.pi*np.arange(1,n+1)[:,None]*k[:N]/(n+1))
    C=U@U.T
    ss=[0.]
    for l in range(1,n):
        ev=np.clip(eigvalsh(C[:l,:l]),1e-15,1-1e-15)
        ss.append(float(-np.sum(ev*np.log(ev)+(1-ev)*np.log1p(-ev))))
    ss.append(0.)
    S=ss[1:-1] if legs==1 else [ss[(l+1)//2]+ss[l//2] for l in range(1,L)]
    return dict(L=L,N=N*legs,E=float(legs*levels[:N].sum()),S=S,
                Sbar={str(w):sbar(S,L,w) for w in (1,6,12)})


def main():
    rows=[]
    paths=list((ROOT/'codex/work_C190/gpu/typhon_out/outG').glob('xx*_L*.json'))
    paths+=list((ROOT/'codex/work_C190/gpu/typhon_out/outX1024').glob('xx*_L*.json'))
    paths.append(ROOT/'codex/work_C190/typhon_out/outB/xx1_L256_P0_rand1.json')
    with threadpool_limits(limits=1):
        for p in sorted(paths):
            r=json.loads(p.read_text()); c=r['cfg']; L=c['L']; legs=1 if c['model']=='xx1' else 2
            ex=exact(L,c['mu'],legs)
            row=dict(file=str(p.relative_to(ROOT)),sha256=sha(p),L=L,model=c['model'],chi=r['hist'][-1]['chi'],
                     exact=ex,dmrg=dict(L=L,Sbar={str(w):sbar(r['S'],L,w) for w in (1,6,12)}),
                     error_Sbar={str(w):sbar(r['S'],L,w)-ex['Sbar'][str(w)] for w in (1,6,12)},
                     error_E=r.get('E_phys',r['hist'][-1]['E_phys'])-ex['E'],
                     status=r['status'])
            # Compact artifact: exact profiles are reproducible, no need to duplicate.
            row['exact'].pop('S'); rows.append(row)
    pairs=[]
    for model in ('xx1','xx2'):
        best={}
        for r in rows:
            if r['model']==model and (r['L'] not in best or r['chi']>best[r['L']]['chi']): best[r['L']]=r
        for L in sorted(best):
            if 2*L not in best: continue
            a,b=best[L],best[2*L]
            pairs.append(dict(model=model,Lsmall=L,Lbig=2*L,dmrg=ctwo(a['dmrg'],b['dmrg']),
                              exact=ctwo(a['exact'],b['exact']),files=[a['file'],b['file']]))
    planned=[]
    (HERE/'control_configs').mkdir(exist_ok=True)
    with threadpool_limits(limits=1):
        for legs in (1,2):
            per=[]
            for L in (256,512):
                n=L//legs; N=round(.1957*n)
                mu=float(-np.cos(np.pi*N/(n+1))-np.cos(np.pi*(N+1)/(n+1)))
                ex=exact(L,mu,legs); ex.pop('S'); per.append(ex)
                stage=lambda chi:dict(chi=chi,res_tol=1e-9,min_sweeps=3,max_sweeps=8,
                                      tolE=4e-9 if L==512 else 2e-9,tolS=5e-4)
                cfg=dict(tag=f'qc_xx{legs}_L{L}_P{N*legs%2}',model=f'xx{legs}',mu=mu,L=L,P=N*legs%2,
                    init='rand:1',sched=[[32,.001,1e-7,2],[64,.0001,1e-8,2],[128,.0001,1e-8,2],
                        [192,1e-5,1e-8,2],[256,1e-5,1e-8,2],[384,0.,1e-8,2],[512,0.,1e-9,2]],
                    conv=[stage(c) for c in (768,1024,1536,2048)],final=stage(2560),
                    maxit=300,mmax=32,variance=False,gwin=48,N_exact=N*legs)
                (HERE/'control_configs'/f'xx{legs}_L{L}.json').write_text(json.dumps(cfg,indent=2)+'\n')
            planned.append(dict(model=f'xx{legs}',points=per,exact_c_two=ctwo(per[0],per[1]),
                                 solver_status='NOT RUN; exact calibration and configs only'))
    out=dict(rows=rows,pairs=pairs,planned_matched_controls=planned,matched_256_512_C195_controls='NOT AVAILABLE',
             note='C198 two-leg data use fibdmrg and are not this matched-solver gate.')
    (HERE/'controls.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in rows: print(r['model'],r['L'],r['chi'],'S errors',r['error_Sbar'],'E error',r['error_E'])
    for p in pairs: print('c_two',p)


if __name__=='__main__': main()
