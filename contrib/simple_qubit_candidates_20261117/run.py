#!/usr/bin/env python3
"""Bounded parity-resolved C190 DMRG campaign; checkpoints stay ignored."""
import argparse,json,os,resource,time,sys
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE,Z,Ana,build,pair_diagnostics


def run(cfg,run_root=HERE):
    out=run_root/'outputs';rt=run_root/'runtime';out.mkdir(parents=True,exist_ok=True);rt.mkdir(parents=True,exist_ok=True)
    tag=cfg['tag']; fout=out/(tag+'.json'); fck=rt/(tag+'.npz')
    if fout.exists() and json.loads(fout.read_text()).get('status')=='done':return json.loads(fout.read_text())
    mpo,gph,shift,b=build(cfg);L=cfg['L'];P=cfg['P'];hist=[];start=time.time()
    if fout.exists() and fck.exists():
        hist=json.loads(fout.read_text())['hist'];psi=Z.MPS.load(str(fck))
    elif cfg.get('from'):
        parent=out/(cfg['from']+'.json')
        if not parent.exists() or json.loads(parent.read_text()).get('status')!='done':
            raise ValueError('Complete the parent run before starting its warm continuation: '+cfg['from'])
        psi=Z.MPS.load(str(rt/(cfg['from']+'.npz')));assert psi.L==L and psi.P==P
    else:psi=Z.random_mps(L,P,D0=8,seed=cfg.get('seed',1),dtype=float)
    assert psi.L==L and psi.P==P
    base=dict(cfg=cfg,bond=b,shift=shift,E_offset=mpo.const-L*shift)
    def write(status,**extra):
        tmp=fout.with_suffix('.tmp');tmp.write_text(json.dumps(dict(base,hist=hist,status=status,**extra),indent=2)+'\n');os.replace(tmp,fout)
    for stage,st in enumerate(cfg['stages']):
        recs=[r for r in hist if r['stage']==stage]
        def stationary(rr):
            if len(rr)<st.get('min',3):return False
            u,v=rr[-2:]
            return (abs(u['E']-v['E'])<2e-9 and
                    max(abs(u[k]-v[k]) for k in ('Sbar1','Sbar6','Sbar12'))<5e-4 and
                    abs(u['n_mid']-v['n_mid'])<2e-4 and v['res_max']<1e-8)
        if len(recs)>=st['max'] or (st.get('alpha',0)==0 and stationary(recs)):continue
        def log(rec):
            rec['stage']=stage;rec['E_phys']=rec['E']+base['E_offset'];hist.append(rec)
            print(json.dumps(dict(tag=tag,stage=stage,**{k:v for k,v in rec.items() if k in ['chi','E_phys','Sbar6','n_mid','trunc','res_max','t']})),flush=True)
        def ckpt(p,_):
            tmp=str(fck)+'.tmp.npz';p.save(tmp,meta=dict(sweeps=len(hist)));os.replace(tmp,fck);write('running')
        specs=[dict(chi=st['chi'],alpha=st.get('alpha',0.),res_tol=st.get('res_tol',1e-9),maxit=300,mmax=32) for _ in range(st['max']-len(recs))]
        psi,_=Z.dmrg(psi,mpo,specs,log=log,ckpt=ckpt,stop=lambda _:stationary([r for r in hist if r['stage']==stage]))
    S=Z.entropies(psi);meas=Z.measure(psi,mpo,gph=gph,gwin=min(24,L-4),cp_rmax=min(24,L//2),h2_list=mpo.h2,h1_list=mpo.h1_in)
    wa,wb=Ana.valley_weights(meas)
    final=dict(S=S.tolist(),meas=meas,Sbar={str(w):Z.sbar(S,L,w) for w in (1,6,12)},
        c_profile=Ana.c_prof(S,L),valley=dict(Wa=wa,Wb=wb,imbalance=(wa-wb)/(wa+wb)),
        pair=pair_diagnostics(psi,gph),parity=Z.parity_expectation(psi),
        runtime_seconds=time.time()-start,peak_RSS_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)/1e9)
    if cfg.get('variance'):
        var,e=Z.variance(psi,mpo)
        final['variance']=var;final['variance_energy']=e+base['E_offset']
        final['runtime_seconds']=time.time()-start
        final['peak_RSS_GB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)/1e9
    write('done',**final)
    return json.loads(fout.read_text())


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('config')
    ap.add_argument('--run-root',type=Path,default=HERE,help='Separate output/checkpoint directory for fresh reproduction.')
    args=ap.parse_args()
    with threadpool_limits(limits=1):run(json.loads(Path(args.config).read_text()),args.run_root)
