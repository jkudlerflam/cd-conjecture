#!/usr/bin/env python3
"""Bounded C190 DMRG with saved accuracy stages and bulk observables."""
import argparse,json,os,time,sys,resource
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from phase_common import HERE,ROOT,Z,build,write_json
from observables import observe


def run(cfg,root=HERE):
    tag=cfg['tag'];out=root/'outputs'/(tag+'.json');rt=root/'runtime';rt.mkdir(parents=True,exist_ok=True)
    ck=rt/(tag+'.npz');hist=[];start=time.time()
    if out.exists():
        previous=json.loads(out.read_text())
        if previous['status']=='done':return previous
        hist=previous['hist']
    if ck.exists():psi=Z.MPS.load(str(ck))
    elif cfg.get('parent_checkpoint'):
        path=Path(cfg['parent_checkpoint']);path=path if path.is_absolute() else root/path
        psi=Z.MPS.load(str(path))
    else:psi=Z.random_mps(cfg['L'],cfg['P'],D0=8,seed=cfg.get('seed',1),dtype=float)
    assert psi.L==cfg['L'] and psi.P==cfg['P'] and psi.check()
    mpo,gp,shift,bond=build(cfg);offset=mpo.const-cfg['L']*shift
    base=dict(cfg=cfg,bond=bond,shift=shift,E_offset=offset)
    def write(status,**extra):write_json(out,dict(base,hist=hist,status=status,**extra))
    def stable(rr,st):
        if len(rr)<st.get('min',3) or 'Sbar6' not in rr[-1] or 'Sbar6' not in rr[-2]:return False
        a,b=rr[-2:]
        return (abs(a['E']-b['E'])<2e-9 and max(abs(a[k]-b[k]) for k in ('Sbar1','Sbar6','Sbar12'))<5e-4
                and abs(a['n_mid']-b['n_mid'])<2e-4 and b['res_max']<1e-8)
    measured=[]
    for stage,st in enumerate(cfg['stages']):
        rr=[r for r in hist if r['stage']==stage]
        if not (len(rr)>=st['max'] or (st.get('alpha',0)==0 and stable(rr,st))):
            def log(rec):
                rec['stage']=stage;rec['E_phys']=rec['E']+offset;hist.append(rec)
                print(json.dumps(dict(tag=tag,stage=stage,**{k:v for k,v in rec.items() if k in
                    ('chi','E_phys','Sbar6','n_mid','res_max','trunc','t')})),flush=True)
            def checkpoint(p,_):
                tmp=str(ck)+'.tmp.npz';p.save(tmp,meta=dict(tag=tag,cfg=cfg,sweeps=len(hist)));os.replace(tmp,ck);write('running')
            specs=[dict(chi=st['chi'],alpha=st.get('alpha',0),res_tol=1e-9,maxit=300,mmax=32)
                   for _ in range(st['max']-len(rr))]
            psi,_=Z.dmrg(psi,mpo,specs,log=log,ckpt=checkpoint,
                         stop=lambda _:stable([r for r in hist if r['stage']==stage],st))
        if st.get('alpha',0)==0:
            stagecp=rt/(tag+f"_chi{st['chi']}.npz")
            if not stagecp.exists():psi.save(str(stagecp),meta=dict(cfg=cfg,stage=stage,chi=st['chi']))
            if st['chi'] in cfg.get('observe_chi',[]):
                of=root/'observations'/(tag+f"_chi{st['chi']}.json")
                if not of.exists():
                    state=Z.MPS.load(str(stagecp));obs=observe(state,gp,mpo.h2)
                    S=Z.entropies(state);E=Z.energy(state,mpo)+offset
                    obs.update(tag=tag,chi=st['chi'],boundary=cfg.get('boundary','full'),
                        cfg=cfg,E=E,S=S.tolist(),Sbar={str(w):Z.sbar(S,cfg['L'],w) for w in (1,6,12)},
                        parity=Z.parity_expectation(state))
                    assert abs(E-[r for r in hist if r['stage']==stage][-1]['E_phys'])<1e-7
                    write_json(of,obs);print('measured',of.name,obs['runtime_seconds'],flush=True)
                measured.append(str(of.relative_to(root)))
    overlap=None
    if cfg.get('compare_checkpoint'):
        refpath=Path(cfg['compare_checkpoint']);refpath=refpath if refpath.is_absolute() else ROOT/refpath
        if refpath.exists():
            ref=Z.MPS.load(str(refpath));assert ref.L==psi.L and ref.P==psi.P
            E=np.ones((1,1))
            for a,b in zip(psi.M,ref.M):E=Z._apply_site(E,a,None,Bc=b)
            ov=complex(E[0,0]);assert abs(ov)<=1+1e-8
            overlap=dict(reference=str(refpath),overlap=[ov.real,ov.imag],absolute_overlap=abs(ov),infidelity=max(0.,1-abs(ov)**2))
        else:overlap=dict(reference=str(refpath),status='checkpoint unavailable; compare saved observables instead')
    write('done',observations=measured,independent_overlap=overlap,runtime_seconds=time.time()-start,
        peak_RSS_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)/1e9)
    return json.loads(out.read_text())


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('config');ap.add_argument('--run-root',type=Path,default=HERE)
    args=ap.parse_args()
    with threadpool_limits(limits=1):run(json.loads(Path(args.config).read_text()),args.run_root)
