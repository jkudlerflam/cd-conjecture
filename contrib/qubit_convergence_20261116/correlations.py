#!/usr/bin/env python3
"""Recompute the frozen C193 baseline; do not infer missing L=512 observables."""
import importlib.util
import json
from pathlib import Path
from threadpoolctl import threadpool_limits
from audit import ROOT,HERE,sha


def main():
    base=ROOT/'codex/work_C193/ckpt_meas'
    spec=importlib.util.spec_from_file_location('ana193',base/'ana193.py')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    fz=json.loads((base/'frozen.json').read_text())
    rows=[]
    files=sorted((base/'typhon_out/out').glob('sym40_n316_L*_P*_meas.json'))
    files += sorted(p for p in (base/'typhon_out/out').glob('xx*_meas.json') if '_L64_' not in p.name)
    with threadpool_limits(limits=1):
        for f in files:
            dat=mod.from_qubit(str(f)); est=mod.frozen_estimates(dat,fz)
            rows.append(dict(file=str(f.relative_to(ROOT)),sha256=sha(f),L=dat['L'],P=dat['P'],E=dat['E'],
                             Dmax=dat['Dmax'],tag=dat['tag'],checkpoint_checks=dat['check'],estimates=est))
    sym=[r for r in rows if r['tag'].startswith('sym40')]
    ground={L:min([r for r in sym if r['L']==L],key=lambda r:r['E']) for L in (64,128,256)}
    drift={}
    for mode in ('ground','P0','P1'):
        rs=[ground[L] if mode=='ground' else min((r for r in sym if r['L']==L and r['P']==int(mode[-1])),key=lambda r:r['E']) for L in (64,128,256)]
        drift[mode]={}
        for label in ('pair','sf'):
            K=[r['estimates'][f'Ks_{label}'] for r in rs]; dK=[r['estimates'][f'dKs_{label}'] for r in rs]
            b1=mod.bound(K[0],K[1],dK[0],dK[1]); b2=mod.bound(K[1],K[2],dK[1],dK[2])
            drift[mode][label]=dict(K=K,errors=dK,change_64_128=b1,change_128_256=b2,
                B2_clauses=[K[2]>=1.1,b2['ypl']>b2['y1u'],b2['dK']>=b1['dK']],
                B3_clauses=[K[2]<1,b2['y1u']>=b2['ypl'],b2['dK']<b1['dK']-(b1['sig']**2+b2['sig']**2)**.5])
    out=dict(frozen_sha256=sha(base/'frozen.json'),analyzer_sha256=sha(base/'ana193.py'),rows=rows,drift=drift,
        new_256_512_Ks_test='NOT EVALUATED: no large-chi checkpoints or full C193 correlators at L=512',
        caution='Old C193 chi<=768 data; not new measurements. K_s^sf at L=64 fails the control calibration.')
    (HERE/'correlations.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(drift,indent=2))


if __name__=='__main__':main()
