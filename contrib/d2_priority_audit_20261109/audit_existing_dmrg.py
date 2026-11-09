"""Read-only audit of committed C181/C187 finite-DMRG records.

No simulations are launched. Implement entropy estimators independently,
retain the actual model config, and flag incomplete sweep schedules.
"""
from pathlib import Path
import json, re, hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
DIRECTORIES=['codex/work_C181/typhon_out_bf/out_bf','codex/work_C187q/desktop_run3/out']

def sbar(S,L,W):
    cut=np.arange(1,L)
    return float(np.asarray(S)[abs(cut-L/2)<=W].mean())

def cprof(S,L,smooth=12,edge=1/6):
    cut=np.arange(1,L,dtype=float); S=np.asarray(S)
    if smooth>1:
        kernel=np.ones(smooth)/smooth
        cut=np.convolve(cut,kernel,'valid'); S=np.convolve(S,kernel,'valid')
    x=np.log(2*L/np.pi*np.sin(np.pi*cut/L))
    ok=(cut>=edge*L)&(cut<=(1-edge)*L)
    return float(6*np.polyfit(x[ok],S[ok],1)[0])

rows=[]; full={}
for directory in DIRECTORIES:
    for path in sorted((ROOT/directory).glob('*.jsonl')):
        raw=path.read_bytes(); records=[json.loads(line) for line in raw.splitlines() if line.strip()]
        r=records[-1]; cfg=r['cfg']; L,N=cfg['L'],cfg['N']; tag=r['tag'].removeprefix('bf_')
        bychi={}
        for rec in records: bychi[rec['chi']]=rec
        phys={k:cfg.get(k,default) for k,default in dict(wa=1.,phi=np.pi/2,t2=-1.,V2=-2.,V3=-2/3,blockade=True).items()}
        row=dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(raw).hexdigest(),tag=tag,L=L,N=N,
                 physical_config=phys,sweeps=len(records),scheduled_sweeps=len(cfg['chis']),
                 complete=len(records)==len(cfg['chis']),chi=r['chi'],energy=r['E'],truncation=r['trunc'],
                 sbar={str(w):sbar(r['S'],L,w) for w in (1,6,12)},
                 cprof={f'smooth{sm}_edge{edge:.4f}':cprof(r['S'],L,sm,edge)
                        for sm in (1,8,12,16) for edge in (1/6,1/4,1/3)},
                 by_chi=[dict(chi=k,E=v['E'],trunc=v['trunc'],sbar6=sbar(v['S'],L,6))
                         for k,v in bychi.items()],
                 final_sweep_dE=r['E']-records[-2]['E'],
                 final_sweep_dSbar6=sbar(r['S'],L,6)-sbar(records[-2]['S'],L,6))
        if 'meas' in r:
            meas=r['meas']; dens=np.asarray(meas['dens']); ds=np.convolve(dens,np.ones(16)/16,'valid')
            bulk=ds[len(ds)//4:3*len(ds)//4]
            row.update(polarization=abs(meas['Wa_f']-meas['Wb_f'])/(meas['Wa_f']+meas['Wb_f']),
                       central_density=float(dens[L//4:3*L//4].mean()),
                       smoothed_bulk_relative_range=float(np.ptp(bulk)/bulk.mean()))
        rows.append(row); full[tag]=(row,r)

twosize=[]; gaps=[]
for label in ('ctl1','ctl1fix','ctl2','w1n316_p1.5300','w1n316_p1.5550'):
    main={L:full[f'{label}_L{L}_N{3*L//16}'] for L in (64,128,192)
          if f'{label}_L{L}_N{3*L//16}' in full}
    for L,(row,r) in main.items():
        N=r['cfg']['N']; neighbors=[full.get(f'{label}_L{L}_N{N+d}') for d in (-1,1)]
        if all(neighbors):
            em,ep=(v[0]['energy'] for v in neighbors)
            gap=ep+em-2*r['E']; chemical=(ep-em)/2
            gaps.append(dict(label=label,L=L,all_complete=row['complete'] and all(v[0]['complete'] for v in neighbors),
                             gap=gap,Lgap=L*gap,chemical_potential=chemical,
                             hz_V40_leading=-20+(3+chemical/2)/80))
    for L1,L2 in ((64,128),(128,192)):
        if L1 not in main or L2 not in main: continue
        a,b=main[L1][0],main[L2][0]
        twosize.append(dict(label=label,L1=L1,L2=L2,both_complete=a['complete'] and b['complete'],
                            matching_hamiltonians=a['physical_config']==b['physical_config'],
                            c_two={str(w):6*(b['sbar'][str(w)]-a['sbar'][str(w)])/np.log(L2/L1) for w in (1,6,12)}))

out=dict(scope='Committed output snapshot; statistical fit slopes are not certified central charges.',
         rows=rows,two_size=twosize,parity_gap=gaps)
(HERE/'existing_dmrg_audit.json').write_text(json.dumps(out,indent=2)+'\n')
for x in twosize: print(json.dumps(x))
for x in gaps: print(json.dumps(x))
for row in rows:
    if row['L']==128 and row['N']==24:
        print(json.dumps({k:row.get(k) for k in ('tag','truncation','polarization','central_density',
                  'smoothed_bulk_relative_range','final_sweep_dSbar6','cprof')}))
