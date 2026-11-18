"""Reanalyze only committed rounded tables; no unavailable MPS is inferred."""
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np

HERE=Path(__file__).resolve().parents[1]
SOURCE=HERE.parent/'d2_xy_yx_chirality_imps_20261118'
NUM=r'[+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|nan)'


def parse(path):
    runs=[]; current=None
    for line in path.read_text().splitlines():
        if line.startswith('== '):
            lab,branch=line[3:].split()
            current=dict(label=lab,branch=branch,rows=[])
            mh=re.search(r'_h(\d+\.\d+)_',lab)
            if mh: current['h']=float(mh[1])
            runs.append(current)
        elif line.lstrip().startswith('chi='):
            vals={k:float(v) for k,v in re.findall(r'(\w+)=\s*('+NUM+r')',line)}
            current['rows'].append({k:vals[k] for k in ('chi','E','S','xi','xi0','xi1','cat','m','B','Z')})
    return runs


def fit(rows,length):
    x=np.log([r[length] for r in rows]); y=np.array([r['S'] for r in rows])
    if len(x)<3 or np.any(~np.isfinite(x)) or np.any(np.diff(x)<=0): return None
    p=np.polyfit(x,y,1)
    out=dict(n=len(x),chi_min=int(rows[0]['chi']),chi_max=int(rows[-1]['chi']),
             c_slope=float(6*p[0]),S_max_residual=float(max(abs(y-np.polyval(p,x)))),
             xi_min=float(np.exp(x[0])),xi_max=float(np.exp(x[-1])),
             cat_changes=len(set(r['cat'] for r in rows))>1,
             xi_ratio_min=min(r['xi1']/r['xi0'] for r in rows),
             xi_ratio_max=max(r['xi1']/r['xi0'] for r in rows))
    for key in ('m','B'):
        val=np.abs([r[key] for r in rows])
        out[key+'_slope']=(float(np.polyfit(x,np.log(val),1)[0])
                           if np.all(np.isfinite(val)) and min(val)>2e-5 else None)
    return out


def main():
    allout={}; csvrows=[]
    for source in ('target','control'):
        path=SOURCE/'tables'/f'{source}_ladders.txt'; runs=parse(path)
        for run in runs:
            run['fits']={}
            for length in ('xi','xi0','xi1'):
                for chimin in (64,96,128,160):
                    sel=[r for r in run['rows'] if r['chi']>=chimin]
                    result=fit(sel,length)
                    run['fits'][f'{length}_chi{chimin}']=result
                    if result: csvrows.append(dict(source=source,label=run['label'],branch=run['branch'],
                                                    length=length,window=f'chi>={chimin}',**result))
                run['fits'][length+'_last4']=fit(run['rows'][-4:],length)
            run['adjacent']=[]
            for a,b in zip(run['rows'],run['rows'][1:]):
                dl=np.log(b['xi']/a['xi'])
                run['adjacent'].append(dict(chi0=a['chi'],chi1=b['chi'],
                                            c=6*(b['S']-a['S'])/dl if dl else None))
        allout[source]=dict(source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),runs=runs)
    groups={}
    for run in allout['target']['runs']:
        for row in run['rows']:
            groups.setdefault((run['h'],row['chi']),[]).append(dict(branch=run['branch'],**row))
    allout['branch_comparisons']=[]
    for (h,chi),rows in sorted(groups.items()):
        if len(rows)<2: continue
        best=min(rows,key=lambda r:r['E'])
        allout['branch_comparisons'].append(dict(h=h,chi=chi,best_branch=best['branch'],
                                                energy_spread=max(r['E'] for r in rows)-best['E'],
                                                S_spread=np.ptp([r['S'] for r in rows]),
                                                m_abs_spread=np.ptp([abs(r['m']) for r in rows]),
                                                B_abs_spread=np.ptp([abs(r['B']) for r in rows]),
                                                rows=rows))
    # h-dependent energy slopes are an independent normalization consistency
    # check, but the rounded table has no solver errors or statistical noise model.
    ladders={r['h']:r for r in allout['target']['runs'] if r['branch']=='_pin'}
    hs=sorted(ladders); hf=[]
    for h0,h1 in zip(hs,hs[1:]):
        r0={r['chi']:r for r in ladders[h0]['rows']}
        r1={r['chi']:r for r in ladders[h1]['rows']}
        for chi in sorted(r0.keys()&r1.keys()):
            a,b=r0[chi],r1[chi]
            hf.append(dict(h0=h0,h1=h1,chi=chi,minus_de_dh=-(b['E']-a['E'])/(h1-h0),
                           mean_Z=(a['Z']+b['Z'])/2))
    allout['hellmann_feynman_secants']=hf
    def clean(x):
        if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
        if isinstance(x,list): return [clean(v) for v in x]
        if isinstance(x,(float,np.floating)) and not np.isfinite(x): return None
        return x
    (HERE/'results'/'ladder_audit.json').write_text(json.dumps(clean(allout),indent=2,allow_nan=False)+'\n')
    with (HERE/'results'/'ladder_fits.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(csvrows[0]),lineterminator='\n');writer.writeheader();writer.writerows(csvrows)
    for run in allout['target']['runs']:
        f=run['fits']['xi_chi96']; f2=run['fits']['xi_chi128']; f0=run['fits']['xi0_chi96']
        if f:
            print(f'{run["h"]:.5f} {run["branch"]:7} n={f["n"]} c96={f["c_slope"]:.4f} '
                  f'c128={f2["c_slope"] if f2 else float("nan"):.4f} c0={f0["c_slope"]:.4f} '
                  f'mexp={f["m_slope"]} Bexp={f["B_slope"]} ratio={f["xi_ratio_max"]:.3f}')


if __name__=='__main__':main()
