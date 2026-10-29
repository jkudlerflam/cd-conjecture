"""Summarize saved C67 results without assigning an unverified central charge."""
from pathlib import Path
import json
import numpy as np

out={}
for filename in ['grid_ed.jsonl','junction_ed.jsonl','controls.jsonl','fixed_ed.jsonl']:
    if not Path(filename).exists(): continue
    rows=[json.loads(s) for s in Path(filename).read_text().splitlines()]
    group={}
    for name in sorted({r['name'] for r in rows}):
        rr=sorted({r['L']:r for r in rows if r['name']==name}.values(),key=lambda r:r['L'])
        record={'rows':[],'size_fits':{},'size_fit_scope':
                ('Hamiltonian is retuned at each L; not a central-charge measurement'
                 if filename=='junction_ed.jsonl' else
                 'Finite-size slope only; inspect sector changes before interpreting')}
        for r in rr:
            L=r['L']; e={(s['p'],s['m']):s['e'] for s in r['levels']}
            odd=e[1,0][0]
            energies=sorted(x for es in e.values() for x in es)
            record['rows'].append(dict(L=L,v=r['v'],c_chord=r['c_chord'],
                half_entropy=r['S'][str(L//2)],p=r['p'],m=r['m'],
                scaled_low_gaps=[L*(x-r['E0'])/(2*np.pi) for x in energies[1:7]],
                odd_descendant_differences=[L*(e[1,m][0]-odd)/(2*np.pi) for m in [1,L-1]],
                E0=r['E0'],residual=r['residual']))
        for label,ss in [('all',rr),('Lmod4_0',[r for r in rr if r['L']%4==0]),
                         ('Lmod4_2',[r for r in rr if r['L']%4==2]),
                         ('Lge12',[r for r in rr if r['L']>=12])]:
            if len(ss)>=2:
                record['size_fits'][label]=float(3*np.polyfit(np.log([r['L'] for r in ss]),
                    [r['S'][str(r['L']//2)] for r in ss],1)[0])
        group[name]=record
    out[filename]=group

grid=Path('../../contrib/d2_grid/data/fes')
name='Z2_1_358_n11'
if (grid/(name+'.json')).exists():
    d=json.loads((grid/(name+'.json')).read_text())
    rs=[r for r in d['results'] if r['err']<1e-3 and np.isfinite(r['xi']) and r['xi']>0]
    S=np.array([np.mean(r['S']) for r in rs]); xi=np.array([r['xi'] for r in rs])
    out['inherited_grid_fes']={'file':name+'.json','chi':[r['chi'] for r in rs],
        'S':S.tolist(),'xi':xi.tolist(),'error':[r['err'] for r in rs],
        'c_fit':float(6*np.polyfit(np.log(xi),S,1)[0]),
        'c_pair':(6*np.diff(S)/np.diff(np.log(xi))).tolist()}
Path('summary.json').write_text(json.dumps(out,indent=2)+'\n')
for fn,groups in out.items():
    if not fn.endswith('jsonl'): print(fn,groups); continue
    for name,g in groups.items():
        print(fn,name,'last chord',g['rows'][-1]['c_chord'],'size',g['size_fits'])
