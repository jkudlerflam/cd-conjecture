#!/usr/bin/env python3
"""Audit archived physical-chain C184/C152 tables with current entropy windows."""
import json
import math
from collections import defaultdict
import analyze_saved as A

def main():
    groups=defaultdict(dict)
    folders=['codex/work_C184/typhon_out/out','codex/work_C184/local_out','codex/work_C184/c152_j2_copy']
    for folder in folders:
        for path in sorted((A.ROOT/folder).glob('dmrg_*.jsonl')):
            for lineno,line in enumerate(path.read_text().splitlines(),1):
                r=json.loads(line)
                if not all(k in r for k in ['p','S','L','chis']):continue
                key=(folder,tuple(round(x,10) for x in r['p']))
                old=groups[key].get(r['L'])
                if old and old['E']<=r['E']:continue
                vals=[abs(x) for x in r['p'] if abs(x)>1e-9]
                groups[key][r['L']]=dict(label=r['label'],L=r['L'],E=r['E'],p=r['p'],R=max(vals)/min(vals),
                    chi=r['chis'][-1],trunc=r['trunc'],dE_last=r.get('dE_last'),
                    stats=A.D.entropy_stats(r['S'],r['L']),source=f'{path.relative_to(A.ROOT)}:{lineno}',
                    S_four_cut=sum(r['S'][r['L']//2-3:r['L']//2+1])/4)
    output=[]
    for (folder,p),dd in groups.items():
        rows=[dd[L] for L in sorted(dd)]
        if len(rows)<2 or rows[0]['R']>10:continue
        slopes=[]
        for a,b in zip(rows,rows[1:]):
            ss={key:6*(b['stats'][key]-a['stats'][key])/math.log(b['L']/a['L'])
                for key in ['S_central3','S_fractional']}
            ss['four_cut']=6*(b['S_four_cut']-a['S_four_cut'])/math.log(b['L']/a['L'])
            slopes.append(dict(sizes=[a['L'],b['L']],**ss))
        output.append(dict(label=rows[0]['label'],folder=folder,R=rows[0]['R'],p=list(p),rows=rows,slopes=slopes))
    output.sort(key=lambda r:(r['R']>3,-max(s['S_fractional'] for s in r['slopes'])))
    (A.HERE/'archived_lead_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    for r in output:
        if max(s['S_fractional'] for s in r['slopes'])<1.1:continue
        print(r['label'],r['folder'],'R',round(r['R'],5),
              [(s['sizes'],round(s['S_central3'],4),round(s['S_fractional'],4)) for s in r['slopes']],flush=True)

if __name__=='__main__':main()
