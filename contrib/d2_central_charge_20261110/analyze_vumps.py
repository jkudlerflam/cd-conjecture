"""Keep independently refined VUMPS states separate from the DMRG fits."""
import json,hashlib,math
from collections import defaultdict
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
inventory=[];groups=defaultdict(list)
for f in sorted((HERE/'vumps_results').glob('*.json')):
    r=json.loads(f.read_text());row={k:r.get(k) for k in ['status','converged','model','L','chi','gauge','E_physical','max_tangent','convergence_checks','error']}
    row.update(record=str(f.relative_to(HERE)),record_sha256=hashlib.sha256(f.read_bytes()).hexdigest())
    if r['status']=='done':
        row['S_renyi']={alpha:float(np.mean(r['entropy'][alpha])) for alpha in ['1','2','3']}
        row['S']=row['S_renyi']['1'];row['xi']=r['correlation_lengths'][0]
        source=ROOT/r['source'];base=json.loads(source.read_text());seen=set()
        while 'vumps_results' in source.parts:
            assert str(source) not in seen;seen.add(str(source));source=ROOT/base['source'];base=json.loads(source.read_text())
        row['original_dmrg_source']=str(source.relative_to(ROOT))
        row['delta_E_from_dmrg']=r['E_physical']-base['E_physical']
        row['delta_S_from_dmrg']=row['S']-float(np.mean(base['entropy']['1']))
        if r['converged']:groups[(r['model'],r['L'],r['gauge'])].append(row)
    inventory.append(row)
inventory.sort(key=lambda r:(r['model'],r['L'],r['chi'],Path(r['record']).stem.count('_vumps')))
families=[]
for key,rows in sorted(groups.items()):
    bychi=defaultdict(list)
    for r in rows:bychi[r['chi']].append(r)
    points=[min(rr,key=lambda r:r['E_physical']) for _,rr in sorted(bychi.items())]
    pairs=[]
    for a,b in zip(points,points[1:]):
        pairs.append(dict(chis=[a['chi'],b['chi']],c_renyi={alpha:12*float(alpha)/(float(alpha)+1)*(b['S_renyi'][alpha]-a['S_renyi'][alpha])/math.log(b['xi']/a['xi']) for alpha in ['1','2','3']}))
    families.append(dict(model=key[0],L=key[1],gauge=key[2],points=points,pairs=pairs))
out=dict(inventory=inventory,families=families,interpretation='Separate fixed-chi exact-environment refinements. Original failed records remain excluded; no CFT is presumed.')
(HERE/'vumps_summary.json').write_text(json.dumps(out,indent=2)+'\n')
lines=['# Independent VUMPS checks','','All rows are retained; only full passes enter the separate slope table. Changes are relative to the original DMRG source, including for warm retries.','',
       '| Model | χ | Stage | Pass | ΔE per site | ΔS₁ | S₁ | ξ | Tangent norm |','|---|---:|---|---|---:|---:|---:|---:|---:|']
for r in inventory:
    stage='Initial' if Path(r['record']).stem.count('_vumps')==1 else 'Warm retry'
    link=f"[{stage}]({r['record']})"
    if r['status']=='done':lines.append(f"| {r['model']} | {r['chi']} | {link} | {r['converged']} | {r['delta_E_from_dmrg']:.5g} | {r['delta_S_from_dmrg']:.5g} | {r['S']:.9f} | {r['xi']:.6f} | {r['max_tangent']:.5g} |")
    else:lines.append(f"| {r['model']} | {r['chi']} | {link} | {r['status']} | | | | | |")
if any(r['pairs'] for r in families):
    lines+=['','| Model | χ pair | Order 1 coefficient | Order 2 | Order 3 |','|---|---|---:|---:|---:|']
else:lines+=['','No family has two fully accepted bond dimensions, so no VUMPS slope is assigned.']
for r in families:
    for p in r['pairs']:
        c=p['c_renyi'];lines.append(f"| {r['model']} | {p['chis'][0]} → {p['chis'][1]} | {c['1']:.6f} | {c['2']:.6f} | {c['3']:.6f} |")
lines+=['','Source links, failed gates and checksums are in `vumps_summary.json`; raw iteration histories are in `vumps_results/`.','']
(HERE/'VUMPS_RESULTS.md').write_text('\n'.join(lines))
print([(r['model'],[(p['chi'],p['S'],p['xi']) for p in r['points']]) for r in families])
