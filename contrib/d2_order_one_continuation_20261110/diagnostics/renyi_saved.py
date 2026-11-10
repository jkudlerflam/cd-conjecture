#!/usr/bin/env python3
"""Rényi entropy cross-check on fixed states, with explicit input hashes."""
import hashlib
import json
import math
import analyze_saved as A
Z=A.Z

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    dest=A.HERE/'renyi_diagnostics.json'
    out=json.loads(dest.read_text()) if dest.exists() else dict(rows=[])
    names=['bridge_r4_u1_k4','bridge_r1.5_u1.5_k3','bridge_r2_u1.5_k4','control_xx','control_two_xx']
    inputs=[]
    for fn in ['saved_diagnostics.json','attractive_diagnostics.json','bridge_diagnostics.json']:
        data=json.loads((A.HERE/fn).read_text())
        inputs += [r for r in (data if isinstance(data,list) else data['rows']) if fn!='saved_diagnostics.json' or r['name'] in names]
    for r in inputs:
        if any(x['record_sha256']==r['record_sha256'] for x in out['rows']):continue
        path=A.ROOT/r['source_record'];ck=A.Path(r['checkpoint'])
        if sha(path)!=r['record_sha256'] or sha(ck)!=r['checkpoint_sha256']:continue
        cuts=A.schmidt_diagnostics(Z.MPS.load(ck));stats={}
        for field,fac in [('S',1.),('S2',8/6),('S3',9/6)]:
            stats[field]=A.D.entropy_stats([x[field] for x in cuts],r['L'])
            for p in stats[field]['profiles'].values():p['c']*=fac
        row={k:r[k] for k in ['name','L','parity','E','record_sha256','checkpoint_sha256','source_record']}
        row.update(stats=stats,cuts=[{k:c[k] for k in ['cut','S','S2','S3']} for c in cuts])
        if sha(path)!=r['record_sha256'] or sha(ck)!=r['checkpoint_sha256']:continue
        out['rows']=[x for x in out['rows'] if x['source_record']!=row['source_record']]+[row]
        out['slopes']=[]
        for name in sorted(set(x['name'] for x in out['rows'])):
            rows=[x for x in out['rows'] if x['name']==name]
            ground=[min([x for x in rows if x['L']==L],key=lambda x:x['E']) for L in sorted(set(x['L'] for x in rows))]
            for a,b in zip(ground,ground[1:]):
                s=dict(name=name,sizes=[a['L'],b['L']],parities=[a['parity'],b['parity']],c={})
                for field,fac in [('S',6),('S2',8),('S3',9)]:
                    s['c'][field]={k:fac*(b['stats'][field][k]-a['stats'][field][k])/math.log(b['L']/a['L'])
                        for k in ['S_central3','S_fractional']}
                out['slopes'].append(s)
        dest.write_text(json.dumps(out,indent=2)+'\n')
        print(r['name'],r['L'],r['parity'],{k:v['S_fractional'] for k,v in stats.items()},flush=True)

if __name__=='__main__':main()
