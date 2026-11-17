#!/usr/bin/env python3
"""Reuse saved evidence with literal coefficients and source hashes."""
import hashlib,json
from common import HERE,ROOT,coefficients


def source(p):return dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def main():
    files=[ROOT/p for p in ['STATUS.md','REPORT.md','codex/work_C190/models.py','codex/work_C190/RESULT.md',
       'codex/work_C195/RESULT.md','codex/work_C177/RESULT.md','codex/work_C177/tune.json',
       'codex/work_C179/REFEREE.md','contrib/qubit_convergence_20261116/audit.json']]
    old=[]
    for V in (8,12,16,24):
        for L in (64,128,192,256):
            p=ROOT/f'codex/work_C177/typhon_out/oc/F{V}_{L}.json'
            if not p.exists():continue
            d=json.loads(p.read_text());b=d['couplings']
            old.append(dict(source=source(p),V=V,L=L,kappa=4*V*b['D'],a=V*(b['hz']+V/2),
                couplings=b,Sbar=d['Sbar'],n_mid=d['n_mid'],c_prof=d['c_prof'],trunc=d['trunc'],hist=d['hist']))
    fes=[]
    for V in (8,12,16,24):
        for n in (2,4):
            p=ROOT/f'codex/work_C177/typhon_out/vc/F{V}_n{n}.json'
            if p.exists():
                d=json.loads(p.read_text());fes.append(dict(source=source(p),V=V,cell=n,results=d['results']))
    oldaudit=json.loads((ROOT/'contrib/qubit_convergence_20261116/audit.json').read_text())
    out=dict(source_files=[source(p) for p in files],production=coefficients(40,hz=-19.976047250613814),
             prior_audit_comparisons=oldaudit['comparisons'],prior_new_chi_test=oldaudit['new_chi_test'],
             old_open=old,old_infinite=fes,
             decision='Start V=8 and V=12, kappa=1, bracket a with small rationals; preserve old convergence caveats.')
    (HERE/'source_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('production',out['production']);print('old open rows',len(old),'infinite rows',len(fes))
    for V in (8,12,16,24):
        r=next(x for x in old if x['V']==V);print(V,'a',r['a'],'kappa',r['kappa'])


if __name__=='__main__':main()
