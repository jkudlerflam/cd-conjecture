#!/usr/bin/env python3
"""Saved-record consistency, provenance and completed resource inventory."""
import ast,hashlib,json,platform
from pathlib import Path
import numpy as np
import scipy
from threadpoolctl import threadpool_limits
from phase_common import HERE,ROOT,write_json


def main():
    checks=[];runs=[];running=[]
    with threadpool_limits(limits=1):
        for p in sorted((HERE/'observations').glob('*.json')):
            d=json.loads(p.read_text());ops=d['operators'];dens=np.array(ops['density']['profile'])
            assert min(dens)>-1e-8 and max(dens)<1+1e-8
            fields=np.full(d['L'],-1.75)
            if d['boundary']=='half_end_fields':fields[[0,-1]]/=2
            e=sum(ops['bond']['profile'])+fields@(1-2*dens);ee=abs(e-d['E']);assert ee<1e-8
            errors={}
            for name,o in ops.items():
                m=np.array(o['profile'])[o['indices']];raw=np.array(o['raw']);C=np.array(o['connected'])
                ce=float(np.max(abs(raw-np.outer(m,m)-C)));sym=float(np.max(abs(C-C.T)))
                ev=float(np.min(np.linalg.eigvalsh(C)));assert ce<1e-10 and sym<1e-10 and ev>-1e-8
                errors[name]=dict(covariance_identity_error=ce,symmetry_error=sym,minimum_covariance_eigenvalue=ev)
            n=ops['density'];ids=n['indices'];raw=np.array(n['raw']);m=dens[ids]
            assert np.max(abs(raw.diagonal()-m))<1e-8
            assert min(raw.flat)>-1e-8 and np.max(raw-np.minimum(m[:,None],m[None,:]))<1e-8
            pc=0.
            for anchor in d['pair']:
                i=anchor['i'];i0=ids.index(i);qi=raw[i0,ids.index(i+2)]
                for r in anchor['correlations']:
                    j=r['j'];qj=raw[ids.index(j),ids.index(j+2)]
                    z=complex(*r['raw']);dis=complex(*r['disconnected']);c=complex(*r['connected'])
                    assert abs(z-dis-c)<1e-10 and abs(z)**2<=qi*qj+1e-8
                    pc=max(pc,abs(z-dis-c))
            if 'parity' in d:assert abs(d['parity']-(1-2*d['P']))<1e-8
            checks.append(dict(file=str(p.relative_to(HERE)),energy_contraction_error=ee,operators=errors,pair_subtraction_error=pc,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        for p in sorted((HERE/'outputs').glob('*.json')):
            d=json.loads(p.read_text());assert d['cfg']==json.loads((HERE/'configs'/(p.stem+'.json')).read_text())
            if d['status']=='done':runs.append(dict(tag=p.stem,sweeps=len(d['hist']),sweep_seconds=sum(r['t'] for r in d['hist']),
                runtime_seconds=d['runtime_seconds'],peak_RSS_GB=d['peak_RSS_GB'],sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
            else:running.append(p.stem)
    scripts=list(HERE.glob('*.py'))+[ROOT/'codex/work_C190/z2dmrg.py',ROOT/'codex/work_C195/z2dmrg_xp.py']
    for p in HERE.glob('*.py'):ast.parse(p.read_text())
    files=[p for p in HERE.rglob('*') if p.is_file() and not any(x in p.parts for x in ('runtime','__pycache__'))]
    assert max(p.stat().st_size for p in files)<50_000_000
    out=dict(checked_observations=len(checks),checks=checks,completed_runs=runs,running=running,
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        completed_sweep_seconds=sum(r['sweep_seconds'] for r in runs),
        maximum_completed_RSS_GB=max((r['peak_RSS_GB'] for r in runs),default=0),
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in scripts},
        scope='Internal consistency and local numerical checks, not a ground-state or CFT certification.')
    write_json(HERE/'audit_results.json',out)
    print('Passed',len(checks),'observation checks;',len(runs),'completed production runs; running',running)


if __name__=='__main__':main()
