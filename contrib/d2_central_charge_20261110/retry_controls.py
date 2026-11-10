"""One explicit warm retry for a completed failed Ising control point."""
import argparse,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--model',choices=['ising_critical','ising_gapped'],required=True);a=ap.parse_args()
for chi in [16,24,32,48]:
    label=f'{a.model}_infinite_L2_P0_s23_chi{chi}';source=HERE/'results'/f'{label}.json';start=time.monotonic()
    while time.monotonic()-start<7200:
        if source.exists():
            r=json.loads(source.read_text())
            if r['status']=='done':break
            if r['status']!='running':raise RuntimeError('Failed source: '+r['status'])
        time.sleep(10)
    else:raise RuntimeError('Source did not complete within two hours')
    if r['converged']:
        print('PASS original',label,flush=True);continue
    assert r['convergence_checks']['energy_agreement'] and r['norm_before']<1e-3
    target=HERE/'results'/f'{label}_a1.json'
    if target.exists():raise RuntimeError('Retry record already exists')
    cmd=[sys.executable,str(HERE/'run.py'),'--geometry','infinite','--model',a.model,
         '--L','2','--seed','23','--chi',str(chi),'--attempt','1','--sweeps','400','--load',r['checkpoint']]
    print('RETRY',label,flush=True)
    with (HERE/'logs'/f'{label}_a1.log').open('w') as f:
        rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:raise RuntimeError(f'{label}: {rc}')
    print('DONE retry',label,flush=True)
