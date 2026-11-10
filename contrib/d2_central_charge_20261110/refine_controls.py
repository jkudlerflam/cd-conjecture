"""Separate one-site Ising calibration after completed two-site states."""
import json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
for chi in [16,24,32,48]:
    label=f'ising_critical_infinite_L2_P0_s23_chi{chi}';source=HERE/'results'/f'{label}.json';start=time.monotonic()
    while time.monotonic()-start<7200:
        if source.exists():
            r=json.loads(source.read_text())
            if r['status']=='done':break
            if r['status']!='running':raise RuntimeError('Failed source: '+r['status'])
        time.sleep(10)
    else:raise RuntimeError('Source did not finish within two hours')
    target=HERE/'results'/f'{label}_one.json';assert not target.exists()
    cmd=[sys.executable,str(HERE/'run.py'),'--geometry','infinite','--model','ising_critical',
        '--L','2','--seed','23','--chi',str(chi),'--one-site','--min-sweeps','32',
        '--sweeps','600','--load',r['checkpoint']]
    print('START',label,flush=True)
    with (HERE/'logs'/f'{label}_one.log').open('w') as f:
        rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:raise RuntimeError(f'{label}: {rc}')
    print('DONE',label,flush=True)
