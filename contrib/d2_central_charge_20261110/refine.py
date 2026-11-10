"""Sequential fixed-chi one-site follow-up of completed two-site states."""
import argparse,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--chis',nargs='+',type=int,required=True)
ap.add_argument('--gauge',action='store_true');ap.add_argument('--sweeps',type=int,default=200)
a=ap.parse_args()
for chi in a.chis:
    label=f'candidate_infinite_L2_P0_s17_chi{chi}'+('_g' if a.gauge else '')
    source=HERE/'results'/f'{label}.json'
    if not source.exists() or json.loads(source.read_text())['status']!='done':
        print('DEFER unavailable completed source',label,flush=True);break
    target=HERE/'results'/f'{label}_one.json'
    if target.exists():
        print('EXISTS',target.name,flush=True);continue
    cmd=[sys.executable,str(HERE/'run.py'),'--geometry','infinite','--L','2','--chi',str(chi),
        '--one-site','--sweeps',str(a.sweeps),'--load',json.loads(source.read_text())['checkpoint']]
    if a.gauge:cmd+=['--gauge']
    print('START',label,flush=True)
    with (HERE/'logs'/f'{label}_one.log').open('w') as f:
        rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:raise RuntimeError(f'{label}: {rc}')
    print('DONE',label,flush=True)
