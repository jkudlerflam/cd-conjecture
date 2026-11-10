"""Fixed-chi one-site audit after a periodic two-site worker completes."""
import argparse,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,required=True)
ap.add_argument('--chi',type=int,required=True);ap.add_argument('--parity',type=int,default=0)
ap.add_argument('--gauge',action='store_true');ap.add_argument('--sweeps',type=int,default=12)
a=ap.parse_args()
label=f'candidate_periodic_L{a.L}_P{a.parity}_s17_chi{a.chi}'+('_g' if a.gauge else '')
source=HERE/'results'/f'{label}.json';start=time.monotonic()
while time.monotonic()-start<7200:
    if source.exists():
        r=json.loads(source.read_text())
        if r['status']=='done':break
        if r['status']!='running':raise RuntimeError('Source failed: '+r['status'])
    time.sleep(10)
else:raise RuntimeError('Source did not complete within two hours')
if r['converged']:
    print('Full source gates already pass; no automatic refinement',flush=True)
    raise SystemExit(0)
cmd=[sys.executable,str(HERE/'run.py'),'--geometry','periodic','--L',str(a.L),
     '--chi',str(a.chi),'--parity',str(a.parity),'--one-site','--sweeps',str(a.sweeps),
     '--min-sweeps','4','--load',r['checkpoint']]
if a.gauge:cmd+=['--gauge']
print('START one-site follow-up',label,flush=True)
with (HERE/'logs'/f'{label}_one.log').open('w') as f:
    rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
if rc:raise RuntimeError(f'{label}: {rc}')
print('DONE',flush=True)
