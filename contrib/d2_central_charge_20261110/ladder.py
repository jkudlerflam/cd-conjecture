"""Sequential, checkpointed ladder. One numerical worker per invocation."""
import argparse,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--geometry',choices=['infinite','periodic'],required=True)
ap.add_argument('--model',default='candidate');ap.add_argument('--L',type=int,required=True)
ap.add_argument('--parity',type=int,default=0);ap.add_argument('--seed',type=int,default=17)
ap.add_argument('--chis',type=int,nargs='+',required=True);ap.add_argument('--sweeps',type=int,default=100)
ap.add_argument('--load');ap.add_argument('--gauge',action='store_true');a=ap.parse_args();load=a.load
for chi in a.chis:
    label=f'{a.model}_{a.geometry}_L{a.L}_P{a.parity}_s{a.seed}_chi{chi}'
    if a.gauge:label+='_g'
    out=HERE/'results'/f'{label}.json'
    if out.exists():
        r=json.loads(out.read_text())
        if r['status']!='done':raise RuntimeError(f'Record is not finished: {out}')
    else:
        cmd=[sys.executable,str(HERE/'run.py'),'--geometry',a.geometry,'--model',a.model,
            '--L',str(a.L),'--parity',str(a.parity),'--seed',str(a.seed),'--chi',str(chi),'--sweeps',str(a.sweeps)]
        if load:cmd+=['--load',str(load)]
        if a.gauge:cmd+=['--gauge']
        print('START',label,flush=True);start=time.time()
        with (HERE/'logs'/f'{label}.log').open('w') as f:
            rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
        if rc:raise RuntimeError(f'{label} exit {rc}')
        r=json.loads(out.read_text())
        print('DONE',label,'converged',r['converged'],'seconds',time.time()-start,flush=True)
    load=r['checkpoint']
    # Never continue from a numerically unstable environment or nonfinite energy.
    if not r['convergence_checks']['energy_agreement'] or r['norm_before']>1e-3:
        raise RuntimeError(f'Unstable environment/state; inspect {out}')
