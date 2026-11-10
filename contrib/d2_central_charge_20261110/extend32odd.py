"""Retain the known serialization failure, then refine its final state."""
import json,pickle,subprocess,sys,time,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
source=HERE/'results/candidate_periodic_L32_P1_s17_chi128.json'
log=HERE/'logs/pbc32_odd.log';start=time.monotonic()
while time.monotonic()-start<7200:
    r=json.loads(source.read_text())
    if r['status']=='done':break
    if 'TypeError: Object of type Config is not JSON serializable' in log.read_text():
        ck=Path(r['checkpoint'])
        with ck.open('rb') as f:s=pickle.load(f)
        assert s['checkpoint_kind']=='final' and s['L']==32 and s['parity']==1
        r.update(status='postprocessing_failed',failure_stage='JSON serialization after final checkpoint',
            failure='TeNPy Config object in the original mutated options dictionary.',
            failure_log='logs/pbc32_odd.log',checkpoint_sha256=hashlib.sha256(ck.read_bytes()).hexdigest(),
            recovery='Fresh chi192 two-site stage and conditional one-site refinement; original remains excluded.')
        tmp=source.with_suffix('.tmp');tmp.write_text(json.dumps(r,indent=2)+'\n');tmp.replace(source)
        break
    if r['status']=='postprocessing_failed':break
    time.sleep(10)
else:raise RuntimeError('Old L32 odd worker did not finish within two hours')
cmd=[sys.executable,str(HERE/'run.py'),'--geometry','periodic','--L','32','--parity','1',
    '--chi','192','--sweeps','4','--min-sweeps','4','--load',r['checkpoint']]
print('START fresh chi192 record',flush=True)
with (HERE/'logs/pbc32_odd_refined.log').open('w') as f:
    rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
if rc:raise RuntimeError(f'chi192 exit {rc}')
cmd=[sys.executable,str(HERE/'refine_periodic.py'),'--L','32','--chi','192','--parity','1']
raise SystemExit(subprocess.call(cmd))
