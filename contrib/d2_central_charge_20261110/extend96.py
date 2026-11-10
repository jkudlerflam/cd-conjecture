"""Registered chi96 extension after the real-gauge chi64 worker finishes."""
import json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
source=HERE/'results/candidate_infinite_L2_P0_s17_chi64_g.json'
start=time.monotonic()
while time.monotonic()-start<7200:
    if source.exists():
        r=json.loads(source.read_text())
        if r['status']=='done':break
        if r['status']!='running':raise RuntimeError('Source failed: '+r['status'])
    time.sleep(10)
else:raise RuntimeError('Source worker did not finish within two hours')
cmd=[sys.executable,str(HERE/'ladder.py'),'--geometry','infinite','--L','2','--seed','17','--gauge','--chis','96','--sweeps','160','--load',r['checkpoint']]
print('START registered chi96 extension',flush=True)
raise SystemExit(subprocess.call(cmd))
