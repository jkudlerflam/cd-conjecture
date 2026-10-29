"""Rebuild the C88 material certificate and finite checks locally."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,runpy,time
checks=[]
for script in ['run_certificate.py','selector_certificate.py','check_lift.py',
               'check_virtual.py','check_energy.py','check_frozen.py','check_sw.py',
               'check_sector.py']:
    start=time.monotonic()
    runpy.run_path(script,run_name='__main__')
    checks.append(dict(script=script,status='passed',seconds=time.monotonic()-start))
json.dump(dict(tag='numerical_and_integer_rational_checks',checks=checks),open('verification.json','w'),indent=2)
print(json.dumps(dict(all_checks_passed=True,checks=checks)),flush=True)
