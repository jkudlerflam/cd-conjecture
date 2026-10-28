"""Final C19 variant: both hopping and chemical fields from exchange."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
import scipy.linalg as la
from uniform_cell import *
from reflection_chain import run
from completion_checks import third_order

out=dict(algebra=[],ed=[],third_order=[])
for n in [1,2]:
    c=reflection_cell(n,True)
    err=float(la.norm(direct_second_order(c)-c['H2']))
    parity_error=float(la.norm(c['H']*c['parity'][None,:]-c['parity'][:,None]*c['H']))
    rec=dict(n=n,D=c['D'],rho=float(c['rho']),bond_error=err,parity_error=parity_error,
             kernel_error=float(la.norm(c['H']@c['G'])),reflection_error=c['reflection_error'],
             min_weight=float(min(la.eigvalsh(m).min() for m in c['Ms'])))
    assert err<1e-11 and parity_error<1e-12 and rec['min_weight']>0
    out['algebra'].append(rec);print(json.dumps(rec),flush=True)
    third=third_order(c);third['n']=n;out['third_order'].append(third)
    print(json.dumps(third),flush=True)
    out['ed']+=run(c,2,[.04,.02] if n==1 else [.04])
    if n==1:out['ed']+=run(c,3,[.04])
    with open('exchange_field_results.json','w') as f:json.dump(out,f,indent=2)
print('PASS: exchange-only fields and hopping',flush=True)
