import json,numpy as np
from d3 import Sector,measure,append,cases
cs=dict(cases('spin1',96));cache={}
for line in open('charge_audit.jsonl'):
 row=json.loads(line);mins=row['charge_minima'];winner=min(mins,key=lambda r:r['E0'])
 if winner['E0']>=mins[0]['E0']-1e-8:continue
 L=row['L'];M=winner['M'];key=(L,M)
 if key not in cache:cache[key]=Sector('spin1',L,[M])
 result=measure(cache[key],cs[row['label']],row['label'],6)
 result['charge_audit']=mins;result['global_ground_M']=M
 append('corrected_u1_ground.jsonl',result)
