import json
import numpy as np
x=json.load(open('checks.json'));out={}
for name in dict.fromkeys(r['name'] for r in x['cases']):
    a=[r for r in x['cases'] if r['name']==name]
    L=np.array([r['L'] for r in a]);E=np.array([r['E0'] for r in a]);S=np.array([r['S_half'] for r in a])
    sf=np.linalg.lstsq(np.array([np.ones(len(L)),np.log(L),1/L**2]).T,S,rcond=None)[0]
    ef=np.linalg.lstsq(np.array([np.ones(len(L)),1/L**2,1/L**4]).T,E/L,rcond=None)[0]
    out[name]=dict(c_half_fit=float(3*sf[1]),c_step=float(3*(S[-1]-S[-2])/np.log(L[-1]/L[-2])),
                   c_velocity_weighted=float(-6*ef[1]/np.pi),max_residual=max(r['residual'] for r in a))
json.dump(out,open('summary.json','w'),indent=2)
print(json.dumps(out,indent=2))
