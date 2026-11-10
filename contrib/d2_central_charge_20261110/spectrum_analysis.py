"""Conditional Casimir/velocity comparison, including identical-size XX controls."""
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
out=[]
for name in ['candidate','xx']:
    rows=json.loads((HERE/f'spectrum_{name}.json').read_text())['rows']
    for n in [3,4]:
        rr=rows[-n:];L=np.array([r['L'] for r in rr],float)
        E=np.array([r['E0'] for r in rr]);v=np.array([np.mean(r['velocities']) for r in rr])
        for degree in [1,2]:
            if len(rr)<=degree:continue
            ef=np.polyfit(1/L**2,E/L,degree);vf=np.polyfit(1/L**2,v,degree)
            cv=-6/np.pi*ef[-2];vinf=vf[-1]
            out.append(dict(model=name,sizes=L.astype(int).tolist(),polynomial_degree=degree,
                e_infinite=ef[-1],cv=cv,velocity_infinite=vinf,c_conditional=cv/vinf,
                e_fit_rms=float(np.sqrt(np.mean((np.polyval(ef,1/L**2)-E/L)**2))),
                velocity_fit_rms=float(np.sqrt(np.mean((np.polyval(vf,1/L**2)-v)**2)))))
(HERE/'spectrum_analysis.json').write_text(json.dumps(dict(rows=out,
    caveat='Assumes the lowest momentum ±2pi/L states are the relevant sound descendants and one velocity exhausts the Casimir term. Several sectors/velocities or long crossovers invalidate this identification. Fit-window/degree scatter is systematic, not a confidence interval.'),indent=2)+'\n')
for r in out:print(r)
