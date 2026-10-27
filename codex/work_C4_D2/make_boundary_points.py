"""Generate targeted critical-line and multicritical-neighborhood grids."""
import json
from d2 import normalized
import numpy as np
points=[]
def add(family,J,field,label):
    J,field=normalized(J,field)
    points.append(dict(id='b'+str(len(points)),family=family,label=label,J=J,field=field))
for gamma in [0,.03,.1,.3,.6,1.]:
    for hz in [1.6,1.8,1.95,2.,2.05,2.2]:
        add('XY_Ising_saturation',np.diag([1+gamma,1-gamma,0.]),[0,0,hz],f'gamma={gamma},hz={hz}')
for delta in [-1.1,-1.02,-1.,-.999,-.99,-.95,-.9,-.8,0.,.8,1.,1.05,1.2]:
    for hx in [0.,.03,.12]:
        add('XXZ_FM_AF_boundaries',np.diag([1.,1.,delta]),[hx,0,0],f'delta={delta},hx={hx}')
for delta in [0.,.3,-.5]:
    for D in [.4,.55,.6,.65,.8,1.2]:
        for hz in [0.,1.,2.]:
            J=np.diag([1.6,.4,delta]); J[0,1]=D; J[1,0]=-D
            add('DM_pairing_threshold',J,[0,0,hz],f'delta={delta},D={D},hz={hz}')
for gamma in [.2,.6]:
    for delta in [-.5,.3,1.]:
        for hz in [1.5,1.8,2.,2.2,2.5]:
            add('interacting_XYZ_near_Ising',np.diag([1+gamma,1-gamma,delta]),[0,0,hz],f'g={gamma},delta={delta},hz={hz}')
json.dump(points,open('results/boundary_points.json','w'),indent=2)
print(len(points))
