import json,numpy as np
from pathlib import Path
from d3 import cases
jobs=[]
def add(family,L,coeff,label,**kw):jobs.append(dict(family=family,L=L,coeff=list(map(float,coeff)),label=label,levels=12,**kw))
for L in [14,16]:
 for t in [.015,np.pi/48,np.pi/24,.16,.1963495408,1.308996939,2.0943951024]:
  add('sp4',L,[np.cos(t),np.sin(t)],f'sp4_theta{t:.11f}',cartan=True,entropy=L<=14)
for L in [12,14]:
 for i in [70,97,105,110,123,131,139,144,157,160,173,178,186]:
  label,c=cases('su2',192)[i];add('su2',L,c,label,entropy=L<=12)
for L in [14,16]:
 for x,y in [(0,.25),(.125,.5),(.125,.375),(.25,.25),(.25,.5),(.375,.375)]:
  add('orbital',L,[1,y,x],f'orb_x{x}_y{y}',entropy=L<=14)
for L in [16,18]:
 for i in [12,16,20,23,84]:
  label,c=cases('spin1',96)[i];add('spin1',L,c,label,entropy=L<=16)
# U(1) candidates must be compared over all charges. Separate array jobs permit this.
for L in [16]:
 for delta,beta,D in [(0,1,2),(.5,1,2),(1,1,2),(1.5,2,1),(1.5,2,2)]:
  for M in range(L+1):
   add('spin1',L,[1,beta,delta-1,D],f'u1_z{delta}_b{beta}_D{D}_M{M}',charge=[M],entropy=M==0)
for family,lengths in [('uls3',[15,18]),('uls4',[16]),('tb1',[16,18]),('tb32',[12,14])]:
 for L in lengths:add(family,L,[1],family+'_calibration',entropy=L<=16)
Path('cluster_manifest.json').write_text(json.dumps(jobs,indent=2)+'\n')
print(f'{len(jobs)} tasks. Submit: sbatch --array=0-{len(jobs)-1}%4 [partition/account options] typhon.sbatch')
