import numpy as np,json
from scipy.optimize import root
C=32.;tp=10/9

def quantities(E,s):
 D=np.sqrt(E*E-16+16*s*s);A=E-D
 D1=16*s/D;D2=16/D-256*s*s/D**3
 se=C*s*s/A
 se1=C*(2*s/A+s*s*D1/A**2)
 se2=C*(2/A+4*s*D1/A**2+s*s*D2/A**2+2*s*s*D1**2/A**3)
 return se,se1,se2

def fun(x):
 E,s=x;se,se1,se2=quantities(E,s)
 return [8*tp*s+se1,8*tp+se2]
r=root(fun,[-3.71,.49],tol=1e-11)
E,s=r.x;se,se1,se2=quantities(E,s)
U=E+2*tp-4*tp*s*s-se
out=dict(success=bool(r.success),E=float(E),s=float(s),Q=float(2*np.arcsin(s)),U=float(U),delta=float(U+14/9),residual=np.linalg.norm(r.fun))
print(json.dumps(out,indent=2))
with open('pair_saddle_results.json','w') as f:json.dump(out,f,indent=2)
