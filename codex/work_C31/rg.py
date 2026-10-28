"""One-loop number-conserving Cooper/forward flow for centered distinct seas.
Conventions and the shell integrals are spelled out in RESULT.md.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from scipy.integrate import solve_ivp
from itertools import combinations
rows=json.load(open('centered_scan_results.json'));out=[]
for row in rows:
 n=row['n'];v=np.array(row['v']);g=np.array(row['rows'][0]['g2']);f=g/(np.pi*(v[:,None]+v[None,:]))
 pairs=list(combinations(range(n),2));y=np.zeros((n,n))
 for r in row['rows'][0]['pairs']:y[r['a'],r['b']]=y[r['b'],r['a']]=r['y']
 def unpack(z):return z[:n*n].reshape(n,n),z[n*n:].reshape(n,n)
 def dimensions(f):return np.diag(f)[:,None]+np.diag(f)[None,:]-2*f
 def fun(t,z):
  f,y=unpack(z);x=dimensions(f)
  df=y*y*(v[:,None]*v[None,:])/(v[:,None]+v[None,:])**2
  np.fill_diagonal(df,-np.sum(y*y,axis=1)/4)
  dy=-x*y-.5*(y@y);np.fill_diagonal(dy,0.)
  return np.concatenate([df.ravel(),dy.ravel()])
 z0=np.concatenate([f.ravel(),y.ravel()]);s=solve_ivp(fun,[0,3],z0,method='Radau',rtol=2e-10,atol=1e-12)
 ff,yy=unpack(s.y[:,-1]);x0=dimensions(f);xf=dimensions(ff)
 rec=dict(n=n,success=s.success,steps=len(s.t),min_x_initial=min(x0[a,b] for a,b in pairs),min_x_final=min(xf[a,b] for a,b in pairs),max_x_loss=max(x0[a,b]-xf[a,b] for a,b in pairs),max_y_initial=np.max(abs(y)),max_y_final=np.max(abs(yy)),minimum_x_along=min(dimensions(unpack(z)[0])[a,b] for z in s.y.T for a,b in pairs))
 out.append(rec);print(json.dumps(rec))
json.dump(out,open('rg_results.json','w'),indent=2)
