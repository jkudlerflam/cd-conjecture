from ks_test import *
rows=[]
for N in range(8,25,2):
 o=Orbits(N,'even');d=density(N,'ising',.35);h=total(d,N)
 H=o.matrix(h,0);e,g=lowest(H,1);Ht=o.matrix(h,2); et,v=lowest(Ht,10)
 hv=o.matrix(total(d,N,2),0,2)@g[:,0];w=abs(v.conj().T@hv)**2
 row=dict(N=N,gaps=(et-e[0]).tolist(),weights=w.tolist(),candidates=(8*w/(et-e[0])**2).tolist())
 rows.append(row);json.dump(rows,open('stress_tracking.json','w'),indent=2)
 print(N,row['candidates'][:2],flush=True)
