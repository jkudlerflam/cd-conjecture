from ks_test import *
from free_ising_check import calc
checks=[]
N=8;s=np.arange(1<<N,dtype=np.uint32)
def full(terms):
 rows=[];cols=[];vals=[]
 for (x,z),c in terms.items():
  rows.extend(s^x);cols.extend(s);vals.extend(c*(1-2*parity(s&z).astype(float)))
 return coo_matrix((vals,(rows,cols)),shape=(len(s),len(s))).tocsr()
for model,c,label in [('ising',0.,'TFIM'),('ising',.2,'OF0.2'),('ising',.35,'OF0.35'),('spin',0.,'Heisenberg'),('spin',.2,'J1J2')]:
 d=density(N,model,c);H=full(total(d,N));e,g=lowest(H,2);hv=full(total(d,N,2))@g[:,0]
 m0=np.vdot(hv,hv).real;m1=np.vdot(hv,H@hv-e[0]*hv).real
 R=run(N,model,c,label);a=R['centered'];o=Orbits(N,'even' if model=='ising' else 'sz0')
 Hk=o.matrix(total(d,N),1);err=Hk-Hk.conj().T
 check=dict(label=label,full_gap=float(e[1]-e[0]),E_error=abs(e[0]-R['E0']),m0_error=abs(m0-a['m0']),m1_error=abs(m1-a['m1']),hermiticity=float(np.max(abs(err.data),initial=0)))
 checks.append(check); print(check)
r=json.load(open('results_extended.json'))
for x in r:
 if x['label']=='TFIM':
  f=calc(x['N']);a=x['centered']
  assert max(abs(f[k]-a[k]) for k in ['c2','m0N2','m1N3','cMV'])<1e-7
json.dump(checks,open('validation.json','w'),indent=2)
print('All centered Ising sparse results match independent free-fermion sums to <1e-7 in the scaled moments.')
