"""C44 independent finite matrix and source consistency checks."""
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
from audit_numeric import core,buscore,contract,norm
import core_model
import wrapper_tools as wt

def run():
 out=[];zero=(np.array([0.]),np.array([1.]))
 for n in (1,2):
  c=core(n);K=np.zeros((4**n,4**n),complex)
  for bus,phase,weight in [(1,np.pi/2,1),(2,np.pi,.6)]:
   mm=wt.measure(buscore(c,bus));K+=weight*contract(mm,mm,zero,zero,phase)
  err=norm(K-core_model.exact_bond(c));assert err<1e-12
  partsmin=min(la.eigvalsh(x)[0] for x in [*c['parts'],c['ends']])
  residual=max(norm(x@c['G']) for x in [*c['parts'],c['ends']])
  R=c['G'].T@(c['R'][:,None]*c['G']);q=c['params']['q'];S=sum(q**(2*j) for j in range(n+2))
  expected=np.diag([sum(q**(2*j+2)/S*((bits>>j)&1) for j in range(n)) for bits in range(2**n)])
  readerr=norm(R-expected);assert partsmin>-1e-12 and residual<1e-12 and readerr<1e-12
  rec=dict(n=n,core_dimension=c['d'],literal_rule_check=c['literal_check'],smallest_local_eigenvalue=float(partsmin),termwise_code_residual=residual,bare_source_comparison_error=err,readout_compression_error=readerr)
  out.append(rec)
 # Direct inverse check at n=1, with every charge-core source state.
 c=core(1);direct=[]
 for p in (1,2):
  ch=wt.charge(3,p,.16);s=.03*wt.charge(3,2,.16)['gap'];bc=buscore(c,1);ex=bc['ex'];de=int(sum(ex));dc=ch['d'];T=ch['V']
  W=np.kron(np.eye(dc),c['ends'][np.ix_(ex,ex)])
  for j,B in enumerate(c['parts']):W+=np.kron(T.T@(ch['D'][:,j,None]*T),B[np.ix_(ex,ex)])
  H=W+np.diag(np.repeat(ch['ev']/s,de));e,A,res=wt.measure(bc,ch,s)
  errors=[]
  for side,which in [('L',0),('R',1)]:
   F=bc['F'][side]['0'];src=np.kron((T.T@(ch['occ'][which]*ch['chi']))[:,None],(F@c['G'])[ex])
   for shift in (.7,2.):
    raw=src.conj().T@la.solve(H+shift*np.eye(len(H)),src,assume_a='her')
    spec=A[side]['0'].conj().T@(A[side]['0']/(e+shift)[:,None]);errors.append(norm(raw-spec))
  direct.append(dict(p=p,dimension=len(H),maximum_shifted_resolvent_error=max(errors),graph_residual=res))
  assert max(errors)<1e-7
 Path('consistency_results.json').write_text(json.dumps(dict(core=out,direct=direct),indent=2)+'\n');print(json.dumps(dict(core=out,direct=direct),indent=2))
if __name__=='__main__':run()
