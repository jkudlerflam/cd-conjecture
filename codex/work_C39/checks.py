"""Full code-space bond extraction and source/Gram checks."""
import json,numpy as np, scipy.linalg as la
import model,pencil,certify

def normal_quadratic(K,n):
 N=2**(2*n);NN=2**n;local=model.fermions(n)
 P=np.diag([(-1)**b.bit_count() for b in range(NN)])
 C=[np.kron(c,np.eye(NN)) for c in local]+[np.kron(P,c) for c in local]
 ids=[(1<<j)*NN for j in range(n)]+[1<<j for j in range(n)]
 M=K[np.ix_(ids,ids)]-K[0,0]*np.eye(2*n)
 quad=K[0,0]*np.eye(N,dtype=complex)
 for i in range(2*n):
  for j in range(2*n):quad+=M[i,j]*(C[i].T@C[j])
 return quad,M

def run():
 out=[]
 for n in [1,2,3]:
  c=model.cell(n,a=3.5);K=model.exact_bond(c);N=c['N']
  np.savez_compressed(f'bond_n{n}.npz',K2=K,field=c['G'].T@(c['R'][:,None]*c['G']),H3_coefficient=np.kron(K,np.eye(N))+np.kron(np.eye(N),K))
  charges=np.array([b.bit_count() for b in range(N*N)])
  T=model.hopping(K,n);quad,M=normal_quadratic(K,n)
  row=dict(n=n,code_bond_dimension=N*N,code_three_dimension=N**3,source_cell_dimension=3*(2*n+2)*N,scalar=float(K[0,0].real),hopping_real=T.real.tolist(),hopping_imag=T.imag.tolist(),hermiticity=float(la.norm(K-K.conj().T)),charge_commutator=float(la.norm((charges[:,None]-charges[None,:])*K)),normal_order_nonGaussian=float(la.norm(K-quad,2)),gram_checks=[])
  # Exact rational scalar inverse versus independent Lyapunov solve.
  GF=certify.gram(n,certify.F(1,4));GN=np.array(GF,float)
  row['rational_gram_error']=float(la.norm(GN-pencil.gram(n,.25)))
  # Direct pair-resolvent g^2 coefficient, no finite differences in g.
  for D in [10.,100.,1000.]:
   HH=model.k2(n,.25,D,0,np.pi/2,True);TT=model.hopping(HH,n)
   approx=-D**(2*n+4)*TT/2
   row['gram_checks'].append(dict(D=D,relative_error=float(la.norm(approx-GN)/la.norm(GN))))
  out.append(row);print(json.dumps(row,indent=2),flush=True)
 with open('bond_checks.json','w') as f:json.dump(out,f,indent=2)
if __name__=='__main__':run()
