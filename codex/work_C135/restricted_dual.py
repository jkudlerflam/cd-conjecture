from core import *
from certificates import embed
def run(h,n,sym):
    D=2**n;m=n-1
    ops=[embed(h,j,m) for j in range(m-1)]
    x=cp.Variable(m-1);t=cp.Variable()
    K=np.kron(h,np.eye(D//4)).real+sum(x[j]*(np.kron(o,I)-np.kron(I,o)).real for j,o in enumerate(ops))
    keys=np.array([i.bit_count() if sym=='u1' else i.bit_count()%2 for i in range(D)])
    cons=[]
    for key in sorted(set(keys)):
        inds=np.flatnonzero(keys==key);cons.append(K[np.ix_(inds,inds)]-t*np.eye(len(inds))>>0)
    p=cp.Problem(cp.Maximize(t),cons);p.solve(solver='SCS',eps=1e-8,max_iters=20000)
    k=np.kron(h,np.eye(D//4))+sum(x.value[j]*(np.kron(o,I)-np.kron(I,o)) for j,o in enumerate(ops))
    return dict(n=n,value=float(p.value),lower=float(la.eigvalsh(k,subset_by_index=[0,0])[0]),status=p.status,coefficients=x.value.tolist(),scs={k:p.solver_stats.extra_stats['info'][k] for k in ['res_pri','res_dual','gap','status']})
out=[]
for name,h,sym in [('XX',xxz(0),'u1'),('Heisenberg',xxz(1),'u1'),('TFIM_1',tfim(1),'parity')]:
    for n in range(4,9):
        rec=run(h,n,sym);rec['name']=name;out.append(rec);write_json('restricted_dual.json',out);print(name,n,rec['value'],rec['status'],flush=True)
