"""Independent small checks and formatted numerical summary."""
import json,math
import numpy as np

def creation(a,m):
    M=np.zeros((2**m,2**m))
    for s in range(2**m):
        if not (s>>a)&1:
            M[s|(1<<a),s]=(-1)**((s&((1<<a)-1)).bit_count())
    return M

C=[creation(a,4) for a in range(4)];n=[c@c.T for c in C];I=np.eye(16)
Da=I-n[0]-n[2]+2*n[0]@n[2];Db=I-n[1]-n[3]+2*n[1]@n[3]
Ta=C[0]@(I-2*n[1])@C[2].T
Tb=C[1]@(I-2*n[2])@C[3].T
Qa=Da+Ta+Ta.T;Qb=Db+Tb+Tb.T
P=np.zeros((16,16))
for s in range(16):P[((s&3)<<2)|(s>>2),s]=1
assert np.linalg.norm(Qa@Qb-P)<1e-13
assert np.linalg.norm(Qa@Qb-Qb@Qa)<1e-13
with open('ed_results.json') as f: ed=json.load(f)
with open('extra_results.json') as f: extra=json.load(f)
summary={'swap_fermionic_identity_error':float(np.linalg.norm(Qa@Qb-P)),'half_chain':[]}
for kind,d,L in [('swap',3,6),('swap',4,4),('tj',3,6)]:
    a=next(r for r in ed['periodic'] if (r['model'],r['d'],r['L'])==(kind,d,L))
    b=next(r for r in ed['periodic'] if (r['model'],r['d'],r['L'])==(kind,d,2*L))
    summary['half_chain'].append({'model':kind,'d':d,'sizes':[L,2*L],'c':3*(b['entropy'][str(L)]-a['entropy'][str(L//2)] )/np.log(2)})
row=next(r for r in extra['hubbard'] if r['L']==8 and r['U']==0)
errs=[]
for ell in range(1,5):
    x=np.arange(ell);ks=np.array([-np.pi/8,np.pi/8])
    corr=sum(np.exp(1j*k*(x[:,None]-x[None,:])) for k in ks)/8
    p=np.linalg.eigvalsh(corr);p=p[(p>1e-14)&(p<1-1e-14)]
    entropy=-2*np.sum(p*np.log(p)+(1-p)*np.log(1-p))
    errs.append(abs(entropy-row['entropy'][str(ell)]))
summary['free_hubbard_L8_entropy_vs_correlation_max_error']=float(max(errs))
assert max(errs)<1e-11
# Shell labels and intentionally unused entropy fits.
for r in extra['hubbard']:
    r['ground_degenerate']=bool(r['energies'][1]-r['energies'][0]<1e-8)
    r['use_for_c_extrapolation']=False
    r['note']='Spectrum check only; no thermodynamic c inferred from these shell-dependent sizes.'
with open('extra_results.json','w') as f:json.dump(extra,f,indent=2)
with open('audit_results.json','w') as f:json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2))
