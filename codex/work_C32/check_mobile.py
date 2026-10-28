"""Run bounded C32 source and mobility checks. Outputs local JSON."""
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import mobile as m

ROOT=Path(__file__).resolve().parent

def norm(A): return float(np.linalg.norm(A,2))
def tl(A): return A-np.trace(A)/len(A)*np.eye(len(A))

def run():
    start=time.monotonic(); out={'status':'numerical','cases':[]}
    for n in (1,2):
        c=m.core(n,g=.45,mass=.5,express=.2,loop=.08)
        G=c['G']; ee=eigh(c['H'],eigvals_only=True)
        bare=m.measure(c); K=m.bond(bare,bare); scale=norm(tl(K))
        rec={'n':n,'core_dimension':c['d'],'core_gap':float(ee[2**n]),
             'code_residual':norm(c['H']@G),'code_orthogonality':norm(G.conj().T@G-np.eye(2**n)),
             'charge_commutator':norm(c['Q'][:,None]*c['H']-c['H']*c['Q'][None,:]),
             'zero_first_order':max(norm(G.conj().T@F@G) for fs in c['F'].values() for F in fs.values()),
             'bare_traceless_norm':scale,'runs':[]}
        for w in (.16,.10):
            cs={p:m.charge(c['k'],p,w) for p in (1,2)}
            b=c['k']+1; phi=1/b
            comp=m.compensator(cs[1]['r'],cs[2]['r'],b,phi)
            avg={p:m.measure(m.flat_model(c,ch)) for p,ch in cs.items()}
            refs={(p,q):m.bond(avg[p],avg[q]) for p in cs for q in cs}
            beta=m.transport(cs[1],cs[2])
            ab=sum(np.sqrt(2/b)*np.sin(np.pi*j/b)*np.sqrt(2/(b+1))*np.sin(np.pi*j/(b+1)) for j in range(1,b))
            recw={'w':w,'charge_dimensions':[cs[p]['d'] for p in cs],
                  'charge_gaps':[float(cs[p]['gap']) for p in cs],
                  'endpoint_weights':[cs[p]['r'] for p in cs],
                  'beta_vacancy':float(beta),'beta_total':float(beta*ab**2*np.cos(phi/2)**(2*b)),
                  'compensator_norm':norm(comp),
                  'compensator_error':0.,'ratios':[]}
            for p,ch in cs.items():
                u=np.array([1,np.exp(-1j*(ch['N']-1)*phi)])/np.sqrt(2)
                recw['compensator_error']=max(recw['compensator_error'],float(abs(np.vdot(u,comp@u)*ch['r']-1)))
            for ratio in (.03,.01):
                s=ratio*min(cs[p]['gap'] for p in cs)
                ms={p:m.measure(c,ch,s) for p,ch in cs.items()}
                hs={p:m.helix(ch['N'],comp,phi,s) for p,ch in cs.items()}
                ks={(p,q):m.bond(ms[p],ms[q],hs[p],hs[q]) for p in cs for q in cs}
                rr={'s_over_charge_gap':ratio,'s':float(s),
                    'source_dimensions_with_helix':[len(ms[p][0])*len(hs[p][0]) for p in cs],
                    'max_graph_residual':max(ms[p][2] for p in cs),
                    'length_mismatch_full':max(norm(ks[p,q]-ks[1,1]) for p in cs for q in cs),
                    'length_mismatch_traceless_over_bare':max(norm(tl(ks[p,q]-ks[1,1]))/scale for p in cs for q in cs),
                    'fast_limit_error_full':max(norm(ks[p,q]-refs[p,q]) for p in cs for q in cs),
                    'fast_limit_error_traceless_over_bare':max(norm(tl(ks[p,q]-refs[p,q]))/scale for p in cs for q in cs),
                    'bare_error_full':max(norm(ks[p,q]-K) for p in cs for q in cs),
                    'scalar_density_interaction':float(np.trace(ks[2,2]-ks[2,1]-ks[1,2]+ks[1,1]).real/len(K))}
                recw['ratios'].append(rr)
                print('source',n,w,ratio,rr,flush=True)
            recw['averaged_length_mismatch_full']=max(norm(refs[p,q]-refs[1,1]) for p in cs for q in cs)
            recw['averaged_length_mismatch_traceless_over_bare']=max(norm(tl(refs[p,q]-refs[1,1]))/scale for p in cs for q in cs)
            rec['runs'].append(recw)
        out['cases'].append(rec)
        (ROOT/'mobile_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    out['seconds']=time.monotonic()-start
    (ROOT/'mobile_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('elapsed',out['seconds'],flush=True)

if __name__=='__main__': run()
