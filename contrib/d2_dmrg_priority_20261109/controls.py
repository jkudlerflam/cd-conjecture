"""Same Davidson solver; XX and two-XX-leg controls versus exact correlations."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(k,'2')
import argparse,json,sys,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'codex/work_C183'))
import dmrg2 as D
import dmrg_c as DC
from ctrl183 import xx_bond,leg2_bond,HZ,SZ,I2


def exact_xx(L):
    h=np.diag(np.full(L,-2*HZ))+np.diag(np.full(L-1,2.),1)+np.diag(np.full(L-1,2.),-1)
    e,U=np.linalg.eigh(h);C=U[:,e<0]@U[:,e<0].T;S=[]
    for l in range(1,L):
        x=np.linalg.eigvalsh(C[:l,:l]);x=x[(x>1e-14)&(x<1-1e-14)]
        S.append(float(-np.sum(x*np.log(x)+(1-x)*np.log(1-x))))
    return HZ*L+float(np.sum(e[e<0])),np.array(S),int(np.sum(e<0))


def main():
    a=argparse.ArgumentParser();a.add_argument('--lengths',type=int,nargs='+',default=[24,48,64]);a.add_argument('--chis',type=int,nargs='+',default=[16,32,64,96,96]);args=a.parse_args()
    for kind in ('xx','leg2'):
        for L in args.lengths:
            sites=L if kind=='xx' else L//2;fac=1 if kind=='xx' else 2
            he=xx_bond() if kind=='xx' else leg2_bond()
            edge=.5*HZ*(SZ if kind=='xx' else np.kron(SZ,I2)+np.kron(I2,SZ))
            ee,se,ne=exact_xx(sites);ee*=fac;se*=fac
            t=time.time();r=D.dmrg(he,sites,chis=args.chis,cutoff=1e-12,res_tol=1e-10,edge_field=edge,verbose=True)
            S=np.array(r['S']);c=sites//2;fits=DC.profile_fits(S,sites)
            out=dict(kind=kind,L=L,sites=sites,HZ=float(HZ),chis=args.chis,E=r['E'],exact_E=ee,energy_error=r['E']-ee,
                     exact_N=fac*ne,S=S.tolist(),exact_S=se.tolist(),max_entropy_error=float(np.max(abs(S-se))),
                     Sbar=float(np.mean(S[c-2:c+1])),exact_Sbar=float(np.mean(se[c-2:c+1])),c_prof=fits[0],c_sm=fits[1],
                     hist=r['hist'],stats=r['stats'],time_s=time.time()-t)
            (HERE/'results'/f'control_{kind}_{L}.json').write_text(json.dumps(out,indent=2)+'\n')
            print({k:v for k,v in out.items() if k not in ['S','exact_S','hist','stats']},flush=True)
if __name__=='__main__':main()
