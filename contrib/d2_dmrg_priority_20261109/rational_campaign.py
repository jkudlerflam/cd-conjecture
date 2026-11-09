"""Independent Davidson tests of the rational simplification of C181 phi=1.53.
Repeat the literal uniform Pauli bond; only the global parity is selected.
"""
import argparse,json,os,subprocess,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent
CP=dict(Jx=2001/4000,Jy=-1999/4000,D=1/160,Jz=10,hz=-199759/10000,Kxz=0.)

def run(L,parity,stage):
    tag='rational_'+parity
    chis=[16,32,32,32] if stage=='validate' else ([16,32,64,96,96,96] if L==32 else [16,32,64,96,96,128,128])
    cmd=[sys.executable,str(HERE/'full_chain.py'),tag,str(L),'--cp',json.dumps(CP),'--chis',*map(str,chis),'--res','1e-8','--parity',parity,'--nocorr','--outdir',str(HERE/'results')]
    with (HERE/'logs'/f'{tag}_L{L}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    path=HERE/'results'/f'{tag}_{L}.json';d=json.loads(path.read_text())
    if stage=='validate':subprocess.run([sys.executable,str(HERE/'check_ed.py'),str(path),'--output',str(HERE/'results'/f'ED_{tag}_{L}.json')],check=True)
    print({k:d[k] for k in ['tag','L','E','Sbar','c_prof','n_mid','N','P_total','Ymid','trunc','time_s']},flush=True)

def main():
    a=argparse.ArgumentParser();a.add_argument('--stage',choices=['validate','production'],default='validate');a.add_argument('--workers',type=int,default=2);args=a.parse_args()
    plan=[(L,p,args.stage) for L in ([10] if args.stage=='validate' else [32,64]) for p in ['even','odd']]
    (HERE/f'rational_{args.stage}_plan.json').write_text(json.dumps(dict(couplings=CP,stage=args.stage,tasks=plan,criterion='independent parity validation, entropy/energy convergence; no CFT claim from short chains'),indent=2)+'\n')
    with ThreadPoolExecutor(max_workers=args.workers) as pool:list(pool.map(lambda z:run(*z),plan))
if __name__=='__main__':main()
