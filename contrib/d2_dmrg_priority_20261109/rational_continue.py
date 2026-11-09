"""Bounded continuation from saved rational-candidate MPS checkpoints."""
import argparse,json,subprocess,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
HERE=Path(__file__).resolve().parent

def run(task):
    L,p,chi=task
    old=HERE/'results'/f'rational_{p}_{L}.json';d=json.loads(old.read_text())
    tag=f'rational_{p}_continued'
    cmd=[sys.executable,str(HERE/'full_chain.py'),tag,str(L),'--cp',json.dumps(d['couplings']),'--chis',str(chi),str(chi),'--res','1e-8','--parity',p,'--nocorr','--resume',d['mps_checkpoint'],'--outdir',str(HERE/'results')]
    with (HERE/'logs'/f'{tag}_L{L}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    n=json.loads((HERE/'results'/f'{tag}_{L}.json').read_text())
    print(dict(L=L,parity=p,E=n['E'],Sbar=n['Sbar'],deltaE=n['E']-d['E'],deltaSbar=n['Sbar']-d['Sbar'],stats=n['stats']),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--L',type=int,required=True);a.add_argument('--chi',type=int,required=True);a.add_argument('--parities',nargs='+',default=['even','odd']);args=a.parse_args()
    plan=dict(L=args.L,chi=args.chi,parities=args.parities,residual=1e-8,sweeps=2)
    (HERE/f'rational_continue_{args.L}_{"_".join(args.parities)}_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[(args.L,p,args.chi) for p in args.parities]))
