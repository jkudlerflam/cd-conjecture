"""Reproducible bounded full-chain DMRG campaigns; no cluster use."""
import argparse,json,os,subprocess,sys,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent
CP=dict(Jx=.4979983,Jy=-.5020017,D=.0152122,Jz=4.,hz=-7.9463,Kxz=.13693063937629152)

def run(tag,L,cp,chis,seed=0):
    cmd=[sys.executable,str(HERE/'full_chain.py'),tag,str(L),'--cp',json.dumps(cp),'--chis',*map(str,chis),'--res','1e-9','--seed',str(seed),'--nocorr','--outdir',str(HERE/'results')]
    log=HERE/'logs'/f'{tag}_L{L}.log'
    with log.open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    d=json.loads((HERE/'results'/f'{tag}_{L}.json').read_text())
    print({k:d[k] for k in ['tag','L','E','Sbar','c_prof','n_mid','N','Xmid','Ymid','trunc','time_s']},flush=True)
    return d

def main():
    a=argparse.ArgumentParser();a.add_argument('--stage',choices=['ladder','sizes','realhop'],default='ladder');a.add_argument('--workers',type=int,default=1);args=a.parse_args()
    if args.stage=='ladder':
        fields=[-7.9500,-7.9480,-7.9470,-7.9463,-7.9455,-7.9450,-7.9440,-7.9420,-7.9400]
        plan=[dict(tag='lad_'+str(h).replace('-','m').replace('.','p'),L=24,cp=dict(CP,hz=h),chis=[16,32,48,48,64,64]) for h in fields]
    elif args.stage=='sizes':
        plan=[dict(tag='size_'+str(h).replace('-','m').replace('.','p'),L=48,cp=dict(CP,hz=h),chis=[16,32,48,64,64,96,96]) for h in [-7.9450,-7.9440]]
    else:
        sums=[-.0065,-.006,-.0055,-.00525,-.005,-.0045,-.0040034,-.0035]
        plan=[dict(tag='hop_'+str(s).replace('-','m').replace('.','p'),L=24,cp=dict(CP,hz=-7.9463,Jx=(1+s)/2,Jy=(s-1)/2),chis=[16,32,48,48,64,64]) for s in sums]
    (HERE/(args.stage+'_plan.json')).write_text(json.dumps(plan,indent=2)+'\n')
    if args.workers==1:
        for x in plan:run(**x)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(lambda x:run(**x),plan))
if __name__=='__main__':main()
