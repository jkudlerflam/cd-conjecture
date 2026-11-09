"""Summarize complete saved outputs without converting finite-size fits into claims."""
import json,math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def read(p):return json.loads(p.read_text())
def slopes(d):
    L=d['L'];S=np.array(d['S']);l=np.arange(1,L);x=np.log((2*L/np.pi)*np.sin(np.pi*l/L))/6
    out={}
    for f in (1/8,1/4,1/3,3/8):
        w=(l>=math.ceil(f*L))&(l<=math.floor((1-f)*L))
        A=np.stack([x[w],np.ones(w.sum()),(-1.)**l[w]],axis=1)
        c=np.linalg.lstsq(A,S[w],rcond=None)[0][0]
        out[str(f)]=float(c)
    return out

def row(d):
    return dict(hz=d['couplings']['hz'],L=d['L'],E=d['E'],Sbar=d['Sbar'],n_mid=d['n_mid'],N=d['N'],
                Xmid=d['Xmid'],Ymid=d['Ymid'],c_prof=d['c_prof'],window_fits=slopes(d),trunc=d['trunc'],
                final_residual=d['stats'][-1]['res_max'],energy_change_last=d['hist'][-1][1]-d['hist'][-2][1],seconds=d['time_s'])

def main():
    ladders=sorted([read(p) for p in (HERE/'results').glob('lad_*.json')],key=lambda d:d['couplings']['hz'])
    sizes=sorted([read(p) for p in (HERE/'results').glob('size_*.json')],key=lambda d:d['couplings']['hz'])
    out=dict(ladder=[row(d) for d in ladders],sizes=[row(d) for d in sizes],two_size=[],controls=[])
    for d in sizes:
        match=[a for a in ladders if a['couplings']==d['couplings']]
        if match:
            a=match[0];out['two_size'].append(dict(hz=d['couplings']['hz'],lengths=[a['L'],d['L']],
                central_three_c=6*(d['Sbar']-a['Sbar'])/math.log(d['L']/a['L'])))
    for kind in ('xx','leg2'):
        ps=[HERE/'results'/f'control_{kind}_{L}.json' for L in (24,48)]
        if all(p.exists() for p in ps):
            ds=list(map(read,ps));out['controls'].append(dict(kind=kind,c_two=6*(ds[1]['Sbar']-ds[0]['Sbar'])/math.log(2),
                exact_c_two=6*(ds[1]['exact_Sbar']-ds[0]['exact_Sbar'])/math.log(2),max_entropy_error=max(d['max_entropy_error'] for d in ds),
                max_energy_error=max(abs(d['energy_error']) for d in ds)))
    (HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['two_size'],indent=2));print(json.dumps(out['controls'],indent=2))
if __name__=='__main__':main()
