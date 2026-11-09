"""Summarize independent rational-parity data and convergence limits in both sectors."""
import json,math
import numpy as np
from pathlib import Path
from analyse import slopes
HERE=Path(__file__).resolve().parent

def fixed_fraction_entropy(d,fraction):
    L=d['L'];S=np.asarray(d['S']);l=np.arange(1,len(S)+1)*L/(len(S)+1);lo=L*(1-fraction)/2;hi=L*(1+fraction)/2
    pts=np.r_[lo,l[(l>lo)&(l<hi)],hi];v=np.interp(pts,l,S)
    return float(np.sum((v[1:]+v[:-1])*np.diff(pts)/2)/(hi-lo))

def fraction_controls():
    from controls import exact_xx
    tables={}
    for mode,lengths in [('saved_dmrg',[24,48]),('exact',[32,64])]:
        rows=[]
        for kind in ('xx','leg2'):
            ds=[]
            for L in lengths:
                if mode=='saved_dmrg':
                    ds.append(json.loads((HERE/'results'/f'control_{kind}_{L}.json').read_text()))
                else:
                    sites=L if kind=='xx' else L//2;factor=1 if kind=='xx' else 2
                    _,S,N=exact_xx(sites);ds.append(dict(L=L,S=(factor*S).tolist(),exact_N=N*factor))
            for f in (.25,.5):
                means=[fixed_fraction_entropy(d,f) for d in ds]
                rows.append(dict(control=kind,fraction=f,lengths=lengths,N=[d['exact_N'] for d in ds],means=means,c=6*(means[1]-means[0])/math.log(lengths[1]/lengths[0])))
        tables[mode]=rows
    return tables

def main():
    out={'couplings':dict(Jx=2001/4000,Jy=-1999/4000,D=1/160,Jz=10,hz=-199759/10000,Kxz=0.),'sectors':[],'ground_branches':[]}
    for L in (10,32,64):
        branches=[]
        for parity in ('even','odd'):
            paths=[HERE/'results'/f'rational_{parity}_{L}.json',HERE/'results'/f'rational_{parity}_continued_{L}.json']
            files=[p for p in paths if p.exists()]
            if not files:continue
            p=files[-1];d=json.loads(p.read_text());h=d['hist'];stats=d['stats'][-1]
            q=dict(L=L,parity=parity,file=str(p.relative_to(HERE)),E=d['E'],E_excess=d['E']-((L-1)*d['couplings']['Jz']+L*d['couplings']['hz']),Sbar=d['Sbar'],n_mid=d['n_mid'],N=d['N'],P_total=d['P_total'],P_hop_mid=d['P_hop_mid'],c_prof=d['c_prof'],window_fits=slopes(d),chi=d['chis'][-1],trunc=d['trunc'],final_residual=stats['res_max'],last_sweep_energy_change=h[-1][1]-h[-2][1],seconds=d['time_s'])
            q['S_middle_quarter']=fixed_fraction_entropy(d,.25)
            q['S_middle_half']=fixed_fraction_entropy(d,.5)
            q['E_expectation']=sum(d['ebond'])+d['couplings']['hz']*(d['Z'][0]+d['Z'][-1])/2+d['couplings'].get('hx',0.)*(d['X'][0]+d['X'][-1])/2
            q['energy_expectation_minus_solver']=q['E_expectation']-q['E']
            q['local_residual_pass']=q['final_residual']<=1.01e-8
            if len(files)>1:
                old=json.loads(files[0].read_text());q['continued_minus_original_Sbar']=d['Sbar']-old['Sbar'];q['continued_minus_original_E']=d['E']-old['E']
            branches.append(q);out['sectors'].append(q)
        if len(branches)==2:
            branches.sort(key=lambda x:x['E_expectation']);best=dict(branches[0]);best['other_parity_gap']=branches[1]['E_expectation']-branches[0]['E_expectation'];out['ground_branches'].append(best)
    gs={d['L']:d for d in out['ground_branches']}
    if 32 in gs and 64 in gs:
        out['c32_to_64']=6*(gs[64]['Sbar']-gs[32]['Sbar'])/math.log(2)
        out['c_middle_quarter_32_to_64']=6*(gs[64]['S_middle_quarter']-gs[32]['S_middle_quarter'])/math.log(2)
        out['c_middle_half_32_to_64']=6*(gs[64]['S_middle_half']-gs[32]['S_middle_half'])/math.log(2)
        out['both_residual_pass']=all(gs[L]['local_residual_pass'] for L in (32,64))
    out['fixed_fraction_controls']=fraction_controls()
    (HERE/'rational_summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
