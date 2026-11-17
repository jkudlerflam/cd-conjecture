#!/usr/bin/env python3
"""Convergence gates, exact controls and all-parity size comparisons."""
import json
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE,Ana


def row(d,file):
    h=d['hist'];last=h[-1];cfg=d['cfg'];L=cfg['L']
    final=[r for r in h if r['stage']==last['stage']]
    prev=[r for r in h if r['chi']<last['chi'] and 'Sbar6' in r]
    if not prev and cfg.get('from'):
        p=HERE/'outputs'/(cfg['from']+'.json')
        if p.exists():prev=[r for r in json.loads(p.read_text())['hist'] if r['chi']<last['chi'] and 'Sbar6'in r]
    de=abs(final[-1]['E']-final[-2]['E']) if len(final)>1 else None
    ds=max(abs(final[-1][k]-final[-2][k]) for k in ('Sbar1','Sbar6','Sbar12')) if len(final)>1 else None
    dn=abs(final[-1]['n_mid']-final[-2]['n_mid']) if len(final)>1 else None
    dc=abs(last['Sbar6']-prev[-1]['Sbar6']) if prev else None
    gates=dict(residual=last['res_max']<=1e-8,stationarity=de is not None and de<=1e-8 and ds<=1e-3 and dn<=5e-4,
               discarded_weight=last['trunc']<=1e-7,chi_step=dc is not None and dc<=.01)
    r=dict(file=file,tag=cfg['tag'],L=L,P=cfg['P'],V=cfg.get('V'),a=cfg.get('a'),kappa=cfg.get('kappa',1),
        legs=cfg.get('legs'),chi=last['chi'],E=last['E_phys'],n_mid=last['n_mid'],Sbar=d['Sbar'],
        c_profile=d['c_profile'],gates=gates,all_gates=all(gates.values()),dE=de,dS=ds,dn=dn,dS_chi=dc,
        res_max=last['res_max'],trunc=last['trunc'],runtime_seconds=d['runtime_seconds'],peak_RSS_GB=d['peak_RSS_GB'],
        valley=d['valley'],pair=d['pair'])
    if cfg.get('legs'):
        ex=Ana.ff_entropies(L,d['bond']['mu'],cfg['legs'])
        r['exact_control']=dict(Sbar={str(w):Ana.sbar(ex,L,w) for w in (1,6,12)},
            max_profile_error=float(np.max(abs(np.array(d['S'])-ex))),
            Sbar_errors={str(w):d['Sbar'][str(w)]-Ana.sbar(ex,L,w) for w in (1,6,12)})
    return r


def two(u,v):
    return {str(w):6*(v['Sbar'][str(w)]-u['Sbar'][str(w)])/np.log(v['L']/u['L']) for w in (1,6,12)}


def main():
    rows=[]
    with threadpool_limits(limits=1):
        for f in sorted((HERE/'outputs').glob('*.json')):
            d=json.loads(f.read_text())
            if d.get('status')=='done':rows.append(row(d,str(f.relative_to(HERE))))
    groups={}
    for r in rows:
        key=(r['V'],r['a'],r['kappa'],r['legs'],r['L'],r['P'])
        if key not in groups or r['chi']>groups[key]['chi'] or (r['chi']==groups[key]['chi'] and r['E']<groups[key]['E']):groups[key]=r
    fams=sorted({k[:4] for k in groups},key=str);pairs=[];parity=[]
    for fam in fams:
        ls=sorted({k[4] for k in groups if k[:4]==fam})
        for L in ls:
            sec=[r for k,r in groups.items() if k[:5]==fam+(L,)]
            if len(sec)==2:
                p0=next(r for r in sec if r['P']==0);p1=next(r for r in sec if r['P']==1);t2=1/(2*fam[0])
                parity.append(dict(family=fam,L=L,gap=p1['E']-p0['E'],L_gap_over_t2=L*(p1['E']-p0['E'])/t2,
                                   sectors=[r['tag'] for r in sec],both_gates=all(r['all_gates'] for r in sec)))
            if 2*L not in ls:continue
            sec2=[r for k,r in groups.items() if k[:5]==fam+(2*L,)]
            small=min(sec,key=lambda r:r['E']);big=min(sec2,key=lambda r:r['E'])
            pair=dict(family=fam,lengths=[L,2*L],ground_tags=[small['tag'],big['tag']],ground_parities=[small['P'],big['P']],
                c_two=two(small,big),ground_gates=small['all_gates'] and big['all_gates'],
                densities=[small['n_mid'],big['n_mid']],all_parity_pairings=[dict(P=[x['P'],y['P']],c_two=two(x,y),gates=x['all_gates'] and y['all_gates']) for x in sec for y in sec2])
            if small['legs']:
                exa=dict(L=L,Sbar=small['exact_control']['Sbar']);exb=dict(L=2*L,Sbar=big['exact_control']['Sbar']);pair['exact_c_two']=two(exa,exb)
            elif len(sec)==len(sec2)==2:
                target=(small['n_mid']+big['n_mid'])/2
                def interp(ss):
                    x,y=sorted(ss,key=lambda r:r['n_mid']);z=(target-x['n_mid'])/(y['n_mid']-x['n_mid'])
                    return dict(L=x['L'],Sbar={str(w):x['Sbar'][str(w)]+z*(y['Sbar'][str(w)]-x['Sbar'][str(w)]) for w in (1,6,12)}),z
                ia,za=interp(sec);ib,zb=interp(sec2)
                pair['density_corrected']=dict(target=target,weights=[za,zb],extrapolates=not(0<=za<=1 and 0<=zb<=1),c_two=two(ia,ib))
            pairs.append(pair)
    out=dict(rows=rows,size_pairs=pairs,parity=parity,
             caution='Finite-size diagnostics only; these gates do not establish a thermodynamic central charge.')
    (HERE/'analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in rows:print(r['tag'],'n',round(r['n_mid'],5),'S6',round(r['Sbar']['6'],6),'gates',r['gates'],'dS_chi',r['dS_chi'])
    for p in pairs:print('pair',p['family'],p['lengths'],p['c_two'],'gates',p['ground_gates'],'n',p['densities'])


if __name__=='__main__':main()
