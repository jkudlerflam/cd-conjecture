#!/usr/bin/env python3
"""All-sector size, accuracy, starting-state and boundary-probe diagnostics."""
import argparse,json
import numpy as np
from threadpoolctl import threadpool_limits
from phase_common import HERE,OLD,write_json


def gates(raw,root=HERE):
    hist=raw['hist'];last=hist[-1];rr=[r for r in hist if r['stage']==last['stage']]
    prev=[r for r in hist if r['chi']<last['chi'] and 'Sbar6' in r]
    cfg=raw['cfg'];a,b=rr[-2:]
    de=abs(a['E']-b['E']);ds=max(abs(a[k]-b[k]) for k in ('Sbar1','Sbar6','Sbar12'));dn=abs(a['n_mid']-b['n_mid'])
    if not prev and cfg.get('parent_output'):
        p=root/cfg['parent_output'];prev=json.loads(p.read_text())['hist']
    dc=abs(last['Sbar6']-prev[-1]['Sbar6']) if prev else None
    checks=dict(residual=last['res_max']<=1e-8,stationarity=de<=1e-8 and ds<=1e-3 and dn<=5e-4,
        discarded_weight=last['trunc']<=1e-7,chi_step=dc is not None and dc<=.01)
    strict=de<2e-9 and ds<5e-4 and dn<2e-4 and last['res_max']<1e-8
    return dict(gates=checks,all_gates=all(checks.values()),dE=de,dS=ds,dn=dn,dS_chi=dc,
        res_max=last['res_max'],trunc=last['trunc'],stopping_rule_met=strict,
        schedule_limit_hit=len(rr)>=cfg['stages'][last['stage']]['max'] and not strict)


def summary(d,path):
    s=dict(file=str(path),tag=d['tag'],L=d['L'],P=d['P'],chi=d['chi'],boundary=d['boundary'],E=d['E'],operators={})
    for name,x in d['operators'].items():
        s['operators'][name]={w:{k:x[w][k] for k in ('M','mean','range','RMS','profile_peak','connected_peak')} for w in ('half','quarter')}
        C=np.array(x['connected']);M=len(C);rr=abs(np.arange(M)[:,None]-np.arange(M)[None,:]);sel=(rr>=M/3)&(rr<=M/2)
        s['operators'][name]['long_distance_RMS']=float(np.sqrt(np.mean(C[sel]**2)))
    return s


def compare(a,b):
    rows=[]
    for name in ('density','bond','kinetic'):
        for window in ('half','quarter'):
            for metric in ('RMS','connected_peak_weight'):
                x=a['operators'][name][window];y=b['operators'][name][window]
                u=x['RMS'] if metric=='RMS' else x['connected_peak']['connected_over_M']
                v=y['RMS'] if metric=='RMS' else y['connected_peak']['connected_over_M']
                tol=max(.05*max(abs(u),abs(v)),1e-4)
                rows.append(dict(operator=name,window=window,metric=metric,values=[u,v],change=abs(u-v),tolerance=tol,pass_gate=abs(u-v)<=tol))
    return dict(tags=[a['tag'],b['tag']],chis=[a['chi'],b['chi']],all_observable_gates=all(r['pass_gate'] for r in rows),rows=rows)


def common_q(d,q,name):
    o=d['operators'][name];ids=np.array(o['indices']);m=np.array(o['profile'])[ids];M=len(ids)
    p=np.exp(1j*q*(ids-ids[0]));C=np.array(o['connected'])
    sc=float(np.real(p.conj()@C@p)/M);sd=float(abs(p@(m-m.mean()))**2/M)
    return dict(q=q,S_connected=sc,S_disconnected=sd,connected_over_M=sc/M,disconnected_over_M=sd/M,total_over_M=(sc+sd)/M)


def supplementary(root,original,obs):
    channels=[];pairing=[];rings=[];controls=[]
    for r in original:
        key=r['file'].split('/')[-1];d=obs[key];M=d['operators']['density']['half']['M']
        anchors=[]
        for a in d['pair']:
            selected=[c for c in a['correlations'] if M/3<=c['j']-a['i']<=M/2]
            anchors.append(dict(i=a['i'],count=len(selected),**{name+'_RMS':
                float(np.sqrt(np.mean([abs(complex(*c[name]))**2 for c in selected])))
                for name in ('raw','disconnected','connected')}))
        pairing.append(dict(tag=r['tag'],L=r['L'],P=r['P'],chi=r['chi'],anchors=anchors,
            scope='Exploratory RMS at separations M/3 through M/2; explicit pair source means raw pair weight alone does not diagnose spontaneous order.'))
        p=root/'channels'/key
        if p.exists():
            c=json.loads(p.read_text())
            C=np.array(c['connected']);M=len(C);rr=abs(np.arange(M)[:,None]-np.arange(M)[None,:])
            far=(rr>=M/3)&(rr<=M/2)
            channels.append(dict(tag=r['tag'],L=r['L'],P=r['P'],chi=r['chi'],windows={w:
                {k:c[w][k] for k in ('M','mean','RMS','S_zero_over_M','uniform_order_weight','connected_peak_including_zero')}
                for w in ('half','quarter')},long_distance_RMS=float(np.sqrt(np.mean(C[far]**2)))))
    for p in sorted((root/'periodic').glob('*.json')):
        d=json.loads(p.read_text());rings.append(dict(L=d['L'],P=d['P'],energies=d['energies'],
            gaps_within_parity=d['gaps_within_parity'],momenta=d['momentum_if_translation_eigenstate'],
            peaks={k:d[k]['peak'] for k in ('density','bond','kinetic','valley_bond')},
            valley_uniform_order_weight=d['valley_bond']['uniform_order_weight'],
            runtime_seconds=d['runtime_seconds'],peak_RSS_GB=d['peak_RSS_GB']))
    p=root/'controls.json'
    if p.exists():
        for d in json.loads(p.read_text())['rows']:
            controls.append(dict(kind=d['kind'],L=d['L'],density={k:d['density'][k]
                for k in ('RMS','profile_peak','connected_peak')}))
    return dict(secondary_valley_channel=channels,pairing=pairing,periodic=rings,exact_controls=controls)


def main(root=HERE):
    obs={};rows=[];numerics=[]
    for p in sorted((root/'observations').glob('*.json')):
        d=json.loads(p.read_text());obs[p.name]=d;rows.append(summary(d,str(p.relative_to(root))))
    for p in sorted((root/'outputs').glob('*.json')):
        d=json.loads(p.read_text())
        if d['status']=='done':numerics.append(dict(tag=d['cfg']['tag'],**gates(d,root)))
    for r in rows:
        p=root/'outputs'/(r['tag']+'.json')
        r['run_status']=json.loads(p.read_text())['status'] if p.exists() else 'reused_completed_state'
    olda=json.loads((OLD/'analysis.json').read_text())
    for P in (0,1):
        r=next(r for r in olda['rows'] if r['tag']==f'legacy_V4_a1_L96_P{P}')
        numerics.append(dict(tag=f'reused_L96_P{P}_chi256',**{k:r[k] for k in
            ('gates','all_gates','dE','dS','dn','dS_chi','res_max','trunc','stopping_rule_met','schedule_limit_hit')}))
    cmp=[]
    for tag in sorted({r['tag'] for r in rows}):
        rs=sorted([r for r in rows if r['tag']==tag],key=lambda r:r['chi'])
        if len(rs)>=2:cmp.append(compare(rs[-2],rs[-1]))
        elif rs:
            d=obs[rs[0]['file'].split('/')[-1]]
            parent=d.get('cfg',{}).get('parent_observation')
            if parent:
                pd=obs[parent.split('/')[-1]];cmp.append(compare(summary(pd,parent),rs[0]))
    best={}
    for r in rows:
        key=(r['boundary'],r['L'],r['P'])
        if key not in best or (r['chi'],-r['E'])>(best[key]['chi'],-best[key]['E']):best[key]=r
    original=sorted([r for k,r in best.items() if k[0]=='full'],key=lambda r:(r['L'],r['P']))
    ref=min([r for r in original if r['L']==max(x['L'] for x in original)],key=lambda r:r['E'])
    fixedq={name:ref['operators'][name]['half']['profile_peak']['q'] for name in ('density','bond','kinetic')}
    fixed=[]
    for r in original:
        d=obs[r['file'].split('/')[-1]]
        fixed.append(dict(tag=r['tag'],L=r['L'],P=r['P'],chi=r['chi'],operators={k:common_q(d,q,k) for k,q in fixedq.items()}))
    start={}
    a=next((r for r in rows if r['tag']=='reused_L96_P0_chi256'),None)
    b=next((r for r in rows if r['tag']=='v4_L96_P0_s7' and r['chi']==256),None)
    if a and b:
        start=dict(dE=b['E']-a['E'],energy_gate=abs(b['E']-a['E'])<=1e-8,**compare(a,b))
        da=obs[a['file'].split('/')[-1]];db=obs[b['file'].split('/')[-1]]
        start['profile_max_differences']={k:float(np.max(abs(np.array(da['operators'][k]['profile'])-np.array(db['operators'][k]['profile'])))) for k in ('density','bond','kinetic')}
        p=root/'outputs'/'v4_L96_P0_s7.json'
        if p.exists():start['state_overlap']=json.loads(p.read_text()).get('independent_overlap')
    size_pairs=[]
    lengths=sorted({r['L'] for r in original})
    for L,L2 in zip(lengths[:-1],lengths[1:]):
        for P in (0,1):
            u=best.get(('full',L,P));v=best.get(('full',L2,P))
            if not u or not v:continue
            ratios={}
            for name in ('density','bond','kinetic'):
                ratios[name]={}
                for w in ('half','quarter'):
                    x=u['operators'][name][w];y=v['operators'][name][w]
                    ratios[name][w]=dict(RMS_ratio=y['RMS']/x['RMS'],profile_amplitude_ratio=y['profile_peak']['amplitude']/x['profile_peak']['amplitude'],
                        normalized_connected_peak_ratio=y['connected_peak']['connected_over_M']/x['connected_peak']['connected_over_M'])
                ratios[name]['long_distance_RMS_ratio']=v['operators'][name]['long_distance_RMS']/u['operators'][name]['long_distance_RMS']
            size_pairs.append(dict(lengths=[L,L2],P=P,tags=[u['tag'],v['tag']],ratios=ratios))
    eligibility=[]
    for P in (0,1):
        u=best.get(('full',96,P));v=best.get(('full',144,P))
        if u and v:
            x=u['operators']['density']['half'];y=v['operators']['density']['half']
            ar=y['profile_peak']['amplitude']/x['profile_peak']['amplitude'];cr=y['connected_peak']['connected_over_M']/x['connected_peak']['connected_over_M']
            eligibility.append(dict(P=P,amplitude_ratio=ar,connected_peak_weight_ratio=cr,condition_met=ar>.85 and cr<=.85))
    out=dict(rows=rows,has_interim_observations=any(r['run_status']=='running' for r in rows),
        numerical_gates=numerics,accuracy_comparisons=cmp,best_original_states=original,
        common_wavevector_reference=ref['tag'],common_wavevectors=fixedq,common_q_comparisons=fixed,
        independent_start=start,size_comparisons=size_pairs,boundary_probe_eligibility=eligibility,
        supplementary=supplementary(root,original,obs),
        long_distance_scope='Exploratory RMS over separations M/3 through M/2 in the central-half covariance; not an additional registered decision gate.',
        caution='Finite-size diagnostics only; no phase or central charge is assigned by this script.')
    write_json(root/'analysis.json',out)
    for r in original:
        d=r['operators']['density']['half'];print(r['tag'],r['chi'],'RMS',d['RMS'],'A',d['profile_peak']['amplitude'],
            'q/pi',d['profile_peak']['q_over_pi'],'Sconn/M',d['connected_peak']['connected_over_M'])
    for x in numerics:print('gates',x['tag'],x['gates'],x['dS_chi'])
    print('boundary',eligibility)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run-root',type=__import__('pathlib').Path,default=HERE);args=ap.parse_args()
    with threadpool_limits(limits=1):main(args.run_root)
