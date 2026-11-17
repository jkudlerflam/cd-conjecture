#!/usr/bin/env python3
"""Post-hoc numerical entropy slopes; no automatic asymptotic c assignment."""
import argparse,json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from phase_common import HERE,OLD,write_json


def fit(x,y,factor):
    x=np.array(x);y=np.array(y);X=np.column_stack([np.ones(len(x)),x])
    beta=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@beta
    return dict(c_ent=float(factor*beta[1]),intercept=float(beta[0]),
        residual_RMS=float(np.sqrt(np.mean(res**2))),residuals=res.tolist())


def summarize(L,S,factor=6):
    S=np.asarray(S);cuts=np.arange(1,L);x=np.log(L/np.pi*np.sin(np.pi*cuts/L))
    assert len(S)==L-1 and np.all(np.isfinite(S));averages={}
    masks={f'W{w}':abs(cuts-L//2)<=w for w in (0,1,6,12)}
    masks['central_quarter']=(cuts>=3*L/8)&(cuts<=5*L/8)
    for k,m in masks.items():
        averages[k]=dict(S=float(S[m].mean()),mean_log_chord=float(x[m].mean()),cuts=cuts[m].tolist())
    profiles=[]
    for edge in (1/8,1/4,1/3):
        for halfwidth in (0,3):
            kernel=np.ones(2*halfwidth+1)/(2*halfwidth+1)
            sx=np.convolve(x,kernel,'valid');sy=np.convolve(S,kernel,'valid')
            sc=cuts[halfwidth:len(cuts)-halfwidth] if halfwidth else cuts
            m=(sc>=edge*L)&(sc<=(1-edge)*L)
            if np.count_nonzero(m)<3 or np.ptp(sx[m])<1e-12:continue
            profiles.append(dict(edge_fraction=edge,smoothing_halfwidth=halfwidth,
                cuts=sc[m].tolist(),**fit(sx[m],sy[m],factor)))
    return dict(averages=averages,profile_fits=profiles)


def size_fits(rows,factor=6):
    rows=sorted(rows,key=lambda r:r['L']);pairs=[];global_fits=[]
    for a,b in zip(rows[:-1],rows[1:]):
        estimates={}
        for window in a['averages']:
            x=a['averages'][window];y=b['averages'][window];ds=y['S']-x['S']
            estimates[window]=dict(c_ent_logL=float(factor*ds/np.log(b['L']/a['L'])),
                c_ent_log_chord=float(factor*ds/(y['mean_log_chord']-x['mean_log_chord'])),
                entropy=[x['S'],y['S']])
        pairs.append(dict(lengths=[a['L'],b['L']],parities=[a.get('P'),b.get('P')],
            tags=[a['tag'],b['tag']],estimates=estimates))
    if len(rows)>=3:
        for window in rows[0]['averages']:
            y=[r['averages'][window]['S'] for r in rows]
            global_fits.append(dict(window=window,logL=fit(np.log([r['L'] for r in rows]),y,factor),
                log_chord=fit([r['averages'][window]['mean_log_chord'] for r in rows],y,factor)))
    return dict(states=[r['tag'] for r in rows],pairs=pairs,global_fits=global_fits)


def main(root=HERE):
    synthetic=[dict(tag=f'synthetic_{L}',L=L,**summarize(L,
        .5*np.log(L/np.pi*np.sin(np.pi*np.arange(1,L)/L))+.7)) for L in (48,96,144)]
    for p in size_fits(synthetic)['pairs']:
        assert max(abs(e['c_ent_log_chord']-3) for e in p['estimates'].values())<1e-12
    for r in synthetic:
        assert max(abs(p['c_ent']-3) for p in r['profile_fits'])<1e-12
    analysis=json.loads((root/'analysis.json').read_text());allrows=[];selected=[]
    chosen={r['file'] for r in analysis['best_original_states']}
    for p in sorted((root/'observations').glob('*.json')):
        d=json.loads(p.read_text());S=d.get('S')
        if S is None:S=json.loads((OLD/'outputs'/f"legacy_V4_a1_L96_P{d['P']}.json").read_text())['S']
        row=dict(tag=p.stem,L=d['L'],P=d['P'],chi=d['chi'],E=d['E'],**summarize(d['L'],S))
        allrows.append(row)
        if str(p.relative_to(root)) in chosen:selected.append(row)
    ground=[min([r for r in selected if r['L']==L],key=lambda r:r['E']) for L in sorted({r['L'] for r in selected})]
    sectors={str(P):size_fits([r for r in selected if r['P']==P]) for P in (0,1)}
    controls=[];controlfits={};p=root/'controls.json'
    if p.exists():
        for d in json.loads(p.read_text())['rows']:
            if 'S' not in d:continue
            controls.append(dict(tag=d['kind']+f"_L{d['L']}",kind=d['kind'],L=d['L'],**summarize(d['L'],d['S'])))
        controlfits={kind:size_fits([r for r in controls if r['kind']==kind]) for kind in sorted({r['kind'] for r in controls})}
    rings=[]
    for p in sorted((root/'periodic').glob('*.json')):
        d=json.loads(p.read_text())
        if 'S' in d:rings.append(dict(tag=p.stem,L=d['L'],P=d['P'],E=d['energies'][0],**summarize(d['L'],d['S'],3)))
    rground=[min([r for r in rings if r['L']==L],key=lambda r:r['E']) for L in sorted({r['L'] for r in rings})]
    # Explicit changes in the primary W6 slope caused only by increasing the
    # largest-size chi. This is a numerical sensitivity, not a confidence bound.
    chi=[]
    if selected:
        largest=max(r['L'] for r in selected);lower=max(r['L'] for r in selected if r['L']<largest)
        for P in (0,1):
            rs=sorted([r for r in allrows if r['L']==largest and r['P']==P],key=lambda r:r['chi'])
            if len(rs)>=2:
                a,b=rs[-2:];ds=b['averages']['W6']['S']-a['averages']['W6']['S']
                chi.append(dict(P=P,chis=[a['chi'],b['chi']],delta_Sbar6=ds,
                    delta_c_ent_logL=float(6*ds/np.log(largest/lower))))
    out=dict(analysis='Post-hoc at the user request; natural logarithms; c_ent denotes a finite-size entropy slope, not c_eff.',
        formula='S(l,L)=(c/6)log[(L/pi)sin(pi*l/L)]+constant for open chains; periodic coefficient c/3.',
        source='https://arxiv.org/abs/hep-th/0405152',asymptotic_c=None,
        has_interim_observations=analysis['has_interim_observations'],all_candidate_states=allrows,
        selected_states=selected,ground_branch=size_fits(ground),fixed_parity=sectors,
        largest_chi_sensitivity=chi,exact_controls=controls,control_size_fits=controlfits,
        periodic_states=rings,periodic_ground_branch=size_fits(rground,3) if rings else None,
        limitations=['Only three large-chain sizes; no asymptotic extrapolation is identified.',
            'Numerical gates and independent-start checks remain in analysis.json.',
            'Profile and size fits are sensitivity diagnostics; no statistically independent fit errors are claimed.',
            'Finite entanglement slopes do not establish positive-velocity conformal modes.'])
    write_json(root/'entropy_analysis.json',out)
    for p in out['ground_branch']['pairs']:print('ground',p['lengths'],p['parities'],p['estimates'])
    for P,x in sectors.items():
        for p in x['pairs']:print('parity',P,p['lengths'],p['estimates']['W6'])
    for x in out['ground_branch']['global_fits']:print('three-size',x['window'],x['logL']['c_ent'],x['log_chord']['c_ent'])
    print('largest chi sensitivity',chi)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run-root',type=Path,default=HERE);args=ap.parse_args()
    with threadpool_limits(limits=1):main(args.run_root)
