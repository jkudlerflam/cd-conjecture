#!/usr/bin/env python3
"""Summarize complete parity pairs; optionally freeze registered refinements."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import numpy as np
from scan import HERE, D, point


def load_points():
    points={}
    optional=[HERE/'targeted_refinement.json',HERE/'r3_extension.json']
    for p in [HERE/'grid.json']+sorted(HERE.glob('refinement_round*.json'))+[p for p in optional if p.exists()]+sorted(HERE.glob('r3_refinement_round*.json')):
        for x in json.loads(p.read_text())['points']:points[x['name']]=x
    return points


def summarize():
    points=load_points();groups={};records={};running={}
    paths=sorted((HERE/'validation_snapshots').glob('*.json'))+sorted((HERE/'dmrg').glob('*.json'))
    for p in paths:
        r=json.loads(p.read_text())
        if r['status']!='done':
            running.setdefault((r['name'],r['L']),[]).append(dict(parity=r['parity'],seed=r['seed'],source=str(p)))
            continue
        key=(r['name'],r['L'],r['parity'],r['seed'])
        if key not in records or r['chi']>=records[key]['chi']:
            r['_record_source']=str(p);records[key]=r
    for r in records.values():groups.setdefault((r['name'],r['L']),[]).append(r)
    rows=[]
    for (name,L),rr in sorted(groups.items()):
        required={0} if name=='control_xx' else {0,1}
        if {r['parity'] for r in rr}!=required:continue
        best=min(rr,key=lambda r:r['E_physical']);cfg=points[name]
        corr=best['measurements']['long_connected_correlations']
        n=np.array(best['measurements']['dens'])
        row=dict(name=name,L=L,parameters={k:cfg.get(k) for k in ['family','r','D','q','hz']},
            couplings=cfg['couplings'],metrics=cfg.get('metrics'),
            selected_parity=best['parity'],energy=best['E_physical'],chi=best['chi'],
            selected_source=best['_record_source'],running_upgrades=running.get((name,L),[]),
            density=best['density_bulk'],density_bulk_range=best['density_bulk_range'],
            density_global=float(np.mean(n)),
            central_density_staggered=float(np.mean(((-1)**np.arange(L))[L//4:3*L//4]*n[L//4:3*L//4])),
            S_fractional=best['S_fractional'],S_central3=best['S_central3'],profiles=best['profiles'],
            convergence=best['convergence'],both_sectors_converged=all(r['convergence']['sweep_converged'] for r in rr),
            long_correlation_rms={k:float(np.sqrt(np.mean(np.array(v)[len(v)//2:,1]**2))) for k,v in corr.items() if k!='origin'},
            sector_states=[dict(parity=r['parity'],seed=r['seed'],chi=r['chi'],source=r['_record_source'],E=r['E_physical'],n=r['density_bulk'],S=r['S_fractional'],convergence=r['convergence']) for r in rr])
        if len(rr)==2:row['parity_splitting']=next(r['E_physical'] for r in rr if r['parity']==1)-next(r['E_physical'] for r in rr if r['parity']==0)
        rows.append(row)
    pairs=[]
    for name in sorted({r['name'] for r in rows}):
        rr=sorted([r for r in rows if r['name']==name],key=lambda r:r['L'])
        for a,b in zip(rr,rr[1:]):
            f=6/math.log(b['L']/a['L'])
            pairs.append(dict(name=name,sizes=[a['L'],b['L']],
                c_fractional=f*(b['S_fractional']-a['S_fractional']),
                c_central3=f*(b['S_central3']-a['S_central3']),
                densities=[a['density'],b['density']],
                both_sizes_sectors_converged=a['both_sectors_converged'] and b['both_sectors_converged']))
    ladders={}
    for r in rows:
        if not r['name'].startswith('control'):
            ladders.setdefault((r['L'],r['parameters']['family'],r['parameters']['q']),[]).append(r)
    concavity=[]
    for key,rr in ladders.items():
        rr.sort(key=lambda r:Fraction(r['parameters']['hz']))
        for a,b in zip(rr,rr[1:]):
            ha,hb=[float(Fraction(r['parameters']['hz'])) for r in [a,b]]
            secant=(b['energy']-a['energy'])/(hb-ha)
            za=a['L']*(1-2*a['density_global']);zb=b['L']*(1-2*b['density_global'])
            violation=max(0.,secant-za,zb-secant)
            concavity.append(dict(names=[a['name'],b['name']],L=a['L'],field_secant=secant,
                expected_interval=[zb,za],violation=violation,pass_tolerance_1e5=violation<=1e-5))
    out=dict(rows=rows,pairs=pairs,field_concavity_checks=concavity)
    D.atomic_json(HERE/'summary.json',out)
    return out


def refine(out,round_number):
    points=load_points();groups={}
    for row in out['rows']:
        if row['L']==24 and not row['name'].startswith('control'):
            key=(row['parameters']['r'],row['parameters']['D'],row['parameters']['q'])
            groups.setdefault(key,[]).append(row)
    selected=[];decisions=[]
    for key,rows in sorted(groups.items()):
        rows=sorted(rows,key=lambda r:Fraction(r['parameters']['hz']))
        if min(abs(r['density']-.275) for r in rows)<=.015:
            decisions.append(dict(parameters=key,status='anchor_already_within_0.015'));continue
        brackets=[(a,b) for a,b in zip(rows,rows[1:]) if (a['density']-.275)*(b['density']-.275)<=0]
        if not brackets:
            decisions.append(dict(parameters=key,status='no_bracket'));continue
        a,b=min(brackets,key=lambda ab:Fraction(ab[1]['parameters']['hz'])-Fraction(ab[0]['parameters']['hz']))
        h=(Fraction(a['parameters']['hz'])+Fraction(b['parameters']['hz']))/2
        p=point(*key,h)
        if p['name'] in points:raise ValueError('Duplicate refinement')
        selected.append(p);decisions.append(dict(parameters=key,status='midpoint',field=str(h),
            bracket_names=[a['name'],b['name']],bracket_densities=[a['density'],b['density']]))
    target=HERE/f'refinement_round{round_number}.json'
    assert not target.exists(),'Refinement already frozen'
    D.atomic_json(target,dict(points=selected,decisions=decisions,round=round_number,
        registration='PROTOCOL.md density refinement',sizes=[24],chi=64))
    print('Frozen refinement',round_number,len(selected),'points',flush=True)


def promote(out):
    groups={};points=load_points()
    for r in out['rows']:
        if r['L']==24 and not r['name'].startswith('control'):groups.setdefault(r['parameters']['family'],[]).append(r)
    chosen=[];decisions=[]
    for family,rows in sorted(groups.items()):
        symmetric=[r for r in rows if r['parameters']['q']=='0']
        a=min(symmetric,key=lambda r:(abs(r['density']-.275),float(Fraction(r['parameters']['hz']))))
        offsets=[r for r in rows if r['parameters']['q']!='0']
        eligible=[r for r in offsets if abs(r['density']-.275)<=.08]
        def score(r):
            c1=r['profiles']['0.25']['c'];c2=r['profiles']['0.375']['c']
            return float(np.clip(min(c1,c2),0,3)-abs(c1-c2)/2-5*abs(r['density']-.275)-2*r['density_bulk_range'])
        b=max(eligible,key=score) if eligible else min(offsets,key=lambda r:abs(r['density']-.275))
        for r,reason in [(a,'symmetric closest density'),(b,'registered score' if eligible else 'unflagged closest-density control')]:
            chosen.append(points[r['name']]);decisions.append(dict(name=r['name'],reason=reason,density=r['density'],score=score(r),
                all_sectors_converged=r['both_sectors_converged']))
    target=HERE/'promoted.json';assert not target.exists(),'Promotion already frozen'
    D.atomic_json(target,dict(points=chosen,decisions=decisions,sizes=[48],chi=128,
        registration='PROTOCOL.md Stage B'))
    print(json.dumps(decisions,indent=2))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refine',type=int);ap.add_argument('--promote',action='store_true')
    args=ap.parse_args();out=summarize()
    if args.refine:refine(out,args.refine)
    elif args.promote:promote(out)
    else:
        for r in out['rows']:
            print(r['name'],r['L'],'P',r['selected_parity'],'n',round(r['density'],5),'S',round(r['S_fractional'],5),
                'cprof',round(r['profiles']['0.25']['c'],3),round(r['profiles']['0.375']['c'],3),'conv',r['both_sectors_converged'])
        for p in out['pairs']:print('PAIR',p)


if __name__=='__main__':main()
