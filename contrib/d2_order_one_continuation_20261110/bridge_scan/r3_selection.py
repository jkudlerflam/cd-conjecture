#!/usr/bin/env python3
"""Freeze only the r3 extension's registered refinement or promotion."""
import argparse
from fractions import Fraction
from scan import HERE, D, point
from summarize import summarize, load_points

ap=argparse.ArgumentParser();ap.add_argument('--refine',type=int);ap.add_argument('--promote',action='store_true')
args=ap.parse_args();out=summarize();points=load_points()
rows=[r for r in out['rows'] if r['L']==24 and r['parameters']['family']=='r3_d1_8']
assert len(rows)>=10,'Complete the initial r3 grid first'
selected=[];decisions=[]
if args.refine:
    assert args.refine in (1,2)
    for q in ['0','1/16']:
        rr=sorted([r for r in rows if r['parameters']['q']==q],key=lambda r:Fraction(r['parameters']['hz']))
        if min(abs(r['density']-.275) for r in rr)<=.015:
            decisions.append(dict(q=q,status='anchor_within_0.015'));continue
        brackets=[(a,b) for a,b in zip(rr,rr[1:]) if (a['density']-.275)*(b['density']-.275)<=0]
        if not brackets:
            decisions.append(dict(q=q,status='no_bracket'));continue
        a,b=min(brackets,key=lambda ab:Fraction(ab[1]['parameters']['hz'])-Fraction(ab[0]['parameters']['hz']))
        h=(Fraction(a['parameters']['hz'])+Fraction(b['parameters']['hz']))/2
        x=point('3','1/8',q,h);assert x['name'] not in points
        selected.append(x);decisions.append(dict(q=q,status='midpoint',h=str(h),bracket=[a['name'],b['name']],densities=[a['density'],b['density']]))
    path=HERE/f'r3_refinement_round{args.refine}.json'
elif args.promote:
    for q in ['0','1/16']:
        rr=[r for r in rows if r['parameters']['q']==q]
        eligible=[r for r in rr if .2<=r['density']<=.35]
        def rank(r):
            a,b=[r['profiles'][k]['c'] for k in ['0.25','0.375']]
            return min(a,b),-abs(a-b),-abs(r['density']-.275)
        x=max(eligible,key=rank) if eligible else min(rr,key=lambda r:abs(r['density']-.275))
        selected.append(points[x['name']]);decisions.append(dict(name=x['name'],reason='largest minimum profile slope in density interval' if eligible else 'unflagged closest-density control'))
    a=min([r for r in rows if r['parameters']['q']=='0'],key=lambda r:abs(r['density']-.275))
    if a['name'] not in {x['name'] for x in selected}:
        selected.append(points[a['name']]);decisions.append(dict(name=a['name'],reason='symmetric closest-density check'))
    path=HERE/'r3_promoted.json'
else:raise SystemExit('Use --refine or --promote')
assert not path.exists(),'Selection already frozen'
D.atomic_json(path,dict(points=selected,decisions=decisions,registration='R3_EXTENSION.md'))
print(path.name,decisions,flush=True)
