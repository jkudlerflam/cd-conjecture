#!/usr/bin/env python3
"""Independently enumerate all neutral quadratic and quartic momenta.

The supplied rational root intervals are input, not independently proved
here. No campaign checking implementation is imported. Exact rational
comparisons certify all distinct pair products by their circular neighbors;
an independent vectorized all-pairs pass counts every transfer.
"""
from fractions import Fraction as F
from functools import cmp_to_key
from itertools import combinations_with_replacement
from pathlib import Path
import json
import math
import time
import numpy as np


def order(a,b):
    x,y=a['z'];u,v=b['z']
    half=lambda p,q:int(q<0 or (q==0 and p<0))
    ha,hb=half(x,y),half(u,v)
    if ha!=hb:return ha-hb
    cross=x*v-y*u
    return -int(cross>0)+int(cross<0)


def exact_minimum(items):
    seq=sorted(items,key=cmp_to_key(order));best=None
    for a,b in zip(seq,seq[1:]+seq[:1]):
        x,y=a['z'];u,v=b['z']
        assert x*v-y*u>=0
        lower=max(abs(x-u),abs(y-v))-a['err']-b['err']
        assert lower>0
        if best is None or lower<best[0]:best=(lower,[a['id'],b['id']])
    return best


def main():
    start=time.monotonic()
    root=Path(__file__).resolve().parents[3]
    data=json.loads((root/'codex/work_C33/root_certificate_exact_n23.json').read_text())
    roots=[]
    for i,ends in enumerate(data['intervals']):
        lo,hi=(F(int(v[0]),int(v[1])) for v in ends)
        mid=(lo+hi)/2
        roots.append(dict(z=((1-mid**2)/(1+mid**2),2*mid/(1+mid**2)),err=hi-lo,id=[i]))
    pair=[]
    for i,j in combinations_with_replacement(range(len(roots)),2):
        x,y=roots[i]['z'];u,v=roots[j]['z']
        pair.append(dict(z=(x*u-y*v,x*v+y*u),err=roots[i]['err']+roots[j]['err'],id=[i,j]))
    d2,w2=exact_minimum(roots);d4,w4=exact_minimum(pair)
    # Explicitly enumerate every ordered pair transfer, distinct or not.
    phase=np.array([math.atan2(float(r['z'][1]),float(r['z'][0])) for r in pair])
    dif=(phase[:,None]-phase[None,:]+math.pi)%(2*math.pi)-math.pi
    mask=~np.eye(len(pair),dtype=bool)
    assert np.min(abs(dif[mask]))>float(d4)
    nonzero=int(np.count_nonzero(abs(dif)>1e-12))
    assert nonzero==len(pair)*(len(pair)-1)
    # At Delta_payload=2, dangerous charge branches obey |Q|<2*pi*b.
    # For each transfer q use its reduced representative; minimum distance
    # from a reciprocal physical vector over every allowed ell is computed.
    q=dif[mask];minimum=math.inf;branches=0;b=24
    for ell in range(-25,26):
        Q=q+2*math.pi*ell
        use=abs(Q)<2*math.pi*b
        if not np.any(use):continue
        k=Q[use]/(2*b)
        distance=abs((k+math.pi)%(2*math.pi)-math.pi)
        minimum=min(minimum,float(np.min(distance)))
        branches+=len(k)
    assert minimum>float(d4)/(2*b)
    answer=dict(status='passed',root_count=len(roots),unordered_pairs_with_replacement=len(pair),quadratic_ordered_nonforward_count=len(roots)*(len(roots)-1),quartic_ordered_nonforward_count=nonzero,quartic_dangerous_harmonic_branches=branches,delta2_exact=[str(d2.numerator),str(d2.denominator)],delta4_exact=[str(d4.numerator),str(d4.denominator)],delta2=float(d2),delta4=float(d4),witness2=w2,witness4=w4,physical_margin_certified=float(min(d2,d4))/(2*b),physical_margin_all_quartic_branches=minimum,elapsed_seconds=time.monotonic()-start,input_scope='Rational Fermi intervals inherited; their matrix-inertia proof is not rerun here.')
    Path(__file__).with_name('vertices_results.json').write_text(json.dumps(answer,indent=2)+'\n')
    print(json.dumps({k:v for k,v in answer.items() if not k.endswith('_exact')},indent=2))


if __name__=='__main__':main()
