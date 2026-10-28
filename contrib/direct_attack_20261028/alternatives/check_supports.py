#!/usr/bin/env python3
"""Exact finite certificate for the NN Pauli-to-Majorana support reduction.

Only the Python standard library is required. No floating point arithmetic,
phase diagram assumption, or graph-library recognition routine is used.
"""
from itertools import combinations, permutations
from pathlib import Path
import json

AXES = "XYZ"
LABELS = list(AXES) + [a+b for a in AXES for b in AXES]
INDEX = {s:i for i,s in enumerate(LABELS)}
TEMPLATES = [("X", "XY"), ("XY", "ZZ"), ("XY", "YZ"),
             ("X", "Y", "XX"), ("X", "XX", "YY"),
             ("XX", "YY", "ZZ")]
WITNESSES = [(("XY", 0), [("X", 1), ("XY", -1), ("XY", 1)]),
             (("XY", 0), [("XY", -1), ("XY", 1), ("ZZ", -1)]),
             (("XY", 0), [("XY", -1), ("XY", 1), ("YZ", -1)]),
             (("Y", 0), [("X", 0), ("XX", -1), ("XX", 0)]),
             (("YY", 0), [("X", 0), ("X", 1), ("XX", -1)]),
             (("XX", 0), [("YY", -1), ("YY", 1), ("ZZ", -1)])]

def support(lab, j, L=None):
    return {(j+k) % L if L else j+k:a for k,a in enumerate(lab)}

def anti(p, q):
    return sum(p[k] != q[k] for k in p.keys() & q.keys()) % 2

def mask(labels):
    return sum(1 << INDEX[a] for a in set(labels))

def names(s):
    return [label for i,label in enumerate(LABELS) if s >> i & 1]

for center, leaves in WITNESSES:
    c=support(*center); ls=[support(*x) for x in leaves]
    assert all(anti(c,x) for x in ls)
    assert all(not anti(x,y) for x,y in combinations(ls,2))

# Generate all axis permutations and reflections of the six templates.
forbidden=set()
for p in permutations(AXES):
    mp=dict(zip(AXES,p))
    for t in TEMPLATES:
        labs=[''.join(mp[a] for a in lab) for lab in t]
        forbidden.add(mask(labs))
        forbidden.add(mask([lab[::-1] for lab in labs]))
assert len(forbidden)==37

# An independent direct anticommutation computation on a seven-site ring.
# Every witness is local and has diameter at most three, so seven sites
# separate all possible overlaps without wraparound ambiguity.
L=7
ops=[support(lab,j,L) for lab in LABELS for j in range(L)]
adj=[set() for _ in ops]
for i,p in enumerate(ops):
    for j,q in enumerate(ops):
        if anti(p,q):adj[i].add(j)
direct=set()
for c in range(0,len(ops),L):
    for a,b,d in combinations(sorted(adj[c]),3):
        if b not in adj[a] and d not in adj[a] and d not in adj[b]:
            direct.add((1<<(c//L))|(1<<(a//L))|(1<<(b//L))|(1<<(d//L)))
minimal=[]
for f in sorted(direct,key=lambda s:(s.bit_count(),s)):
    if not any(f&g==g for g in minimal):minimal.append(f)
assert set(minimal)==forbidden
allowed=[s for s in range(1<<len(LABELS)) if not any(s&f==f for f in forbidden)]
maximal=[s for s in allowed if not any(s!=t and s&t==s for t in allowed)]

expected={mask(AXES)}
for a in AXES:
    other=[b for b in AXES if b!=a]
    expected.add(mask([a,a+a]))
    expected.add(mask(other+[a+a]))
    expected.add(mask([a+b for b in AXES]))
    expected.add(mask([b+a for b in AXES]))
    expected.add(mask([a]+[b+c for b in other for c in other]))
assert len(expected)==16
assert set(maximal)==expected

out={"status":"passed", "arithmetic":"exact integer and Boolean",
     "ring_length":L, "total_supports":1<<len(LABELS),
     "allowed_supports":len(allowed), "minimal_forbidden_count":len(forbidden),
     "maximal_allowed_count":len(maximal),
     "minimal_forbidden":[names(s) for s in sorted(forbidden)],
     "maximal_allowed":[names(s) for s in sorted(maximal)],
     "claw_witnesses":[{"support":t,"center":c,"leaves":ls}
                       for t,(c,ls) in zip(TEMPLATES,WITNESSES)]}
Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if not isinstance(v,list)},indent=2))
