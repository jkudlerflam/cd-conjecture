#!/usr/bin/env python3
"""Exact transport-count checks for the length-cutoff lemma.

Checks an abstract hole-transfer graph, a superset of physical marker
moves. Does not assert that every path has nonzero microscopic amplitude.
"""
from pathlib import Path
import json

def cost(r):
    return sum((j+3)//4 for j in range(1,r+1))

def closed_cost(r):
    q,t=divmod(r,4)
    return 2*q*(q+1)+t*(q+1)

def minimum_inbound(M,r):
    distances=sorted(min(j,M-j) for j in range(1,M) for _ in range(2))
    if r>len(distances):return None
    return sum(distances[:r])

checks=[]
for M in range(2,31):
    for r in range(1,2*(M-1)+1):
        actual=minimum_inbound(M,r)
        assert actual>=cost(r)==closed_cost(r)
        checks.append((M,r,actual,cost(r)))

# Explicit paths realize the bound in the unrestricted nearest-neighbor
# hole graph when enough initially occupied donor sites exist on both sides.
paths=[]
for r in range(1,25):
    radius=(r+3)//4
    sites=list(range(-radius,radius+1))
    holes={j:2 for j in sites}
    donors=[]
    for d in range(1,radius+1):
        donors.extend([-d,-d,d,d])
    moves=[]
    for source in donors[:r]:
        here=source
        while here:
            there=here+(1 if here<0 else -1)
            assert holes[here]>0
            holes[here]-=1;holes[there]+=1
            moves.append((here,there))
            here=there
    assert holes[0]==2+r
    assert len(moves)==cost(r)
    for here,there in reversed(moves):
        assert holes[there]>0
        holes[there]-=1;holes[here]+=1
    assert all(v==2 for v in holes.values())
    paths.append({"extra_holes":r,"inbound_moves":cost(r),
                  "round_trip_moves":2*cost(r),
                  "first_possible_lambda_degree":4*cost(r)})

out={"status":"passed", "finite_ring_checks":len(checks),
     "explicit_path_checks":len(paths),
     "scope":"abstract nearest-neighbor hole-transfer graph only",
     "thresholds":paths}
Path(__file__).with_name('length_cutoff_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='thresholds'},indent=2))
