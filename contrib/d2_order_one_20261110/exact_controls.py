#!/usr/bin/env python3
"""Exact free-fermion OBC XX controls in the same physical cut coordinates."""
import json
import math
from pathlib import Path
import numpy as np
from dmrg_screen import entropy_stats


def xx_profile(L):
    h=np.diag(-np.ones(L-1),1)+np.diag(-np.ones(L-1),-1)
    e,v=np.linalg.eigh(h)
    C=v[:,e<0]@v[:,e<0].T
    out=[0.]
    for l in range(1,L):
        p=np.linalg.eigvalsh(C[:l,:l]);p=p[(p>1e-14)&(p<1-1e-14)]
        out.append(float(-np.sum(p*np.log(p)+(1-p)*np.log(1-p))))
    return np.array(out+[0.])


def main():
    rows=[]
    for legs in [1,2]:
        for L in [24,48,96,192]:
            if legs==1:S=xx_profile(L)[1:-1]
            else:
                single=xx_profile(L//2)
                S=np.array([single[(l+1)//2]+single[l//2] for l in range(1,L)])
            rows.append(dict(legs=legs,L=L,Sprof=S.tolist(),**entropy_stats(S,L)))
    pairs=[]
    for legs in [1,2]:
        r=[x for x in rows if x['legs']==legs]
        for a,b in zip(r,r[1:]):
            pairs.append(dict(legs=legs,sizes=[a['L'],b['L']],
                c_central3=6*(b['S_central3']-a['S_central3'])/math.log(2),
                c_fractional=6*(b['S_fractional']-a['S_fractional'])/math.log(2)))
    Path(__file__).with_suffix('.json').write_text(json.dumps(dict(rows=rows,pairs=pairs),indent=2)+'\n')
    print(json.dumps(pairs,indent=2))


if __name__=='__main__':main()
