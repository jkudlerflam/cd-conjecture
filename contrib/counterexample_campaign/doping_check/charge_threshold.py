#!/usr/bin/env python3
"""Momentum-resolved free-charge threshold between internal twist sectors.

Dynamic programming minimizes over every M-particle occupation of L orbitals.
No restriction to a single particle-hole pair or low-energy excitations.
Source sector is theta=0 and physical translation momentum zero.
"""
from itertools import combinations
from pathlib import Path
import json
import math
import numpy as np


def spectrum_minimum(L,M,r):
    theta=2*np.pi*r/M
    k=(2*np.pi*np.arange(L)-theta-np.pi*(M-1))/L
    eps=-2*np.cos(k)
    dp=np.full((M+1,L),np.inf)
    dp[0,0]=0
    for n,e in enumerate(eps):
        for m in range(min(M,n+1),0,-1):
            dp[m]=np.minimum(dp[m],np.roll(dp[m-1],n)+e)
    # Physical translation phase is exp[2pi i (r+M(M-1)/2-sum n)/L].
    required=(r+M*(M-1)//2)%L
    return float(dp[M,required]), float(np.sort(eps)[:M].sum()), eps, dp[M]


def run(L,frac):
    M=L//2
    assert L%2==0 and (M*frac[0])%frac[1]==0
    r=M*frac[0]//frac[1]
    theta=2*np.pi*r/M
    energy,own_sea,eps,dp=spectrum_minimum(L,M,r)
    reference=-2/np.sin(np.pi/L)
    threshold=energy-reference
    distance=abs((theta+np.pi)%(2*np.pi)-np.pi)
    bound=(2/np.pi)*distance-2*np.tan(np.pi/(2*L))
    assert threshold>=bound-1e-9,(L,r,threshold,bound)
    if L<=16:
        brute=np.full(L,np.inf)
        for c in combinations(range(L),M):
            residue=sum(c)%L
            brute[residue]=min(brute[residue],sum(eps[i] for i in c))
        assert np.max(np.abs(brute-dp))<1e-11
    return {"L":L,"M":M,"r":r,"theta_over_pi":float(theta/np.pi),"fixed_momentum_threshold":float(threshold),"proved_lower_bound":float(bound),"own_sea_energy_difference":float(own_sea-reference),"candidate_thermodynamic_edge":float(2*np.sin(distance/2)),"all_occupations_optimized":math.comb(L,M)}


if __name__=='__main__':
    rows=[run(L,frac) for L in (16,32,64,128,256) for frac in ((1,8),(1,4),(3,8),(1,2))]
    Path(__file__).with_name('charge_threshold_results.json').write_text(json.dumps({"status":"passed","t":1,"density":0.5,"rows":rows},indent=2)+'\n')
    for r in rows:print(r['L'],r['theta_over_pi'],r['fixed_momentum_threshold'],r['proved_lower_bound'])
