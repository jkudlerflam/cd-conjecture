#!/usr/bin/env python3
"""All-occupation threshold checks with bounded rounding away from half density."""
from pathlib import Path
import json
import math
import numpy as np
from charge_threshold import spectrum_minimum


def run(L,M,frac):
    assert (M*frac[0])%frac[1]==0
    r=M*frac[0]//frac[1]
    theta=2*np.pi*r/M
    energy,own_sea,eps,dp=spectrum_minimum(L,M,r)
    reference=-2*np.sin(np.pi*M/L)/np.sin(np.pi/L)
    threshold=energy-reference
    distance=abs((theta+np.pi)%(2*np.pi)-np.pi)
    a=abs(M-L/2)
    A=math.ceil(a)+1
    bound=(2/np.pi)*distance-(8*a+(8+4*np.pi)*A*A)/L-2*np.tan(np.pi/(2*L))
    assert threshold>=bound-1e-8
    # Sharper directly evaluated finite-ring inequality from mu=0 reference.
    sea=set(np.argsort(eps)[:M]);zero_sea={n for n,e in enumerate(eps) if e<0}
    mismatch=sea.symmetric_difference(zero_sea)
    k=(2*np.pi*np.arange(L)-theta-np.pi*(M-1))/L
    distances=np.arcsin(np.minimum(1,np.abs(np.cos(k))))
    D=float(sum(distances[n] for n in mismatch))
    W=float(sum(abs(np.cos(k[n])) for n in mismatch))
    K=(r+M*(M-1)//2-sum(sea))%L
    P=2*np.pi*K/L
    pdist=abs((P+np.pi/2)%np.pi-np.pi/2)
    sharp=(4/np.pi)*(pdist-D)-2*W+own_sea-reference
    assert threshold>=sharp-1e-8
    return {"L":L,"M":M,"theta_over_pi":float(theta/np.pi),"density":M/L,"rounding_abs":float(a),"fixed_momentum_threshold":float(threshold),"proved_uniform_rounding_bound":float(bound),"proved_evaluated_bound":float(sharp),"reference_mismatch_count":len(mismatch),"reference_deviation_sum":D}


if __name__=='__main__':
    rows=[run(L,M,f) for L,M in [(33,16),(65,32),(129,64),(257,128),(193,96),(197,96),(205,96)] for f in [(1,4),(1,2)]]
    Path(__file__).with_name('rounded_threshold_results.json').write_text(json.dumps({"status":"passed","t":1,"rows":rows},indent=2)+'\n')
    for r in rows:print(r['L'],r['M'],r['theta_over_pi'],r['fixed_momentum_threshold'],r['proved_uniform_rounding_bound'])
