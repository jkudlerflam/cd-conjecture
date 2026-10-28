"""Explicit number-conserving target with polynomial momentum margins.

All resonance checks use integers. No CFT measurement is performed.
"""
import json
from itertools import combinations_with_replacement
import numpy as np


def scaffold(m):
    assert m >= 1
    N = 2*m
    base = 2*N+1
    a = [i+base*i*i for i in range(1, N+1)]
    A = a[-1]
    # p_i/pi = numerator_i / denominator, in (1/4,3/8].
    den = 8*A
    nums = [2*A+x for x in a]
    return dict(m=m,base=base,A=A,denominator=den,numerators=nums,
                nu_numerator=1,nu_denominator=2)


def target(m):
    s = scaffold(m)
    p = np.pi*np.array(s['numerators'],dtype=float)/s['denominator']
    left,right = p[:m],p[m:]
    widths = (right-left)/2
    return dict(left=left,right=right,width=widths,
                theta=-(right+left)/2,mu=2*np.cos(widths),
                velocity=2*np.sin(widths))


def checks(m):
    s = scaffold(m)
    nums,den = s['numerators'],s['denominator']
    # Equal sums are allowed only for equal unordered pairs. This includes
    # repeated branch indices, a stronger test than the fermionic quartic list.
    sums = {}
    for i,j in combinations_with_replacement(range(2*m),2):
        v = nums[i]+nums[j]
        assert v not in sums, (m,sums.get(v),(i,j))
        sums[v]=(i,j)
    ordered = sorted(sums)
    gap = min(b-a for a,b in zip(ordered,ordered[1:]))
    assert gap >= 1
    # Pair-sum differences are smaller than den/4, so no 2pi wrap exists.
    assert ordered[-1]-ordered[0] < den/4
    min_width_num = min(nums[m+i]-nums[i] for i in range(m))
    # Half width > pi/64, upper bound pi/16.
    assert 32*min_width_num > den
    assert max(nums[m+i]-nums[i] for i in range(m))*8 < den
    # Mixed h=+/-1 plus a number-conserving bilinear has margin >7pi/8.
    maxdiff = nums[-1]-nums[0]
    assert 8*maxdiff < den
    t = target(m)
    return dict(m=m,branch_count=2*m,pair_sum_count=len(sums),
                denominator=den,minimum_pair_sum_spacing_numerator=gap,
                certified_payload_margin_over_pi=f'1/{den}',
                minimum_velocity=float(min(t['velocity'])),
                maximum_velocity=float(max(t['velocity'])),
                velocity_lower_bound=float(2*np.sin(np.pi/64)),
                mixed_margin_lower_bound_over_pi='7/8')


if __name__ == '__main__':
    out = [checks(m) for m in (1,2,3,8,18,32,64,128,256)]
    with open('rational_target_checks.json','w') as f:
        json.dump(out,f,indent=2)
        f.write('\n')
    print(json.dumps(out,indent=2))
