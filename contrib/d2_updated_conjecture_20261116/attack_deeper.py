"""Further nested cancellation: a finite-ring diagnostic, not a c estimate.

H = D1 + eps F + eps^2(-T+2N-Pair) - 2 eps^4 N.
Extract the lowest band in the space excluding separations 1 through 4.
All operators in the microscopic Hamiltonian remain uniform nearest neighbor.
"""
import json
import time
from pathlib import Path
import numpy as np
from attack_nested import operators


def run(L, eps):
    E,_,F,W = operators(L)
    N = np.array([s.bit_count() for s in range(1 << L)])
    P4 = []
    for s in range(1 << L):
        ns = [(s >> j)&1 for j in range(L)]
        if not any(ns[j]*ns[(j+r)%L] for j in range(L) for r in range(1,5)):
            P4.append(s)
        for j in range(L):
            k = (j+1)%L
            if ns[j] == ns[k]: W[s^(1 << j)^(1 << k),s] -= 1
    H = np.diag(E-2*eps**4*N)+eps*F+eps**2*W
    ev,vec = np.linalg.eigh(H)
    m = len(P4)
    u,sv,vh = np.linalg.svd(vec[P4,:m])
    polar = u@vh
    He = (polar*ev[:m])@polar.T
    pos = {s:i for i,s in enumerate(P4)}
    numbers = N[P4]
    result = {'L':L,'eps':eps,'dim_P4':m,'band_gap_over_eps4':float((ev[m]-ev[m-1])/eps**4),
              'min_overlap':float(min(sv)),
              'onsite_over_eps6':float((He[pos[1],pos[1]]-He[pos[0],pos[0]])/eps**6),
              'hops_over_eps6':{str(r):float(He[pos[1 << r],pos[1]]/eps**6) for r in range(1,L//2+1)},
              'number_mixing_over_eps6':float(np.max(abs(He*(numbers[:,None] != numbers[None,:])))/eps**6)}
    print(json.dumps(result),flush=True)
    return result


if __name__ == '__main__':
    start = time.perf_counter()
    results = [run(L,e) for L in (10,11,12) for e in (0.075,0.05,0.035)]
    Path(__file__).with_name('deeper_results.json').write_text(json.dumps({
        'purpose':'Higher-order nested-constraint attack; no central charge assigned',
        'numpy_version':np.__version__,'results':results,
        'runtime_seconds':time.perf_counter()-start},indent=2)+'\n')
