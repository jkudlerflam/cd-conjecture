"""C36 finite checks. Python standard library, NumPy and SciPy only.
All generated files are local to this script's folder.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
from itertools import combinations, product
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal

HERE = Path(__file__).resolve().parent

def hexagons(p):
    reps = [(1,a,b) for a in range(p) for b in range(p)]
    reps += [(0,1,a) for a in range(p)] + [(0,0,1)]
    vec = np.asarray(reps, dtype=int)
    incidence = (vec @ vec.T) % p == 0
    edges = list(zip(*np.nonzero(incidence)))
    edge_id = {edge:i+1 for i,edge in enumerate(edges)}
    size = len(vec); n = len(edges)
    common = {}
    for a,b in combinations(range(size),2):
        ids = np.flatnonzero(incidence[a] & incidence[b])
        assert len(ids) == 1
        common[a,b] = int(ids[0])
    # Count unordered 3+3 partitions once, by requiring position 0 in the first triple.
    choices = [(0,)+pair for pair in combinations(range(1,6),2)]
    details = []; count = 0; partitions = 0; histogram = {}
    for a,b,c in combinations(range(size),3):
        ab,ac,bc = common[a,b],common[a,c],common[b,c]
        if len({ab,ac,bc}) != 3:
            continue
        count += 1
        pairs = [(a,ab),(b,ab),(b,bc),(c,bc),(c,ac),(a,ac)]
        labels = [edge_id[e] for e in pairs]
        # Exact species neutrality, without a floating-point nullspace.
        balance = np.zeros(2*size,dtype=int)
        for i,(point,line) in enumerate(pairs):
            sign = (-1)**i
            balance[point] += sign; balance[size+line] -= sign
        assert np.all(balance == 0)
        equal = []
        for pos in choices:
            other = [i for i in range(6) if i not in pos]
            left = [labels[i] for i in pos]
            right = [labels[i] for i in other]
            if sum(left) == sum(right):
                equal.append([left,right])
        histogram[len(equal)] = histogram.get(len(equal),0)+1
        partitions += len(equal)
        if equal:
            details.append(dict(labels=labels,equal_partitions=equal))
        # Independently compare the exact momentum congruence, all 64 chiral signs.
        if p <= 3:
            direct = 0
            for signs in product([-1,1],repeat=6):
                numerator = sum(sign*(3*(n+1)+label) for sign,label in zip(signs,labels))
                assert numerator % (24*(n+1)) != 12*(n+1)
                if numerator % (24*(n+1)) == 0:
                    direct += 1
                    assert sum(signs) == 0
                    assert sum(sign*label for sign,label in zip(signs,labels)) == 0
            assert direct == 2*len(equal)
    assert count == p**3*(p+1)*(p*p+p+1)//6
    summary = dict(p=p,q=2*size,N=n,hexagons=count,
                   resonant_hexagons=len(details),unordered_equal_partitions=partitions,
                   chiral_assignments_one_z_orientation=2*partitions,
                   partitions_per_hexagon=histogram)
    if p <= 3:
        (HERE/f'hexagons_p{p}.json').write_text(json.dumps(details,indent=2)+'\n')
    return summary

def all_band_kernel():
    # Legal domain-wall basis |r>=|1^r 0^(L-r)>, r=0,...,L.
    # A local range-three spin Hamiltonian restricts to this Jacobi matrix.
    # Endpoint components are nonzero analytically by the eigenvector recurrence.
    out=[]
    delta=1.; s=.1; t=.4
    for length in [4,8,16,32]:
        vals,vecs=eigh_tridiagonal(np.arange(length+1,dtype=float),
                                 np.full(length,t),lapack_driver='stev')
        spread=float(vals[-1]-vals[0])
        pref=abs(vecs[0,-1])*abs(vecs[0,0])**2*s*s
        assert spread >= length-1e-10
        assert pref > 0
        # Asymptotic norm lower bound from one nonzero matrix entry.
        out.append(dict(L=length,matrix_dimension=length+1,
                        spectral_spread=spread,delta=delta,s=s,
                        exponential_rate=s*spread-delta,
                        log_prefactor=float(np.log(pref)),
                        ground_kernel_decay_rate=delta))
    return out

def active_excursion_majorant():
    # An explicit finite positive majorant generator with branching proportional
    # to defect number. Its column sums obey the proved absorption inequality.
    rng=np.random.default_rng(36)
    out=[]
    for states in [8,32,128]:
        weights=rng.integers(1,8,size=states)
        matrix=rng.random((states,states))
        matrix *= weights / matrix.sum(axis=0)
        delta=1.; sb=.1
        generator=sb*matrix-np.diag(delta*weights)
        bound=-(delta-sb)
        assert np.max(generator.sum(axis=0)) <= bound+1e-12
        spectral=float(np.max(np.real(np.linalg.eigvals(generator))))
        assert spectral <= bound+1e-10
        out.append(dict(states=states,max_column_sum=float(np.max(generator.sum(axis=0))),
                        spectral_abscissa=spectral,upper_bound=bound))
    return out

def main():
    result=dict(hexagons=[hexagons(p) for p in [2,3,5,7]],
                all_band_kernel=all_band_kernel(),
                active_majorants=active_excursion_majorant())
    (HERE/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
