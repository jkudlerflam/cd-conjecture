"""Exact comparison with a convenient uniform crossing bound."""
import numpy as np
from pathlib import Path
from positivity_C33 import bareiss
import json
records=[]
for n in range(1,24):
    filename=f'certificate_n{n}_D1000.npz'
    path=Path(filename) if Path(filename).exists() else Path('../work_C33')/filename
    with np.load(path,allow_pickle=True) as ar:
        L=ar['lower_numerators'];E=ar['radius_numerators'];d=int(ar['detC']);bits=int(ar['bits'])
    err=max(sum(row) for row in E);det=bareiss(L,positive=True)
    assert det is not None
    minors=0
    for j in range(n):
        idx=[i for i in range(n) if i!=j];minors+=bareiss(L[np.ix_(idx,idx)])
    num=det-err*minors
    assert num>0
    # S is the Hermitian part of -T, after removing alpha_n^4.
    # Its lower bound is 2*(det-err*minors)/(minors*2^(2bits)*detC^2).
    assert 2*num*10**250 > minors*(1<<(2*bits))*d*d
    records.append(dict(n=n,normalized_hermitian_part_exceeds_10_minus_250=True,
        exact_integer_comparison=True,source=str(path)))
json.dump(records,open('crossing_bounds.json','w'),indent=2)
print('Exact crossing bound verified for n=1,...,23')
