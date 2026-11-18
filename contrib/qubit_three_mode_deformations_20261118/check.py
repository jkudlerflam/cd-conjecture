"""Reproduce the normal map and the all-length local saturation certificate."""
from __future__ import annotations
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import block_certificate as C
from normal_map import A, eliminate, expected, leading

HERE=Path(__file__).resolve().parent


def embedded_block_sum(L, tau, mu):
    """Independently embed every 12-site compressed block in a full ring."""
    states = leading.basis(L, 3)
    lookup = {int(s):i for i, s in enumerate(states)}
    blocks = {}
    for n in range(5):
        local, mat = C.block(12, n, tau, mu)
        blocks[n] = (local, {s:i for i, s in enumerate(local)}, mat)
    rows, cols, data = [], [], []
    for col, ss in enumerate(states):
        state = int(ss)
        for start in range(L):
            sites = [(start+j) % L for j in range(12)]
            local = sum(((state >> j) & 1) << k for k, j in enumerate(sites))
            local_states, local_lookup, matrix = blocks[local.bit_count()]
            column = matrix[:, local_lookup[local]]
            outside = state
            for j in sites:
                outside &= ~(1 << j)
            for row in np.flatnonzero(abs(column) > 1e-15):
                dest = outside | sum(((local_states[row] >> k) & 1) << j for k, j in enumerate(sites))
                if dest in lookup:
                    rows.append(lookup[dest]); cols.append(col); data.append(column[row])
    return states, sparse.coo_matrix((data,(rows,cols)),shape=(len(states),len(states))).tocsr()


def rational_psd_certificates():
    reports = {}
    # Negative imaginary hopping is complex conjugation of the positive case.
    for real, imag in (("1/10", "0"), ("-1/10", "0"), ("0", "1/10")):
        label=f'{real}+i({imag})'
        reports[label]=[]
        for n in range(5):
            matrix = C.exact_block(12, n, real, imag, "7/100")
            report = C.exact_ldl(matrix)
            report['N'] = n
            if n==4:
                words=C.words(12,n)
                crystals=[sum(1<<j for j in range(r,12,3)) for r in range(3)]
                indices=[words.index(s) for s in crystals]
                if imag!='0':indices += [i+len(words) for i in indices]
                assert all(matrix[row][col]==0 for row in range(len(matrix)) for col in indices)
                assert report['zero_pivots']==len(indices)
                report['crystal_kernel_verified']=True
            # Hash the exact input; fractions have a canonical textual form.
            report['rational_matrix_sha256'] = hashlib.sha256(
                '\n'.join(','.join(map(str,row)) for row in matrix).encode()).hexdigest()
            reports[label].append(report)
    return reports


def main():
    start=time.monotonic()
    reports={}
    with threadpool_limits(1):
        cases=[dict(a=0j,zeta=(1+2j)/100,eta=(2-1j)/200,mu=1/7,
                    beta=-2+1j,pair5=1+3j,pair6=-2j),
               dict(a=(2+1j)/8,zeta=(-2+1j)/100,eta=(-1+2j)/200,mu=-1/9,
                    beta=1-2j,pair5=-2+1j,pair6=3+1j)]
        comparisons=[]
        for L in (10,12,13):
            for index, pars in enumerate(cases):
                states, matrix, diag = eliminate(L, **pars)
                other_states, target = expected(L, **pars)
                assert np.array_equal(states, other_states)
                error=A.maxabs(matrix-target)
                assert error < 1e-10, (L,index,error)
                assert diag['third_order_norm'] < 1e-10
                assert diag['fifth_order_norm'] < 1e-10
                comparisons.append(dict(L=L,case=index,dimension=len(states),full_matrix_error=error,**diag))
        reports['full_normal_map']=comparisons
        reports['rational_psd']=rational_psd_certificates()
        embedding=[]
        for L in (15,18):
            states, actual=embedded_block_sum(L,(1+1j)/200,7/1000)
            _, target=expected(L,states=states,zeta=(1+1j)/200,mu=7/1000)
            target=target+sparse.eye(len(states))*float(Fraction(7*L,3000))
            error=A.maxabs(actual-target)
            assert error < 1e-12
            embedding.append(dict(L=L,dimension=len(states),block_sum_error=error))
        reports['block_embedding']=embedding
    reports['source_sha256']={str(p.relative_to(HERE.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (HERE.parent/'qubit_cancellation_classification_20261116/algebra.py',
                  HERE.parent/'qubit_resonant_liquids_20261116/leading.py')}
    reports['runtime_seconds']=time.monotonic()-start
    (HERE/'checks.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps({key:val for key,val in reports.items() if key!='rational_psd'},indent=2))
    print('Exact rational LDL certificates: 15 blocks passed; stored in checks.json.')


if __name__=='__main__':
    main()
