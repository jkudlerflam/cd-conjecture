"""Small ground-kernel checks, no spectral classification inferred."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np

def matrix(h,L,periodic):
    H=np.zeros((2**L,2**L),complex)
    for s in range(2**L):
        for i in range(L if periodic else L-1):
            j=(i+1)%L; x=2*((s >> i)&1)+((s >> j)&1)
            for y in range(4):
                target=(s & ~(1 << i) & ~(1 << j)) | ((y >> 1) << i) | ((y & 1) << j)
                H[target,s]+=h[y,x]
    return H

def projector(v):
    v=np.array(v,complex); v/=np.linalg.norm(v)
    return np.outer(v,v.conj())
parents={
    'singlet':projector([0,1,-1,0]),
    'generic_entangled':projector([1,.3+.4j,-.7,.2j]),
    'aligned_product':projector([0,0,0,1]),
    'independent_product':projector([1,1,0,0]),
    'rank_two_tangent':projector([0,1,-1,0])+projector([0,0,0,1]),
}
output={}
for label,h in parents.items():
    rows=[]
    for L in range(3,8):
        dims=[]
        for periodic in (False,True):
            e=np.linalg.eigvalsh(matrix(h,L,periodic))
            dims.append(int(np.count_nonzero(np.abs(e)<1e-9)))
        rows.append({'L':L,'open_dimension':dims[0],'periodic_dimension':dims[1]})
    output[label]=rows
with open('ff_checks_results.json','w') as f:
    json.dump(output,f,indent=2)
print(json.dumps(output,indent=2))
