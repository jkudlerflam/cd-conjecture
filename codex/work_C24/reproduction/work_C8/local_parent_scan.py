"""Linear parent-Hamiltonian search for the smallest explicit Choi code.
A negative result is about these six code embeddings only.
"""
import itertools
import json
import numpy as np
from scipy.linalg import svd
from c8_checks import tensor
from c8_resume_checks import joint_code


def constraints(ops,V):
    D=V.shape[1]
    columns=[]
    for H in ops:
        HV=H@V
        energy=np.trace(V.conj().T@HV)/D
        residual=HV-energy*V
        columns.append(np.concatenate([residual.real.ravel(),residual.imag.ravel()]))
    return np.array(columns).T


def main():
    q,m=2,2
    *_,V,kl,kr=joint_code(q,m)
    ps=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    L=5
    nonuniform=[]
    for j in range(L):
        for a in range(1,4):
            ids=[0]*L; ids[j]=a
            nonuniform.append(tensor([ps[x] for x in ids]))
    for j in range(L-1):
        for a,b in itertools.product(range(1,4),repeat=2):
            ids=[0]*L; ids[j]=a; ids[j+1]=b
            nonuniform.append(tensor([ps[x] for x in ids]))
    uniform=[]
    # Exclude only the identity bond. Includes oriented one-site terms.
    for a,b in itertools.product(range(4),repeat=2):
        if a==b==0: continue
        H=np.zeros((2**L,2**L),complex)
        for j in range(L-1):
            ids=[0]*L; ids[j]=a; ids[j+1]=b
            H+=tensor([ps[x] for x in ids])
        uniform.append(H)
    rows=[]
    # Choi order is L,R,envL,envR,logical. Keep the ports at the endpoints.
    for middle in itertools.permutations([2,3,4]):
        order=[0,*middle,1]
        vp=V.reshape([2]*L+[2]).transpose(order+[L]).reshape(2**L,2)
        row={'order':order}
        for name,ops in [('uniform',uniform),('nonuniform',nonuniform)]:
            C=constraints(ops,vp)
            _,s,vh=svd(C,full_matrices=False)
            null=vh[s<1e-10]
            row[name]={'number_of_parameters':len(ops),'nullity':len(null),
                       'smallest_singular_values':[float(x) for x in s[-6:]]}
            # Does any null direction actually act nontrivially outside the code?
            if len(null):
                mats=[sum(c*H for c,H in zip(vec,ops)) for vec in null]
                eigs=[np.linalg.eigvalsh(H) for H in mats]
                row[name]['null_basis_spectra']=[e.tolist() for e in eigs]
                row[name]['null_basis_coefficients']=[vec.tolist() for vec in null]
                row[name]['null_code_energies']=[float(np.trace(vp.conj().T@H@vp).real/2) for H in mats]
        rows.append(row)
    with open('local_parent_scan.json','w') as f: json.dump(rows,f,indent=2)
    for row in rows:
        print(row['order'], {n:{k:row[n][k] for k in ['nullity','smallest_singular_values']}
                             for n in ['uniform','nonuniform']})


if __name__=='__main__': main()
