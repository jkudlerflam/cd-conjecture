"""C86 scalar mechanism tests. NumPy/SciPy, local outputs only."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import splu


def coefficient(H, left, right, active=None):
    s = len(H)
    if active is None:
        active = list(range(s))
    n = len(active)
    lu = splu(sp.kronsum(sp.csc_matrix(H), sp.csc_matrix(H), format='csc'))
    x = lu.solve(np.outer(right, left).ravel()).reshape(s, s)
    U = np.zeros((s*s, n), dtype=complex)
    V = U.copy()
    for j, pos in enumerate(active):
        U[pos*s:(pos+1)*s, j] = x[pos, :]
        V[pos::s, j] = x[:, pos]
    C = U.conj().T @ lu.solve(V)
    A = U.conj().T @ lu.solve(U) + V.conj().T @ lu.solve(V)
    return C.real, -2j*A.imag


def path(s, mass=4., end_mass=0.):
    H = mass*np.eye(s, dtype=complex)
    for j in range(s-1):
        H[j, j+1] = H[j+1, j] = -1
    H[-1, -1] += end_mass
    return H


def main():
    records = []
    for n in [2, 3, 4, 6, 9, 10]:
        for end_mass in [0., 1., 10., 1000.]:
            H = path(n+1, end_mass=end_mass)
            F, M = coefficient(H, np.eye(n+1)[0], np.eye(n+1)[-1], list(range(n)))
            S = (F+F.T)/2
            ev = la.eigvalsh(S)
            records.append(dict(status='numerical',family='single_rail',n=n,
                vertices=n+1,end_mass=end_mass,
                hermitian_eigenvalues_scaled=(ev/la.norm(S, 2)).tolist(),
                onsite_norm=float(la.norm(M)),
                diagonal_positive=bool(np.all(np.diag(S)>0))))
    # A reflection-invariant graph with complex edges and nontrivial loops.
    reflection_checks = []
    for s in [4, 5, 6, 10]:
        J = np.eye(s)[::-1]
        H = path(s, mass=6)
        for i in range(s-2):
            H[i, i+2] += (1+1j*(2*i-s+3))/8
            H[i+2, i] = H[i, i+2].conjugate()
        H = (H+J@H@J)/2
        assert la.eigvalsh(H)[0]>0 and la.norm(H-J@H@J)<1e-12
        F, M = coefficient(H, np.eye(s)[0], np.eye(s)[-1])
        G = F@J
        assert la.norm(G-G.T)<1e-12
        reflection_checks.append(dict(status='numerical',vertices=s,
            reflection_residual=float(la.norm(H-J@H@J)),
            gram_symmetry_residual=float(la.norm(G-G.T)),
            gram_min_eigenvalue=float(la.eigvalsh((G+G.T)/2)[0]),
            hermitian_inertia=[int(sum(la.eigvalsh((F+F.T)/2)>1e-16)),
                              int(sum(la.eigvalsh((F+F.T)/2)<-1e-16))],
            onsite_norm=float(la.norm(M)), onsite_reflection_residual=float(la.norm(M-J@M@J))))
    port_checks=[]
    amplitudes=[0, 1, -1, 1j, -1j, (1+1j)/2, (1-1j)/2]
    for n in [9,10]:
        H=path(n+1)
        for c in amplitudes:
            for d in amplitudes:
                left=np.eye(n+1,dtype=complex)[0];left[1]=c
                right=np.eye(n+1,dtype=complex)[-1];right[-2]=d
                F,M=coefficient(H,left,right,list(range(n)))
                S=(F+F.T)/2;ev=la.eigvalsh(S)
                port_checks.append(dict(status='numerical',n=n,
                    left_amplitude=[float(c.real),float(c.imag)],
                    right_amplitude=[float(d.real),float(d.imag)],
                    min_scaled=float(ev[0]/la.norm(S,2)),
                    max_scaled=float(ev[-1]/la.norm(S,2)),
                    endpoint_minor_determinant=float(S[0,0]*S[-1,-1]-S[0,-1]**2)))
    json.dump(dict(single_rail=records,reflection=reflection_checks,two_site_sources=port_checks),open('exploration.json','w'),indent=2)
    print(json.dumps(dict(single_rail_cases=len(records),reflection=reflection_checks,
        two_site_source_cases=len(port_checks),
        two_site_source_indefinite=sum(r['min_scaled']<0<r['max_scaled'] for r in port_checks))))


if __name__ == '__main__':
    main()
