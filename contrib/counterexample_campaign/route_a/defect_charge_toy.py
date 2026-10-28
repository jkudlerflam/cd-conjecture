"""Tagged material-bond defects in the exact periodic vacancy charge gas.

Gap variables g_j count vacancies after material particle j. In total
physical momentum zero their hopping transfers one vacancy between adjacent
gaps, with matrix element -t (no bosonic square-root factors).
The potential J sum_{j in defects} 1[g_j=0] is an exact tagged-contact toy.
It probes charge relaxation around defects; it is not the full d=24 model.
"""
from pathlib import Path
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import itertools
import json
import math
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent


def compositions(total, length):
    # Stars and bars, with no recursive generator overhead.
    for separators in itertools.combinations(range(total+length-1), length-1):
        p = (-1,) + separators + (total+length-1,)
        yield tuple(p[j+1]-p[j]-1 for j in range(length))


def build(L, N):
    gaps = list(compositions(L-N, N))
    lookup = {g:i for i,g in enumerate(gaps)}
    dim = len(gaps)
    rows, cols = [], []
    for i,g in enumerate(gaps):
        for j in range(N):
            if not g[j]:
                continue
            for dest in ((j-1) % N, (j+1) % N):
                h = list(g); h[j] -= 1; h[dest] += 1
                rows.append(lookup[tuple(h)]); cols.append(i)
    H = sp.csr_matrix((-np.ones(len(rows)), (rows,cols)), shape=(dim,dim))
    assert (H-H.T).nnz == 0
    G = np.array(gaps, dtype=np.int16)
    contact = (G == 0).astype(float)
    # The positive Slater-Vandermonde ground state in relative coordinates.
    x = np.zeros((dim,N), dtype=np.int16)
    x[:,1:] = np.cumsum(G[:,:-1]+1, axis=1)
    logpsi = np.zeros(dim)
    for i in range(N):
        for j in range(i+1,N):
            logpsi += np.log(2*np.sin(np.pi*(x[:,j]-x[:,i])/L))
    psi = np.exp(logpsi-logpsi.max()); psi /= np.linalg.norm(psi)
    E0 = -2*math.sin(math.pi*N/L)/math.sin(math.pi/L)
    residual = float(np.linalg.norm(H@psi-E0*psi))
    assert residual < 1e-10, (L,residual)
    alpha = N/L - math.sin(math.pi*N/L)**2/(N*L*math.sin(math.pi/L)**2)
    alpha_ed = (psi*psi)@contact
    assert np.max(abs(alpha_ed-alpha)) < 1e-11
    return H, psi, E0, contact, alpha, residual


def case(L, N):
    start = time.monotonic()
    H, psi, E0, contact, alpha, residual = build(L,N)
    dim = len(psi)
    b = (contact[:,0]-alpha)*psi
    def matvec(x):
        return H@x-E0*x+psi*np.dot(psi,x)
    A = spla.LinearOperator((dim,dim), matvec=matvec, dtype=float)
    x, info = spla.cg(A,b,rtol=2e-11,atol=1e-13,maxiter=1500)
    assert info == 0, (L,info)
    susceptibility = (contact-alpha).T@(psi*x)
    response_residual = float(np.linalg.norm(matvec(x)-b))
    assert response_residual < 1e-9
    rho = N/L
    aprime = (-rho*rho+rho*math.sin(2*math.pi*rho)/math.pi
              -math.sin(math.pi*rho)**2/math.pi**2)
    eprime2 = 2*math.pi*rho**3*math.sin(math.pi*rho)
    pedestal = aprime**2/(2*eprime2)
    runs = []
    for sites, couplings in [([0],[.04,.02,.01]),
                            ([0,1],[.02]),
                            (list(range(0,N,2)),[.02]),
                            (list(range(N)),[.02])]:
        potential = contact[:,sites].sum(axis=1)
        for J in couplings:
            Adef = H+sp.diags(J*potential)
            w,U = spla.eigsh(Adef,k=1,which="SA",v0=psi,tol=2e-12,maxiter=2500)
            E = float(w[0]); state=U[:,0]
            err = float(np.linalg.norm(Adef@state-E*state))
            assert err < 1e-9
            delta = E-E0; leading = J*len(sites)*alpha
            runs.append(dict(defect_sites=sites,J=J,energy_cost=delta,
                first_order=leading,relative_cost=delta/leading,
                relaxation_per_defect_over_J2=(leading-delta)/(len(sites)*J*J),
                eigenpair_residual=err))
    result = dict(L=L,N=N,density=N/L,dimension=dim,
        alpha_exact=alpha,free_energy=E0,free_state_residual=residual,
        susceptibility_kernel=susceptibility.tolist(),
        susceptibility_abs_row_sum=float(np.sum(abs(susceptibility))),
        susceptibility_abs_first_moment=float(sum(
            min(j,N-j)*abs(susceptibility[j]) for j in range(N))),
        canonical_pedestal_coefficient=pedestal,
        susceptibility_pedestal_subtracted_abs_first_moment=float(sum(
            min(j,N-j)*abs(susceptibility[j]+pedestal/N) for j in range(N))),
        susceptibility_local=float(susceptibility[0]),
        response_residual=response_residual,runs=runs,
        seconds=time.monotonic()-start)
    print(json.dumps({k:v for k,v in result.items() if k not in ("runs","susceptibility_kernel")}), flush=True)
    return result


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-L", type=int, default=20)
    args = ap.parse_args()
    results = []
    for L in range(8,args.max_L+1,2):
        results.append(case(L,L//2))
        (HERE/"defect_charge_results.json").write_text(json.dumps(results,indent=2)+"\n")
