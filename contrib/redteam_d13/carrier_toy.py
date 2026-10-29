"""Attack (i)/(iv): is the mobile payload amplitude alpha times the frozen one?

Toy of the virtual rail traversal that produces E_unit.  One carrier excitation
rides on material letters 0..m-1 of a gap ring of N letters with N_v = N vacancies
(half occupancy).  In C88 every material term, hops AND the D-scale diagonal, is a
contact term, so
    H = t H_gap + J [ -sum_k c_k (|k+1><k| + h.c.) + (D/2) sum_i (c_{i-1}+c_i)|i><i| ].
The traversal amplitude is A_r = <vac, r|(H - E_vac)^-1|vac, 0>.  The fast-vacancy
dictionary used by C88/C93 replaces every c_k by alpha, giving
A_r^avg = <r|(J alpha h)^-1|0>.  We measure the per-step deviation
    eps_step = (A_{r+1}/A_r)/(A^avg_{r+1}/A^avg_r) - 1
and the overall deviation A_0/A_0^avg - 1 as functions of D J/t.
"""
import json
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

sys.path.insert(0, '.')
from gapfluid import compositions, hop_matrix, ground

N, Nv, m = 8, 8, 5
G = compositions(N, Nv)
Hg = hop_matrix(G)
Ev, Vv = ground(Hg)
psi = Vv[:, 0]
c = [(G[:, j] == 0).astype(float) for j in range(N)]
alpha = float(psi ** 2 @ c[0])
M = len(G)

res = []
for D in (10.0, 30.0, 100.0):
    # frozen-averaged carrier on the segment
    h = np.zeros((m, m))
    for i in range(m):
        h[i, i] = D
        if i + 1 < m:
            h[i, i + 1] = h[i + 1, i] = -1.0
    for t_over_J in (1.0, 10.0, 100.0, 1e3, 1e4, 1e5):
        J, t = 1.0, t_over_J
        Havg = J * alpha * h
        Aavg = np.linalg.solve(Havg, np.eye(m)[:, 0])
        # full model on (gap config) x (carrier position)
        Hc = sp.kron(t * Hg, sp.identity(m), format='csr')
        for i in range(m):
            P = sp.csr_matrix(([1.0], ([i], [i])), shape=(m, m))
            cl = c[(i - 1) % N]
            cr = c[i]
            Hc = Hc + J * (D / 2) * sp.kron(sp.diags(cl + cr), P, format='csr')
            if i + 1 < m:
                Hop = sp.csr_matrix(([-1.0, -1.0], ([i, i + 1], [i + 1, i])), shape=(m, m))
                Hc = Hc + J * sp.kron(sp.diags(c[i]), Hop, format='csr')
        Hs = Hc - t * Ev[0] * sp.identity(M * m, format='csr')
        rhs = np.kron(psi, np.eye(m)[:, 0])
        x, info = spl.cg(Hs, rhs, rtol=1e-12, maxiter=200000)
        assert info == 0
        lmin = spl.eigsh(Hs, k=1, which='SA', tol=1e-8)[0][0]
        A = np.array([np.kron(psi, np.eye(m)[:, r]) @ x for r in range(m)])
        ratio = A / Aavg
        eps_step = [float((A[r + 1] / A[r]) / (Aavg[r + 1] / Aavg[r]) - 1) for r in range(m - 1)]
        rec = dict(lowest_eigenvalue_of_H_minus_Evac=float(lmin), D=D, t_over_J=t_over_J, DJ_over_t=D / t_over_J, alpha=alpha,
                   ratio_r=ratio.tolist(), eps_step=eps_step,
                   eps_step_over_DJt=[e / (D / t_over_J) for e in eps_step])
        print(D, t_over_J, round(float(lmin),4), np.round(ratio, 5), np.round(eps_step, 6), flush=True)
        res.append(rec)
json.dump(dict(N=N, Nv=Nv, m=m, results=res), open('carrier_toy.json', 'w'), indent=1)
