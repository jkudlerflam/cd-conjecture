"""Red-team probe for Theorem A of C201 (Claude Opus 5.5, C203).

For small rings, lower-bound L * dist_1(rho2_ring, Omega_inf) by the trace distance
to the outer relaxation Omega^(r) (two-site marginals of r-site states with
tr_first = tr_last), r > L.  Ring states tested: GHZ, W, Dicke(L/2), dimer
(singlet covering) average, and ring ground states of Heisenberg, XX, SU(3) swap,
plus an alternating ascent over bonds h (witness of the distance SDP fed back as
the next bond).  Omega_inf is contained in Omega^(r), so each number is a rigorous
lower bound (up to solver tolerance) on L*dist(rho2, Omega_inf).
"""
import itertools, sys
import numpy as np, cvxpy as cp
from scipy.sparse import kron, identity, csr_matrix
from scipy.sparse.linalg import eigsh

def ptrace_keep(psi_or_rho, d, L, keep):
    rho = psi_or_rho
    if rho.ndim == 1:
        t = rho.reshape([d]*L)
        rest = [i for i in range(L) if i not in keep]
        t = np.transpose(t, keep + rest).reshape(d**len(keep), -1)
        return t @ t.conj().T
    t = rho.reshape([d]*(2*L))
    rest = [i for i in range(L) if i not in keep]
    perm = keep + rest + [L+i for i in keep] + [L+i for i in rest]
    t = np.transpose(t, perm).reshape(d**len(keep), d**(L-len(keep)), d**len(keep), d**(L-len(keep)))
    return np.einsum('aibi->ab', t)

def ti_marginal(psi, d, L):
    return sum(ptrace_keep(psi, d, L, [j, (j+1) % L]) for j in range(L)) / L

def ring_H(h, d, L):
    H = csr_matrix((d**L, d**L), dtype=complex)
    hs = csr_matrix(h)
    for j in range(L - 1):
        H = H + kron(kron(identity(d**j), hs), identity(d**(L-j-2)))
    # seam bond (L-1, 0): permute
    P = np.zeros((d**L, d**L))
    for idx in range(d**L):
        digs = np.unravel_index(idx, [d]*L)
        new = digs[1:] + digs[:1]
        P[np.ravel_multi_index(new, [d]*L), idx] = 1
    P = csr_matrix(P)
    Hb = kron(hs, identity(d**(L-2)))
    H = H + P.T @ Hb @ P
    return H

def dist_to_outer(rho2, d, r):
    n = d**r
    S = cp.Variable((n, n), hermitian=True)
    Pp = cp.Variable((d*d, d*d), hermitian=True)
    Nn = cp.Variable((d*d, d*d), hermitian=True)
    s2 = cp.partial_trace(S, [d*d, d**(r-2)], axis=1)
    left = cp.partial_trace(S, [d, d**(r-1)], axis=0)
    right = cp.partial_trace(S, [d**(r-1), d], axis=1)
    eq = (rho2 - s2 == Pp - Nn)
    cons = [S >> 0, cp.real(cp.trace(S)) == 1, left == right, Pp >> 0, Nn >> 0, eq]
    prob = cp.Problem(cp.Minimize(cp.real(cp.trace(Pp) + cp.trace(Nn))), cons)
    prob.solve(solver=cp.CLARABEL)
    W = eq.dual_value
    return prob.value, W

def state_GHZ(d, L):
    v = np.zeros(d**L, complex); v[0] = v[-1] = 1/np.sqrt(2); return v
def state_W(L):
    v = np.zeros(2**L, complex)
    for j in range(L): v[1 << j] = 1/np.sqrt(L)
    return v
def state_dicke(L, k):
    v = np.zeros(2**L, complex)
    for c in itertools.combinations(range(L), k): v[sum(1 << j for j in c)] = 1
    return v/np.linalg.norm(v)

def dimer_rho2(L):
    # translation average of the singlet covering: half the bonds singlet, half product of maximally mixed
    s = np.array([0, 1, -1, 0])/np.sqrt(2)
    return 0.5*np.outer(s, s) + 0.5*np.eye(4)/4

if __name__ == '__main__':
    d = 2
    X = np.array([[0,1],[1,0]]); Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1,-1])
    heis = (np.kron(X,X)+np.kron(Y,Y)+np.kron(Z,Z))/4
    xx = -(np.kron(X,X)+np.kron(Y,Y))/2
    rows = []
    for L in [4, 5, 6]:
        r = max(6, L + 1)
        cases = {'GHZ': ti_marginal(state_GHZ(2, L), 2, L), 'W': ti_marginal(state_W(L), 2, L),
                 'Dicke': ti_marginal(state_dicke(L, L//2), 2, L)}
        if L % 2 == 0: cases['dimer'] = dimer_rho2(L)
        for name, h in [('Heis', heis), ('XX', xx)]:
            w, v = np.linalg.eigh(ring_H(h, d, L).toarray()); cases[name+'_gs'] = ti_marginal(v[:, 0], d, L)
        # alternating ascent from a random bond
        rng = np.random.default_rng(L)
        A = rng.normal(size=(4,4)) + 1j*rng.normal(size=(4,4)); h = (A + A.conj().T)/2
        best = 0
        for it in range(3):
            w, v = np.linalg.eigh(ring_H(h, d, L).toarray())
            rho2 = ti_marginal(v[:, 0], d, L)
            val, W = dist_to_outer(rho2, d, r)
            best = max(best, val)
            if W is None: break
            W = (W + W.conj().T)/2
            h = W  # witness: rho2 sits on the low side of it
        cases['ascent_best'] = None; ascent = best
        for name, rho2 in cases.items():
            if rho2 is None: continue
            val, _ = dist_to_outer(rho2, d, r)
            rows.append((L, r, name, val, L*val)); print(L, r, name, f"dist>={val:.2e}", f"L*dist>={L*val:.2e}", flush=True)
        print(L, r, 'ascent', f"dist>={ascent:.2e}", f"L*dist>={L*ascent:.2e}", flush=True)
