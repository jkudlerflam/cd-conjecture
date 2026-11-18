"""Referee implementation of the README's scalar formula for T, mu^L, mu^R,
written from the README text (section "Second order"), not from scalar9.py.

Vertices of the token graph of one species: X_0..X_n (index j), D_0..D_n (index n+1+j).
Orientation of the complex edges (the README does not state it) is taken from the
bond: a token moving one step to the right on the express rail has amplitude
-e^{i phi_s}, and the closing link  -t_d e^{i psi_s} |X_s y><s y|  takes D_n to X_n.

Two back ends: numpy (double) and mpmath (arbitrary precision), both by dense
Kronecker-sum linear solves -- deliberately different from scalar9's
eigen-decomposition Sylvester solver.
"""
import numpy as np


def Hs_np(n, a, D, phi, td, psi):
    N = 2 * n + 2
    H = np.zeros((N, N), complex)
    for j in range(n + 1):
        H[j, j] = a + D
        H[n + 1 + j, n + 1 + j] = a
    for j in range(n):
        H[j + 1, j] = -np.exp(1j * phi); H[j, j + 1] = -np.exp(-1j * phi)
        H[n + 2 + j, n + 1 + j] = -1; H[n + 1 + j, n + 2 + j] = -1
    H[0, n + 1] = -1; H[n + 1, 0] = -1
    H[n, 2 * n + 1] += -td * np.exp(1j * psi); H[2 * n + 1, n] += -td * np.exp(-1j * psi)
    return H


def blocks_np(n, a, D, td, phiA=None, psiA=np.pi, psiB=np.pi / 2):
    if phiA is None:
        phiA = np.pi / (2 * (n + 1))
    N = 2 * n + 2
    H = {"A": Hs_np(n, a, D, phiA, td, psiA), "B": Hs_np(n, a, D, 0.0, td, psiB)}
    I = np.eye(N)
    Rinv = {(s, t): np.kron(H[s], I) + np.kron(I, H[t]) for s in "AB" for t in "AB"}
    R = {k: np.linalg.inv(v) for k, v in Rinv.items()}
    src = np.zeros(N * N, complex); src[n * N + (n + 1)] = 1.0        # |X_n> x |D_0>
    def PiL(j):
        P = np.zeros((N, N)); P[n + 1 + j, n + 1 + j] = 1; return np.kron(P, I)
    def PiR(j):
        P = np.zeros((N, N)); P[n + 1 + j, n + 1 + j] = 1; return np.kron(I, P)
    T = np.zeros((n, n), complex); muL = np.zeros((n, n), complex); muR = np.zeros((n, n), complex)
    bra_AB = src.conj() @ R[("A", "B")]
    bra_BA = src.conj() @ R[("B", "A")]
    ket_BA = R[("B", "A")] @ src
    ket_AB = R[("A", "B")] @ src
    for j in range(n):
        for jp in range(n):
            T[j, jp] = 1j * (bra_AB @ PiL(j) @ R[("B", "B")] @ PiR(jp) @ ket_BA
                             + bra_AB @ PiR(jp) @ R[("A", "A")] @ PiL(j) @ ket_BA)
            muL[j, jp] = -(bra_AB @ PiL(j) @ R[("B", "B")] @ PiL(jp) @ ket_AB
                           - bra_BA @ PiL(jp) @ R[("A", "A")] @ PiL(j) @ ket_BA)
            # "mu^R likewise on the right factor": the right token has the opposite species,
            # so on the right factor the roles of AB and BA are exchanged (referee's reading)
            muR[j, jp] = -(bra_BA @ PiR(j) @ R[("B", "B")] @ PiR(jp) @ ket_BA
                           - bra_AB @ PiR(jp) @ R[("A", "A")] @ PiR(j) @ ket_AB)
    return T, muL, muR


def blocks_mp(n, a, D, td, mp, phiA=None, psiA=None, psiB=None):
    """same formula in mpmath; Kronecker-sum systems solved by LU (dimension (2n+2)^2)."""
    if phiA is None: phiA = mp.pi / (2 * (n + 1))
    if psiA is None: psiA = mp.pi
    if psiB is None: psiB = mp.pi / 2
    N = 2 * n + 2
    def Hs(phi, psi):
        H = mp.zeros(N, N)
        for j in range(n + 1):
            H[j, j] = a + D; H[n + 1 + j, n + 1 + j] = a
        for j in range(n):
            H[j + 1, j] = -mp.expj(phi); H[j, j + 1] = -mp.expj(-phi)
            H[n + 2 + j, n + 1 + j] = -1; H[n + 1 + j, n + 2 + j] = -1
        H[0, n + 1] = -1; H[n + 1, 0] = -1
        H[n, 2 * n + 1] += -td * mp.expj(psi); H[2 * n + 1, n] += -td * mp.expj(-psi)
        return H
    H = {"A": Hs(phiA, psiA), "B": Hs(mp.mpf(0), psiB)}
    def kronsum(Hs_, Ht_):
        M = mp.zeros(N * N, N * N)
        for i in range(N):
            for k in range(N):
                if Hs_[i, k] != 0:
                    for j in range(N):
                        M[i * N + j, k * N + j] += Hs_[i, k]
        for j in range(N):
            for l in range(N):
                if Ht_[j, l] != 0:
                    for i in range(N):
                        M[i * N + j, i * N + l] += Ht_[j, l]
        return M
    LU = {}
    for s in "AB":
        for t in "AB":
            LU[(s, t)] = mp.mp.LU_decomp(kronsum(H[s], H[t]), overwrite=True, use_cache=False)   # factor once
    def solve(st, b):
        A_, p_ = LU[st]
        return mp.mp.U_solve(A_, mp.mp.L_solve(A_, b, p_))
    src = mp.zeros(N * N, 1); src[n * N + (n + 1)] = 1
    ket_BA = solve(("B", "A"), src); ket_AB = solve(("A", "B"), src)
    # R is Hermitian: <src| R_st = (R_st |src>)^dagger
    def projL(vec, j):
        out = mp.zeros(N * N, 1)
        for c in range(N):
            out[(n + 1 + j) * N + c] = vec[(n + 1 + j) * N + c]
        return out
    def projR(vec, j):
        out = mp.zeros(N * N, 1)
        for r in range(N):
            out[r * N + (n + 1 + j)] = vec[r * N + (n + 1 + j)]
        return out
    def dot(u, v):
        return mp.fsum(mp.conj(u[i]) * v[i] for i in range(N * N))
    T = mp.zeros(n, n); muL = mp.zeros(n, n); muR = mp.zeros(n, n)
    BB_R_BA = [solve(("B", "B"), projR(ket_BA, jp)) for jp in range(n)]
    AA_L_BA = [solve(("A", "A"), projL(ket_BA, j)) for j in range(n)]
    BB_L_AB = [solve(("B", "B"), projL(ket_AB, jp)) for jp in range(n)]
    AA_R_AB = [solve(("A", "A"), projR(ket_AB, j)) for j in range(n)]
    for j in range(n):
        for jp in range(n):
            T[j, jp] = 1j * (dot(projL(ket_AB, j), BB_R_BA[jp]) + dot(projR(ket_AB, jp), AA_L_BA[j]))
            muL[j, jp] = -(dot(projL(ket_AB, j), BB_L_AB[jp]) - dot(projL(ket_BA, jp), AA_L_BA[j]))
            muR[j, jp] = -(dot(projR(ket_BA, j), BB_R_BA[jp]) - dot(projR(ket_AB, jp), AA_R_AB[j]))
    return T, muL, muR
