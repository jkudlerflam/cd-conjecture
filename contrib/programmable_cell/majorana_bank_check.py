"""Dressed Bell-pair "Majorana bank" cell for the first-order route (Claude's check).

Cell registers: left pair (port pL, partner sL), right pair (pR, sR), logical Clifford
module of 2n Majoranas gamma_1..gamma_2n on n qubits, parity P, eta_a = i gamma_a P.
Pair states: |0> = Omega (maximally entangled), |a> = (F_a x 1)|Omega>, a = 1..q^2-1,
with F_a an orthonormal traceless Hermitian basis, tr(F_a F_b) = q delta_ab.

  H_L = Delta (1 - |Omega><Omega|)_L + g sum_{a<=r} (|a><0| + |0><a|)_L x gamma_a
  H_R = same on the right pair, with eta_a in place of gamma_a.

Claim checked here: H_L alone has an EXACTLY 2^n-fold degenerate ground space for every g,
because B^dag B = r |0><0| x 1 for B = sum_a |a><0| x gamma_a, and its left-port compression
is  Phi_L(F_b) = -(sin 2theta / sqrt r) gamma_b + (sin^2 theta / r) * (Majorana bilinears),
tan 2theta = 2 g sqrt(r) / Delta.  We also test H_L + H_R (same Clifford algebra at both
ends, the C8 joint form) for splitting of the ground space.

Run: python3 majorana_bank_check.py
"""
import numpy as np
from functools import reduce

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
kron = lambda *o: reduce(np.kron, o)


def gellmann(q):
    """Traceless Hermitian basis with tr(F_a F_b) = q delta_ab."""
    F = []
    for j in range(q):
        for k in range(j + 1, q):
            A = np.zeros((q, q), complex); A[j, k] = A[k, j] = 1; F.append(A)
            B = np.zeros((q, q), complex); B[j, k] = -1j; B[k, j] = 1j; F.append(B)
    for l in range(1, q):
        D = np.zeros((q, q), complex)
        D[:l, :l] = np.eye(l); D[l, l] = -l
        F.append(D)
    return [f * np.sqrt(q / np.trace(f @ f).real) for f in F]


def majoranas(n):
    g = []
    for j in range(n):
        pre = [Z] * j; post = [np.eye(2)] * (n - j - 1)
        g.append(kron(*(pre + [X] + post)))
        g.append(kron(*(pre + [Y] + post)))
    P = kron(*([Z] * n))
    return g, P


def pair_states(q, F):
    om = np.eye(q).reshape(q * q) / np.sqrt(q)
    S = [om] + [np.kron(f, np.eye(q)) @ om for f in F]
    return np.array(S).T            # columns: |0>, |1>, ..., |q^2-1>


def cell_terms(q, n, r, Delta, g):
    F = gellmann(q)
    gam, P = majoranas(n)
    eta = [1j * x @ P for x in gam]
    S = pair_states(q, F)
    Dp = q * q
    Pom = np.outer(S[:, 0], S[:, 0].conj())
    Id_pair, Id_log = np.eye(Dp), np.eye(2 ** n)

    def side(ops):
        H = Delta * np.kron(Id_pair - Pom, Id_log)
        for a in range(r):
            K = np.outer(S[:, a + 1], S[:, 0].conj())
            H = H + g * np.kron(K + K.conj().T, ops[a])
        return H                     # acts on pair x logical

    return F, gam, eta, P, side(gam), side(eta)


def check_one_side(q, n, r, Delta=1.0, g=0.2):
    F, gam, eta, P, HL, HR = cell_terms(q, n, r, Delta, g)
    E, V = np.linalg.eigh(HL)
    N = 2 ** n
    theta = 0.5 * np.arctan2(2 * g * np.sqrt(r), Delta)
    G = V[:, :N]
    Dp = q * q
    out = {'split': E[N - 1] - E[0], 'gap': E[N] - E[N - 1],
           'E0_pred': (Delta - np.sqrt(Delta ** 2 + 4 * g * g * r)) / 2, 'E0': E[0]}
    # compression of port operator F_b (port is the first tensor factor of the pair)
    errs, alphas = [], []
    for b in range(r):
        Fb = np.kron(np.kron(F[b], np.eye(q)), np.eye(N))
        C = G.conj().T @ Fb @ G
        # express in the logical basis: G = (pair x logical) frame, find logical image
        # use the exact ground-state isometry Gx = cos|0> - sin r^-1/2 sum_a |a> gamma_a
        alphas.append(-np.sin(2 * theta) / np.sqrt(r))
    # build the analytic isometry and compare projectors, then compress in that frame
    S = pair_states(q, F)
    Gx = np.cos(theta) * np.kron(S[:, [0]], np.eye(N))
    for a in range(r):
        Gx = Gx - np.sin(theta) / np.sqrt(r) * np.kron(S[:, [a + 1]], gam[a])
    out['frame_err'] = np.linalg.norm(Gx @ Gx.conj().T - G @ G.conj().T)
    lin, bil = [], []
    for b in range(r):
        Fb = np.kron(np.kron(F[b], np.eye(q)), np.eye(N))
        C = Gx.conj().T @ Fb @ Gx
        lin.append(np.linalg.norm(C - alphas[b] * gam[b]) / np.linalg.norm(alphas[b] * gam[b]))
        # residual should be the bilinear term only: check it is parity even
        R = C - alphas[b] * gam[b]
        bil.append(np.linalg.norm(R - P @ R @ P))
    out['alpha'] = alphas[0]
    out['rel_resid'] = max(lin)
    out['resid_parity_odd_part'] = max(bil)
    return out


def check_two_sides(q, n, r, Delta=1.0, g=0.2):
    F, gam, eta, P, HL, HR = cell_terms(q, n, r, Delta, g)
    Dp, N = q * q, 2 ** n
    # full cell: pairL x pairR x logical ; HL acts on pairL x logical, HR on pairR x logical
    def embed_L(H):
        H4 = H.reshape(Dp, N, Dp, N)
        return np.einsum('aibj,cd->acibdj', H4, np.eye(Dp)).reshape(Dp * Dp * N, Dp * Dp * N)
    def embed_R(H):
        H4 = H.reshape(Dp, N, Dp, N)
        return np.einsum('cidj,ab->acibdj', H4, np.eye(Dp)).reshape(Dp * Dp * N, Dp * Dp * N)
    H = embed_L(HL) + embed_R(HR)
    E = np.linalg.eigvalsh(H)
    e1 = (Delta - np.sqrt(Delta ** 2 + 4 * g * g * r)) / 2
    return E[:N + 2] - 2 * e1


if __name__ == '__main__':
    for (q, n, r) in [(2, 1, 2), (2, 2, 3), (3, 2, 4), (3, 4, 8)]:
        for g in (0.05, 0.2, 1.0):
            o = check_one_side(q, n, r, g=g)
            print(f'q={q} n={n} r={r} g={g}: split={o["split"]:.1e} gap={o["gap"]:.3f} '
                  f'E0-E0pred={o["E0"] - o["E0_pred"]:.1e} frame_err={o["frame_err"]:.1e} '
                  f'alpha={o["alpha"]:.4f} rel.resid={o["rel_resid"]:.3f} odd part of resid={o["resid_parity_odd_part"]:.1e}')
    print('\nTwo ports on the same Clifford algebra (H_L + H_R), levels minus 2*E0(one side):')
    for (q, n, r) in [(2, 1, 2), (2, 2, 3), (3, 2, 4)]:
        for g in (0.05, 0.1, 0.2):
            lv = check_two_sides(q, n, r, g=g)
            N = 2 ** n
            print(f'q={q} n={n} r={r} g={g}: lowest {N} levels shifted: {np.round(lv[:N], 7)}, next: {np.round(lv[N:N+2], 4)}')


# ---------------------------------------------------------------------------
# Joint cell (both ends on one Clifford algebra) with GRADED port channels:
# F_a = off-diagonal Gell-Mann matrices on edges between {0..q+-1} and {q+..q-1}.
# Checks exact degeneracy, Phi_L(F_b) = aL gamma_b and Phi_R(F_b) = aR eta_b with no
# residual, and that the first-order logical bond is aL*aR * sum_b eta_b (x) gamma_b.
def graded_channels(q, qp):
    F = []
    for j in range(qp):
        for k in range(qp, q):
            A = np.zeros((q, q), complex); A[j, k] = A[k, j] = 1; F.append(A)
            B = np.zeros((q, q), complex); B[j, k] = -1j; B[k, j] = 1j; F.append(B)
    return [f * np.sqrt(q / np.trace(f @ f).real) for f in F]


def joint_graded(q, qp, n, Delta=1.0, g=0.2):
    F = graded_channels(q, qp)
    r = len(F)
    gam, P = majoranas(n)
    gam = gam[:r]
    eta = [1j * x @ P for x in gam]
    S = pair_states(q, F)            # columns |0>, |F_1 Omega>, ..., |F_r Omega>
    Dp, N = q * q, 2 ** n
    Pom = np.outer(S[:, 0], S[:, 0].conj())

    def side(ops):
        H = Delta * np.kron(np.eye(Dp) - Pom, np.eye(N))
        for a in range(r):
            K = np.outer(S[:, a + 1], S[:, 0].conj())
            H = H + g * np.kron(K + K.conj().T, ops[a])
        return H
    HL, HR = side(gam), side(eta)
    eL = lambda H: np.einsum('aibj,cd->acibdj', H.reshape(Dp, N, Dp, N), np.eye(Dp)).reshape(Dp * Dp * N, -1)
    eR = lambda H: np.einsum('cidj,ab->acibdj', H.reshape(Dp, N, Dp, N), np.eye(Dp)).reshape(Dp * Dp * N, -1)
    H = eL(HL) + eR(HR)
    H = (H + H.conj().T) / 2
    import scipy.linalg as sla
    E, V = sla.eigh(H, driver='evr')
    G = V[:, :N]
    # logical frame: fix it by the |0>_L|0>_R component, which is c0 * psi for every psi
    blk = G.reshape(Dp, Dp, N, N)[0, 0]          # (logical, code-index)
    U, s, Wh = np.linalg.svd(blk)
    Gfix = G @ (U @ Wh).conj().T                 # now Gfix[0,0,:,:] proportional to identity
    out = {'r': r, 'split': E[N - 1] - E[0], 'gap': E[N] - E[N - 1]}
    aL, aR, resL, resR = [], [], [], []
    for b in range(r):
        FL = np.kron(np.kron(np.kron(F[b], np.eye(q)), np.eye(Dp)), np.eye(N))
        FR = np.kron(np.kron(np.eye(Dp), np.kron(F[b], np.eye(q))), np.eye(N))
        CL = Gfix.conj().T @ FL @ Gfix
        CR = Gfix.conj().T @ FR @ Gfix
        al = np.trace(gam[b] @ CL).real / N
        ar = np.trace(eta[b] @ CR).real / N
        aL.append(al); aR.append(ar)
        resL.append(np.linalg.norm(CL - al * gam[b])); resR.append(np.linalg.norm(CR - ar * eta[b]))
    out['alphaL'] = (min(aL), max(aL)); out['alphaR'] = (min(aR), max(aR))
    out['resL'] = max(resL); out['resR'] = max(resR)
    return out


if __name__ == '__main__':
    print('\nJoint graded cells (one Clifford algebra at both ends):')
    for (q, qp, n) in [(2, 1, 1), (3, 1, 2), (4, 2, 4)]:
        for g in (0.1, 0.5):
            o = joint_graded(q, qp, n, g=g)
            print(f"q={q} q+={qp} r={o['r']} n={n} g={g}: split={o['split']:.1e} gap={o['gap']:.3f} "
                  f"alphaL in [{o['alphaL'][0]:.5f},{o['alphaL'][1]:.5f}] alphaR in [{o['alphaR'][0]:.5f},{o['alphaR'][1]:.5f}] "
                  f"resid L={o['resL']:.1e} R={o['resR']:.1e}")
