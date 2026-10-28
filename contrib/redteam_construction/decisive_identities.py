"""Independent rechecks of the identities the construction rests on (red team).

1. G1 (contrib/unbounded_c/G1, Eq. 20): build both blocks from the realization lemma by my
   own code, diagonalize them, and compute the second-order bond from PHYSICAL eigenvectors
   and denominators. Target: H2 = -7177/8192 I + 27/8192 P_3.
2. Dressed Bell-pair closure (programmable_cell Sec. 3B): exact 2^n degeneracy for several
   (q, r, g), including r up to q^2 - 1 (ungraded), and the analytic ground energy.
3. Port-idle identity (programmable_cell Sec. 2): F eta = eta F for a history isometry whose
   gates avoid the port, and failure of the compression identity with one port gate.

Run: OPENBLAS_NUM_THREADS=1 python3 decisive_identities.py   (seconds)
"""
import numpy as np
from functools import reduce

rng = np.random.default_rng(7)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
kron = lambda *o: reduce(np.kron, o)


def psd_sqrt(M):
    w, U = np.linalg.eigh(M)
    return (U * np.sqrt(np.clip(w, 0, None))) @ U.conj().T


def complete_unitary(W):
    """Unitary whose first columns are the orthonormal columns of W."""
    m, k = W.shape
    Q, _ = np.linalg.qr(np.hstack([W, rng.normal(size=(m, m - k)) + 0j]))
    Q[:, :k] = W
    # re-orthonormalize the completion against W
    for j in range(k, m):
        v = Q[:, j] - Q[:, :j] @ (Q[:, :j].conj().T @ Q[:, j])
        Q[:, j] = v / np.linalg.norm(v)
    return Q


def g1_check():
    d, N, r = 2, 3, 3
    K = 11
    D = d * K
    Fp = [X, Y, Z]                                # tr(F_a F_b) = 2 delta
    # code G|i> = d^-1/2 sum_a |a>_port |a,i>_K, K-states |a,i> = a*N + i
    G = np.zeros((D, N), complex)
    for i in range(N):
        for a in range(d):
            G[a * K + (a * N + i), i] = 1 / np.sqrt(d)
    S = np.hstack([np.kron(F, np.eye(K)) @ G for F in Fp])    # D x rN, columns (a, i)
    assert np.allclose(S.conj().T @ S, np.eye(r * N)) and np.allclose(G.conj().T @ S, 0)
    # weights of Eq. (20); syndrome index identified with logical index (r = N = 3)
    I9 = np.eye(9)
    P3 = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            P3[i * 3 + j, j * 3 + i] = 1
    om = np.eye(3).reshape(9) / np.sqrt(3)
    Fm, Jm, eta = P3 - I9 / 3, 3 * np.outer(om, om) - I9 / 3, 9 / 64
    MA = [I9 / 2 + eta * Fm, I9 / 2 - eta * Fm]
    MB = [I9 / 2 - eta * Jm, I9 / 2 + eta * Jm]
    Es = [1.0, 3.0]

    def block(M):
        q0 = D - N
        W = np.vstack([psd_sqrt(M[0]), psd_sqrt(M[1]), np.zeros((q0 - 18, 9))])
        U = complete_unitary(W)
        # orthonormal frame of G-perp whose first rN columns are S
        FQ = complete_unitary(np.hstack([S, G]))[:, :]    # columns: S, G, rest
        FQ = np.hstack([FQ[:, :9], FQ[:, 12:]])           # drop the G columns
        DE = np.diag([Es[0]] * 9 + [Es[1]] * 9 + [Es[1]] * (q0 - 18))
        return FQ @ U.conj().T @ DE @ U @ FQ.conj().T

    HA, HB = block(MA), block(MB)

    def amps(H):
        E, V = np.linalg.eigh(H)
        assert np.allclose(E[:N], 0, atol=1e-12) and E[N] > 0.999
        A = [V[:, N:].conj().T @ (np.kron(F, np.eye(K)) @ G) for F in Fp]
        return E[N:], A
    EA, AA = amps(HA)
    EB, AB = amps(HB)
    H1 = sum(np.kron(G.conj().T @ np.kron(F, np.eye(K)) @ G, G.conj().T @ np.kron(F, np.eye(K)) @ G) for F in Fp)
    H2 = np.zeros((N * N, N * N), complex)
    C = 1.0 / (EA[:, None] + EB[None, :])
    for a in range(r):
        for b in range(r):
            # sum_{e,f} (A_a[e]^dag A_b[e]) x (B_a[f]^dag B_b[f]) C_ef
            T = np.einsum('ei,ej,ef,fk,fl->ikjl', AA[a].conj(), AA[b], C, AB[a].conj(), AB[b])
            H2 -= T.reshape(N * N, N * N)
    Psw = P3
    target = -7177 / 8192 * np.eye(9) + 27 / 8192 * Psw
    print('1. G1 (d,N)=(2,3): |H1| = %.1e, |H2 - target|_F = %.2e, H2 hermitian err %.1e'
          % (np.linalg.norm(H1), np.linalg.norm(H2 - target), np.linalg.norm(H2 - H2.conj().T)))


def bell_pair_closure():
    print('2. Dressed Bell-pair cell, one side: splitting of the 2^n ground space, E0 error')
    for (q, n, r) in [(2, 1, 2), (2, 2, 3), (3, 2, 4), (3, 4, 8), (4, 3, 6), (3, 3, 5)]:
        # Hermitian traceless basis tr = q delta
        F = []
        for j in range(q):
            for k in range(j + 1, q):
                A = np.zeros((q, q), complex); A[j, k] = A[k, j] = 1; F.append(A)
                B = np.zeros((q, q), complex); B[j, k] = -1j; B[k, j] = 1j; F.append(B)
        for l in range(1, q):
            Dg = np.zeros((q, q), complex); Dg[:l, :l] = np.eye(l); Dg[l, l] = -l; F.append(Dg)
        F = [f * np.sqrt(q / np.trace(f @ f).real) for f in F]
        # random orthogonal mixing of the channels: closure must not depend on the basis
        O, _ = np.linalg.qr(rng.normal(size=(len(F), len(F))))
        F = [sum(O[i, j] * F[j] for j in range(len(F))) for i in range(len(F))]
        gam = []
        for j in range(n):
            gam.append(kron(*([Z] * j + [X] + [np.eye(2)] * (n - j - 1))))
            gam.append(kron(*([Z] * j + [Y] + [np.eye(2)] * (n - j - 1))))
        om = np.eye(q).reshape(q * q) / np.sqrt(q)
        Nl = 2 ** n
        for g in (0.07, 0.4, 1.3):
            Delta = 1.0
            H = Delta * np.kron(np.eye(q * q) - np.outer(om, om), np.eye(Nl))
            for a in range(r):
                Kt = np.outer(np.kron(F[a], np.eye(q)) @ om, om)
                H = H + g * np.kron(Kt + Kt.conj().T, gam[a])
            E = np.linalg.eigvalsh(H)
            E0 = (Delta - np.sqrt(Delta ** 2 + 4 * g * g * r)) / 2
            print('   q=%d n=%d r=%d g=%.2f: split %.1e, E0 err %.1e, gap %.4f'
                  % (q, n, r, g, E[Nl - 1] - E[0], E[0] - E0, E[Nl] - E[0]))


def port_idle():
    q, w, T = 2, 3, 6          # port qubit, partner qubit, w workspace qubits, T gates
    dims = [q, q] + [2] * w
    Dt = int(np.prod(dims))

    def rand_u(k):
        A = rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k))
        Q, R = np.linalg.qr(A)
        return Q * (np.diag(R) / abs(np.diag(R)))

    def on(Uloc, first, nsites):
        left = int(np.prod(dims[:first])); right = int(np.prod(dims[first + nsites:]))
        return kron(np.eye(left), Uloc, np.eye(right))
    gates = [on(rand_u(4), 2 + rng.integers(0, w - 1), 2) for _ in range(T)]
    Ws = [np.eye(Dt)]
    for U in gates:
        Ws.append(U @ Ws[-1])
    L = T + 1
    eta = np.vstack([Wt for Wt in Ws]) / np.sqrt(L)       # (L*Dt) x Dt history isometry
    Fp = on(X, 0, 1)
    lhs = np.kron(np.eye(L), Fp) @ eta
    rhs = eta @ Fp
    print('3. Port-idle: |F eta - eta F| = %.1e, |eta^dag eta - 1| = %.1e'
          % (np.linalg.norm(lhs - rhs), np.linalg.norm(eta.conj().T @ eta - np.eye(Dt))))
    # one gate touching the port at step 3, undone at step 4
    Up = on(rand_u(4), 0, 2)
    gates2 = gates[:2] + [Up, Up.conj().T] + gates[2:]
    Ws = [np.eye(Dt)]
    for U in gates2:
        Ws.append(U @ Ws[-1])
    eta2 = np.vstack(Ws) / np.sqrt(len(Ws))
    comp = eta2.conj().T @ np.kron(np.eye(len(Ws)), Fp) @ eta2
    print('   with one port gate (undone a step later): |eta^dag F eta - F| = %.3f'
          % np.linalg.norm(comp - Fp))


if __name__ == '__main__':
    g1_check()
    bell_pair_closure()
    port_idle()
