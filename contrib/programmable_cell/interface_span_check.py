"""Span of port-idle interface terms built from ONE port-touching operator (Claude's check).

Port p (dim q), partner s (dim q), rest r (dim m).  The only operator that touches the port
is the Bell projector Pi_Omega on (p, s).  The computation may conjugate by any unitary W
on (s, r) and multiply by any Y on r.  We compute the real dimension of
    span{ (1_p x W)^dag (Pi_Omega x Y) (1_p x W) }
and compare with the dimension of {X Hermitian : tr_{s,r} X proportional to 1_p}, which is
(q m q)^2 - (q^2 - 1).  With no ancilla (m = 1) both partial traces must be scalar, giving
(q q)^2 - 2(q^2 - 1); one ancilla qubit (m >= 2) lifts the partner-side constraint.  Pure-port operators F x 1 can never appear, because partial traces
over (s, r) are invariant under conjugation by W; they are not needed by the gadgets since
the logical factors (gamma_a, T_mu) are traceless.
Run: python3 interface_span_check.py
"""
import numpy as np
rng = np.random.default_rng(3)


def haar(n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    qm, r = np.linalg.qr(a)
    return qm * (np.diag(r) / abs(np.diag(r)))


def herm_basis_vec(H):
    return np.concatenate([H.real.ravel(), H.imag.ravel()])


for q, m in [(2, 1), (2, 2), (3, 1), (2, 4)]:
    D = q * q * m
    om = np.eye(q).reshape(q * q) / np.sqrt(q)
    Pom = np.outer(om, om)
    vecs = []
    for _ in range(4 * D * D):
        W = np.kron(np.eye(q), haar(q * m))                  # acts on (s, r)
        Y = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)); Y = Y + Y.conj().T
        T = W.conj().T @ np.kron(Pom, Y) @ W
        vecs.append(herm_basis_vec(T))
    rank = np.linalg.matrix_rank(np.array(vecs), tol=1e-8)
    pred = D * D - (q * q - 1) * (2 if m == 1 else 1)
    print(f'q={q} m={m}: span dim = {rank}, predicted {pred}, all Hermitian = {D * D}')
