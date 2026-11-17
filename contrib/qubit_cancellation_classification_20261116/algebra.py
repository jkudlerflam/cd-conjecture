"""Full-sector spin matrices and nested Schur complements. No CFT inference."""
import numpy as np
from scipy import sparse


def occupations(L):
    return np.array([[(s >> j) & 1 for j in range(L)] for s in range(1 << L)])


def bond(v=0., u=0., t=0j, p=0j, f=0j, a=0j, b=0j):
    """t moves right; a/b create with the left/right neighbour occupied."""
    I = np.eye(2); n = np.diag([0., 1.])
    bd = np.array([[0., 0.], [1., 0.]]); low = bd.T
    h = v*np.kron(n, n) + u/2*(np.kron(n, I)+np.kron(I, n))
    z = (t*np.kron(low, bd) + p*np.kron(bd, bd)
         + f/2*(np.kron(bd, I)+np.kron(I, bd))
         + a*np.kron(n, bd) + b*np.kron(bd, n))
    return h+z+z.conj().T


def chain(L, h):
    """Literal 4x4 bond embedding, independent of effective local rules."""
    rows, cols, vals = [], [], []
    for s in range(1 << L):
        for j in range(L):
            k = (j+1) % L
            col = 2*((s >> j) & 1)+((s >> k) & 1)
            base = s & ~((1 << j) | (1 << k))
            for row in range(4):
                if h[row, col] != 0:
                    rows.append(base | ((row >> 1) << j) | ((row & 1) << k))
                    cols.append(s); vals.append(h[row, col])
    return sparse.coo_matrix((vals, (rows, cols)), shape=(1 << L,)*2).tocsr()


def sub(A, rows, cols=None):
    return A[rows][:, rows if cols is None else cols]


def maxabs(A):
    if sparse.issparse(A):
        return float(np.max(np.abs(A.data), initial=0.))
    return float(np.max(np.abs(np.asarray(A)), initial=0.))


def projectors(L, blockade=1):
    ns = occupations(L)
    energies = sum(np.sum(ns*np.roll(ns, -r, axis=1), axis=1)
                   for r in range(1, blockade+1))
    return np.flatnonzero(energies == 0)


def second_local(L, a, b, p, u=0., t=0j, f=0j):
    """Closed-form second order, including all number-changing matrix entries."""
    P = projectors(L); pos = {s:i for i, s in enumerate(P)}
    ns = occupations(L); pp = abs(p)**2
    out = np.zeros((len(P),)*2, complex)
    for s in P:
        i = pos[s]; n = ns[s]; N = sum(n)
        D2 = sum(n*np.roll(n, -2)); D3 = sum(n*np.roll(n, -3))
        out[i, i] = (-L*pp + (u+3*pp-abs(a)**2-abs(b)**2)*N
                     + (abs(a-b)**2/2-pp)*D2 - pp/3*D3)
        for j in range(L):
            for r, amp in [(1, t-a*np.conj(b)), (2, -pp/2)]:
                k = (j+r) % L
                if n[j] == 1 and n[k] == 0:
                    dest = s ^ (1 << j) ^ (1 << k)
                    if dest in pos:
                        out[pos[dest], i] += amp
                        out[i, pos[dest]] += np.conj(amp)
            if n[j] == 0:
                dest = s | (1 << j)
                if dest in pos:
                    amp = (f-p*np.conj(a+b)
                           + p*np.conj(a-b)/2*(n[(j+2)%L]-n[(j-2)%L]))
                    out[pos[dest], i] += amp
                    out[i, pos[dest]] += np.conj(amp)
    return out


def nested_fourth(L, q, dense=False, source3=0j, pair3=0j):
    """Exact coefficient arithmetic in doubles, not an epsilon fit.

    H=D1+eps F+eps^2(2N-T+q Pair+h.c.). The first Schur complement
    removes Q=1-P. The second removes P-P2 with denominator 2D2.
    Its P2 fourth coefficient agrees with the Hermitian effective one:
    all metric terms have zero compression because K2 P2=0.
    """
    ns = occupations(L)
    E = np.sum(ns*np.roll(ns, -1, axis=1), axis=1)
    D2 = np.sum(ns*np.roll(ns, -2, axis=1), axis=1)
    P = np.flatnonzero(E == 0); Q = np.flatnonzero(E > 0)
    constraint = L-3*np.sum(ns, axis=1)+D2 if dense else D2
    I = np.flatnonzero(constraint[P] == 0); J = np.flatnonzero(constraint[P] > 0)
    F = chain(L, bond(a=1, b=-1))
    W = chain(L, bond(u=-4 if dense else 2, t=-1, p=q))
    if dense:
        W = W+2*L*sparse.eye(1 << L, dtype=complex, format='csr')
    V3 = chain(L, bond(f=source3, p=pair3))
    R = sparse.diags(1/E[Q])
    A = sub(F, Q, P); B = sub(W, Q, P)
    C = sub(F, Q); D = sub(W, Q)
    RA = R@A; RB = R@B
    S2 = sub(W, P)-A.conj().T@RA
    S3 = sub(V3, P)+RA.conj().T@C@RA-A.conj().T@RB-B.conj().T@RA
    CRA = C@RA
    S4 = (-B.conj().T@RB + CRA.conj().T@RB + RB.conj().T@CRA
          + RA.conj().T@D@RA - CRA.conj().T@R@CRA
          - RA.conj().T@sub(V3, Q, P)-sub(V3, P, Q)@RA)
    R2 = sparse.diags(1/(2*constraint[P[J]]))
    K4 = sub(S4, I)-sub(S3, I, J)@R2@sub(S3, J, I)
    return P[I], K4.toarray(), {
        'second_order_error': maxabs(S2-sparse.diags(2.*constraint[P])),
        'third_order_P2_norm': maxabs(sub(S3, I)),
        'fourth_hermiticity_error': maxabs(K4-K4.conj().T),
    }


def schur_series(coefficients, keep):
    """Static Schur series for diagonal H[0], kernel 'keep'.

    Exact formal finite-order recursion, evaluated in complex doubles.
    This yields the leading effective coefficient when all lower kept
    coefficients vanish. Higher coefficients need the energy-dependent metric.
    """
    dim = coefficients[0].shape[0]; order = len(coefficients)-1
    remove = np.setdiff1d(np.arange(dim), keep)
    R = sparse.diags(1/coefficients[0].diagonal()[remove])
    B = [sub(c, remove, keep) for c in coefficients]
    C = [sub(c, remove) for c in coefficients]
    X = [sparse.csr_matrix((len(remove), len(keep)), dtype=complex)]
    result = [sub(coefficients[0], keep)]
    for n in range(1, order+1):
        rhs = B[n].copy()
        for k in range(1, n):
            rhs -= C[k]@X[n-k]
        X.append(R@rhs)
        s = sub(coefficients[n], keep).copy()
        for k in range(1, n):
            s -= B[k].conj().T@X[n-k]
        result.append(s)
    return result


def dense_coefficients(q, f):
    qq = abs(q)**2; ff = abs(f)**2; x=q.real
    C = .5-x+(qq+ff)/6
    B = -3.5+3*x-5*qq/6-2*ff/3
    A = abs(q-1)**2/2-ff/6
    t = (abs(q-1)**2-ff)/4 + .5j*np.imag((3+np.conj(q))*f)
    return C, B, A, t


def dense_higher(L, q, f, order=6):
    """Fourth-scalar dense locus, with all unlisted microscopic jets zero."""
    C, B, A, t = dense_coefficients(q, f)
    assert abs(A) < 1e-10
    dim=1 << L; eye=sparse.eye(dim, dtype=complex, format='csr')
    H = [chain(L,bond(v=1)), chain(L,bond(a=1,b=-1)),
         chain(L,bond(u=-4,t=-1,p=q))+2*L*eye,
         chain(L,bond(f=f)), chain(L,bond(u=-B,t=-t))-C*L*eye]
    H += [sparse.csr_matrix((dim,dim), dtype=complex) for _ in range(order-4)]
    P=projectors(L); S=schur_series(H,P)
    keep=np.flatnonzero(abs(S[2].diagonal())<1e-10)
    K=schur_series(S[2:],keep)
    return P[keep], K, H


def dilute_higher(L, q, order=6):
    """Three-stage elimination on the tuned dilute path, all unused jets zero."""
    dim=1 << L;eye=sparse.eye(dim,dtype=complex,format='csr')
    H=[chain(L,bond(v=1)),chain(L,bond(a=1,b=-1)),
       chain(L,bond(u=2,t=-1,p=q)),sparse.csr_matrix((dim,dim),dtype=complex),
       chain(L,bond(u=-(4+4*q.real+2*abs(q)**2),t=2+2*q.real))+L*abs(q)**2*eye]
    H += [sparse.csr_matrix((dim,dim),dtype=complex) for _ in range(order-4)]
    P=projectors(L);S=schur_series(H,P)
    keep2=np.flatnonzero(abs(S[2].diagonal())<1e-10)
    K=schur_series(S[2:],keep2)
    assert maxabs(K[2]-sparse.diags(K[2].diagonal()))<1e-9
    keep4=np.flatnonzero(abs(K[2].diagonal())<1e-9)
    final=schur_series(K[2:],keep4)
    return P[keep2][keep4],final,H
