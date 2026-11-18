"""Markerless 2n+2-vertex carrier and finite-delay response.

The Sylvester method follows eleven_state_loop_splitting_20261118/loop.py,
with graph endpoints and physical clock independently checked here.
No marker: n+1 express vertices, eta=t_d D^n, scale D^(2n+2).
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import splu


def carrier(n, D, w, a=4., td=0.):
    p = n+1
    s = 2*p
    H = np.diag([a]*p+[a+D]*p).astype(complex if td else float)
    for j in range(n):
        H[j, j+1] = H[j+1, j] = -w
        H[p+j, p+j+1] = H[p+j+1, p+j] = -1.
    H[0, p] = H[p, 0] = -1.
    J = np.zeros((s, s))
    J[-1, n], J[n, -1] = -1., 1.
    if td:
        H += 1j*td*J
    return H, J


def restrictions(x, n, s):
    U = np.zeros((s*s, n), dtype=x.dtype)
    V = U.copy()
    x = x.reshape(s, s)
    for j in range(n):
        U[j*s:(j+1)*s, j] = x[j, :]
        V[j::s, j] = x[:, j]
    return U, V


def coefficients(n, D, w, a=4., td=0.):
    H, J = carrier(n, D, w, a, td)
    s = len(H)
    lu = splu(sp.kronsum(sp.csc_matrix(H), sp.csc_matrix(H), format='csc'))
    b = np.zeros(s*s)
    b[(s-1)*s] = 1.
    x = lu.solve(b)
    U, V = restrictions(x, n, s)
    RU, RV = lu.solve(U), lu.solve(V)
    T = -2*(U.conj().T@RV).real
    M = -2j*(U.conj().T@RU+V.conj().T@RV).imag
    if td:
        return T, M
    J2 = sp.kronsum(sp.csc_matrix(J), sp.csc_matrix(J), format='csc')
    y = -lu.solve(J2@x)
    U1, V1 = restrictions(y, n, s)
    B = U.T@lu.solve(U1)-U1.T@RU+V.T@lu.solve(V1)-V1.T@RV
    B -= RU.T@(J2@RU)+RV.T@(J2@RV)
    assert la.norm(B+B.T) < 1e-10*max(la.norm(B), 1e-100)
    return T, -2j*(B-B.T)/2


def low_limit(n, w, a=4.):
    p = n+1
    A = a*np.eye(p)-w*(np.eye(p, k=1)+np.eye(p, k=-1))
    R = la.inv(np.kron(A, np.eye(p))+np.kron(np.eye(p), A))
    x0, xn = R[:, 0], R[:, n*p]
    U, V = restrictions(x0, n, p)
    Un, Vn = restrictions(xn, n, p)
    G = 2*U.T@R@V
    B = U.T@R@Un-Un.T@R@U+V.T@R@Vn-Vn.T@R@V
    y = R@np.kron(A, np.eye(p))[:, 0]
    Yu, Yv = restrictions(y, n, p)
    Q = -Yu.T@R@V-U.T@R@Yv
    Q[0, :] += x0@V
    assert la.norm(Q+Q.T) < 1e-10*la.norm(G)
    return G, B, Q


def roots_and_response(T, M1, F):
    n = len(T)
    C = la.inv(la.cholesky(-(T+T.T)/2, lower=True))
    T, M1, F = [C@X@C.T for X in (T, M1, F)]
    I, Z = np.eye(n), np.zeros((n, n))
    zs = la.eigvals(np.block([[-F, -T.T], [I, Z]]), np.block([[T, Z], [Z, I]]))
    assert np.max(abs(abs(zs)-1)) < 1e-7
    result = []
    for k in sorted(float(np.angle(z)) for z in zs if np.angle(z) > 0):
        e, U = la.eigh(T*np.exp(1j*k)+T.T*np.exp(-1j*k)+F)
        z = U[:, np.argmin(abs(e))]
        speed = np.vdot(z, (1j*T*np.exp(1j*k)-1j*T.T*np.exp(-1j*k))@z).real
        result.append((k, -2*np.vdot(z, M1@z).real/speed))
    return result
