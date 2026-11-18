"""Exact checks of the all-n sign theorem in MOMENT_SIGN.md.

Reverse Gram-Schmidt is implemented directly, independently of the earlier
eleven-state inverse-minor implementation. Fractions are exact; the proof,
not a finite run of this file, establishes arbitrary n.
"""
from fractions import Fraction as F
from math import factorial


def dot(x, A, y):
    return sum(x[i]*A[i][j]*y[j]
               for i in range(len(x)) for j in range(len(y)))


def matmul(A, B):
    return [[sum(a*b for a, b in zip(row, col)) for col in zip(*B)]
            for row in A]


def transpose(A):
    return list(map(list, zip(*A)))


def congruence(V, A):
    return matmul(matmul(V, A), transpose(V))


def unit_upper_inverse(V):
    n = len(V)
    W = [[F(i == j) for j in range(n)] for i in range(n)]
    for i in range(n-2, -1, -1):
        for j in range(i+1, n):
            W[i][j] = -sum(V[i][k]*W[k][j] for k in range(i+1, j+1))
    return W


def criterion(n):
    C = [[F(factorial(2*i+2*j+2), factorial(2*i+1)*factorial(2*j+1))
          for j in range(n)] for i in range(n)]
    Q = [[C[i][j]*F(i-j, 2*(i+j+1)) for j in range(n)] for i in range(n)]
    J = [[F(0) for _ in range(n)] for _ in range(n)]
    if n >= 2:
        J[-2][-1], J[-1][-2] = F(1), F(-1)
    V = [[F(i == j) for j in range(n)] for i in range(n)]
    norms = [F(0)]*n
    for j in range(n-1, -1, -1):
        for l in range(j+1, n):
            coeff = dot(V[j], C, V[l])/norms[l]
            V[j] = [a-coeff*b for a, b in zip(V[j], V[l])]
        norms[j] = dot(V[j], C, V[j])
    QQ, JJ = congruence(V, Q), congruence(V, J)
    L = [sum(QQ[j][l]*JJ[j][l]/norms[l] for l in range(j+1, n))
         for j in range(n-1)]
    return C, Q, J, V, norms, L


def proof_identities(n):
    C, Q, J, V, norms, L = criterion(n)
    M = [[F(factorial(2*i+2*j+1), factorial(2*i+1)*factorial(2*j+1))
          for j in range(n)] for i in range(n)]
    E = [[F(factorial(2*i+2*j+3), factorial(2*i+1)*factorial(2*j+1))
          for j in range(n)] for i in range(n)]
    A = matmul([[v*(2*i+1) for i, v in enumerate(row)] for row in V],
               unit_upper_inverse(V))
    N, QQ, JJ, EE = [congruence(V, X) for X in (M, Q, J, E)]
    NA = matmul(N, transpose(A))
    AN = matmul(A, N)
    for j in range(n):
        assert norms[j] > 0 and A[j][j] == 2*j+1
        for l in range(n):
            assert dot(V[j], C, V[l]) == (norms[j] if j == l else 0)
            assert AN[j][l]+NA[j][l] == (norms[j] if j == l else 0)
            assert QQ[j][l] == (AN[j][l]-NA[j][l])/2
            assert A[j][l]*norms[l]+A[l][j]*norms[j] == EE[j][l]-(norms[j] if j == l else 0)
        for l in range(j+1, n):
            assert (-1)**(j+l)*A[j][l] > 0
            assert (-1)**(j+l+1)*N[j][l] > 0
            assert (-1)**(j+l)*QQ[j][l] > 0
            assert (-1)**(j+l+1)*JJ[j][l] > 0
            assert QQ[j][l]*JJ[j][l] < 0
    assert all(x < 0 for x in L)
    if n >= 2:
        assert L[-1] == -F(1, 16*n-12)
    return dict(n=n, L=list(map(str, L)), exact=True,
                all_summands_strictly_negative=True,
                identities=['orthogonality', 'AN+NA^T=diag(d)',
                            'Q=(AN-NA^T)/2', 'dilation integration by parts'])


if __name__ == '__main__':
    for n in range(2, 7):
        print(proof_identities(n))
