"""How the g^4 one-cell and non-Gaussian terms scale with D relative to the g^2 hopping."""
import json, numpy as np
from ref_gauss import series, ferm
rows = []
for n in (1, 2, 3):
    N = 2 ** n; c = ferm(n); P = np.diag([(-1) ** bin(b).count('1') for b in range(N)])
    C = [np.kron(x, np.eye(N)) for x in c] + [np.kron(P, x) for x in c]
    vac = np.zeros(N * N); vac[0] = 1; sp1 = [Ci.T @ vac for Ci in C]
    for D in (2., 8., 32.):
        K0, K2, K4 = series(n, w=.25, v=.5, D=D, q=.5, a=4., eta=.5)
        def parts(K):
            e0 = vac @ K @ vac
            M = np.array([[sp1[a] @ K @ sp1[b] for b in range(2 * n)] for a in range(2 * n)]) - e0 * np.eye(2 * n)
            Kq = e0 * np.eye(N * N, dtype=complex)
            for a in range(2 * n):
                for b in range(2 * n): Kq += M[a, b] * (C[a].T @ C[b])
            return np.linalg.norm(M[:n, n:]), max(np.abs(M[:n, :n]).max(), np.abs(M[n:, n:]).max()), np.linalg.norm(K - Kq, 2)
        h2, o2, ng2 = parts(K2); h4, o4, ng4 = parts(K4)
        rows.append(dict(n=n, D=D, hop2=h2, onsite4_over_hop2=o4 / h2, nonGauss4_over_hop2=ng4 / h2))
        print(rows[-1], flush=True)
json.dump(rows, open('ref_g4scale.json', 'w'), indent=1)
