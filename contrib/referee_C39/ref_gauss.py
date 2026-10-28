"""Referee check: g-expansion of the bond, K = K0 + g^2 K2 + g^4 K4 + ...
K0 is scalar, K2 is an exact fermionic quadratic form whose one-cell part is
scalar, and K4 carries the first non-Gaussian terms (PROOFS Section 2)."""
import json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import ref_bond

def ferm(n):
    N = 2 ** n; out = []
    for j in range(n):
        c = np.zeros((N, N))
        for b in range(N):
            if b >> j & 1: c[b ^ (1 << j), b] = (-1) ** bin(b & ((1 << j) - 1)).count('1')
        out.append(c)
    return out

def series(n, w=.8, v=.65, D=.2, q=.8, a=3.5, eta=.6, phis=(np.pi / 2, np.pi), order=4):
    H0, states, ix = ref_bond.cell_hamiltonian(n, w, v, D, 0., q, a)
    H1, _, _ = ref_bond.cell_hamiltonian(n, w, v, D, 1., q, a)
    Cc = (H1 - H0).tocsr()
    N = 2 ** n; k = n + 2; S = sum(q ** (2 * l) for l in range(k))
    vloc = [('A',)] + [('d', j) for j in range(n)] + [('B',)]
    G = np.zeros((len(states), N))
    for l in range(k):
        for b in range(N): G[ix[('v', vloc[l], b)], b] = q ** l / np.sqrt(S)
    Ks = [np.zeros((N * N, N * N), complex) for _ in range(order // 2 + 1)]
    for r, wt in ((1, 1.), (2, np.sqrt(eta))):
        FB = ref_bond.port_maps(n, states, ix, r, 'B'); FA = ref_bond.port_maps(n, states, ix, r, 'A')
        exc = [i for i, st in enumerate(states) if st[0] != 'v' and st[0][0] == r]
        Pe = sp.csr_matrix((np.ones(len(exc)), (np.arange(len(exc)), exc)), shape=(len(exc), len(states)))
        He = (Pe @ H0 @ Pe.T).tocsc(); Ce = (Pe @ Cc @ Pe.T).tocsc(); I = sp.identity(len(exc), format='csc')
        H2 = (sp.kron(He, I) + sp.kron(I, He)).tocsc(); C2 = (sp.kron(Ce, I) + sp.kron(I, Ce)).tocsr()
        Y = np.zeros((len(exc) ** 2, N * N), complex)
        for sB, sA, amp in [(0, 0, 1.), (1, -1, np.exp(1j * phis[r - 1])), (-1, 1, np.exp(-1j * phis[r - 1]))]:
            Y += amp * np.kron(np.asarray((Pe @ (FB[sB] @ G))), np.asarray((Pe @ (FA[sA] @ G))))
        lu = sla.splu(H2); X = lu.solve(Y)
        cur = X
        for m in range(order + 1):
            if m % 2 == 0:
                Ks[m // 2] += -wt ** 2 * (-1) ** m * (Y.conj().T @ cur)
            cur = lu.solve(C2 @ cur)
    return Ks

if __name__ == '__main__':
  rows = []
  for n in (1, 2, 3):
      N = 2 ** n; c = ferm(n); P = np.diag([(-1) ** bin(b).count('1') for b in range(N)])
      C = [np.kron(x, np.eye(N)) for x in c] + [np.kron(P, x) for x in c]
      vac = np.zeros(N * N); vac[0] = 1; sp1 = [Ci.T @ vac for Ci in C]
      Ks = series(n)
      row = dict(n=n)
      for name, K in zip(['K0', 'K2', 'K4'], Ks):
          e0 = vac @ K @ vac
          M = np.array([[sp1[a] @ K @ sp1[b] for b in range(2 * n)] for a in range(2 * n)]) - e0 * np.eye(2 * n)
          Kq = e0 * np.eye(N * N, dtype=complex)
          for a in range(2 * n):
              for b in range(2 * n): Kq += M[a, b] * (C[a].T @ C[b])
          hop = np.linalg.norm(M[:n, n:])
          row[name] = dict(scalar=float(e0.real), hop=float(hop),
                           onsite_max=float(max(np.abs(M[:n, :n]).max(), np.abs(M[n:, n:]).max())),
                           nonGaussian=float(np.linalg.norm(K - Kq, 2)),
                           traceless=float(np.linalg.norm(K - np.trace(K) / (N * N) * np.eye(N * N), 2)))
      rows.append(row); print(json.dumps(row), flush=True)
  json.dump(rows, open('ref_gauss.json', 'w'), indent=1)
