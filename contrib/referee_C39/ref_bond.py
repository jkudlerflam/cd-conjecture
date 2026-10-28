"""Referee reimplementation of the C39 second-order bond, physical basis.

Built only from the rule table in codex/work_C39/PROOFS.md Section 1-2.
No tree gauge, no eigendecomposition of the excited cell: the bond is
-P V (H0)^-1 V P evaluated by sparse LU on the two-cell excited space.
"""
import itertools, sys, json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def cell_basis(n):
    """Legal single-cell states: (token, location, bits).
    token: 'v' or (r,s), s in (0,1,-1). location: ('A',), ('d',j), ('x',j), ('B',)."""
    N = 2 ** n
    locs_v = [('A',)] + [('d', j) for j in range(n)] + [('B',)]
    locs_e = [('A',)] + [('d', j) for j in range(n)] + [('x', j) for j in range(n)] + [('B',)]
    states = []
    for loc in locs_v:
        for b in range(N):
            states.append(('v', loc, b))
    for r in (1, 2):
        for s in (0, 1, -1):
            for loc in locs_e:
                for b in range(N):
                    states.append(((r, s), loc, b))
    return states, {st: i for i, st in enumerate(states)}


def cell_hamiltonian(n, w, v, D, g, q, a):
    states, ix = cell_basis(n)
    d = len(states)
    H = sp.lil_matrix((d, d), dtype=complex)
    N = 2 ** n
    zbit = lambda b, j: -1 if (b >> j) & 1 else 1
    ph = np.exp(1j * np.pi / (2 * n))
    # vacuum walker: (q|l>-|l+1>)(q<l|-<l+1|)
    vloc = [('A',)] + [('d', j) for j in range(n)] + [('B',)]
    for l in range(n + 1):
        for b in range(N):
            i, j = ix[('v', vloc[l], b)], ix[('v', vloc[l + 1], b)]
            H[i, i] += q * q; H[j, j] += 1; H[i, j] += -q; H[j, i] += -q
    for r in (1, 2):
        wr = w * (v if r == 2 else 1.)
        for s in (0, 1, -1):
            tok = (r, s)
            for b in range(N):
                # onsite: a on data rail, a+D on A,B,x
                for loc in [('A',), ('B',)] + [('x', j) for j in range(n)]:
                    H[ix[(tok, loc, b)], ix[(tok, loc, b)]] += a + D
                for j in range(n):
                    H[ix[(tok, ('d', j), b)], ix[(tok, ('d', j), b)]] += a
                hops = []  # (from, to, amplitude of |to><from|)
                hops.append((('A',), ('d', 0), 1.))
                hops.append((('d', 0), ('x', 0), 1.))
                for j in range(n - 1):
                    hops.append((('d', j), ('d', j + 1), wr * (zbit(b, j) if s == 0 else 1)))
                    hops.append((('x', j), ('x', j + 1), (ph * zbit(b, j)) if s == 0 else 1.))
                hops.append((('x', n - 1), ('B',), (ph * zbit(b, n - 1)) if s == 0 else 1.))
                for fr, to, t in hops:
                    i, j = ix[(tok, fr, b)], ix[(tok, to, b)]
                    H[j, i] += -t; H[i, j] += -np.conj(t)
        # conversions at data sites, bus r:  g|d_+,0><d_0,1| + g|d_-,1><d_0,0| + h.c.
        for j in range(n):
            for b in range(N):
                if (b >> j) & 1:
                    i0 = ix[((r, 0), ('d', j), b)]; ip = ix[((r, 1), ('d', j), b ^ (1 << j))]
                    H[ip, i0] += g; H[i0, ip] += g
                else:
                    i0 = ix[((r, 0), ('d', j), b)]; im = ix[((r, -1), ('d', j), b | (1 << j))]
                    H[im, i0] += g; H[i0, im] += g
    return H.tocsr(), states, ix


def port_maps(n, states, ix, r, end):
    """F_{r,s} at endpoint 'A' or 'B' as sparse matrices on the cell space."""
    d = len(states); N = 2 ** n
    F = {}
    def m(pairs):
        M = sp.lil_matrix((d, d), dtype=complex)
        for i, j in pairs:
            M[i, j] += 1
        return M.tocsr()
    loc = (end,)
    for s in (0, 1, -1):
        pass
    # F_0 = |r,0><v| + |v><r,0| ; F_+ = |r,+><v| + |v><r,-| ; F_- = F_+^dag
    pairs0, pairsp = [], []
    for b in range(N):
        iv = ix[('v', loc, b)]
        pairs0 += [(ix[((r, 0), loc, b)], iv), (iv, ix[((r, 0), loc, b)])]
        pairsp += [(ix[((r, 1), loc, b)], iv), (iv, ix[((r, -1), loc, b)])]
    F[0] = m(pairs0); F[1] = m(pairsp); F[-1] = F[1].getH().tocsr()
    return F


def bond(n, w=.8, v=.65, D=.2, g=.4, q=.8, a=3.5, eta=.6, phis=(np.pi / 2, np.pi), buses=(1, 2)):
    H, states, ix = cell_hamiltonian(n, w, v, D, g, q, a)
    N = 2 ** n; k = n + 2
    S = sum(q ** (2 * l) for l in range(k))
    vloc = [('A',)] + [('d', j) for j in range(n)] + [('B',)]
    G = np.zeros((len(states), N))
    for l in range(k):
        for b in range(N):
            G[ix[('v', vloc[l], b)], b] = q ** l / np.sqrt(S)
    K = np.zeros((N * N, N * N), complex)
    weights = {1: 1., 2: np.sqrt(eta)}
    for r in buses:
        phi = phis[r - 1]
        FB = port_maps(n, states, ix, r, 'B'); FA = port_maps(n, states, ix, r, 'A')
        exc = [i for i, st in enumerate(states) if st[0] != 'v' and st[0][0] == r]
        Pe = sp.csr_matrix((np.ones(len(exc)), (np.arange(len(exc)), exc)), shape=(len(exc), len(states)))
        He = (Pe @ H @ Pe.T).tocsc()
        I = sp.identity(len(exc), format='csc')
        H2 = (sp.kron(He, I) + sp.kron(I, He)).tocsc()
        # V P : code (left B-port cell, right A-port cell) -> excited pair
        Y = np.zeros((len(exc) ** 2, N * N), complex)
        for sB, sA, amp in [(0, 0, 1.), (1, -1, np.exp(1j * phi)), (-1, 1, np.exp(-1j * phi))]:
            L = Pe @ (FB[sB] @ G)   # left cell
            R = Pe @ (FA[sA] @ G)   # right cell
            Y += amp * np.kron(np.asarray(L.todense() if sp.issparse(L) else L), np.asarray(R.todense() if sp.issparse(R) else R))
        lu = sla.splu(H2)
        X = lu.solve(Y)
        K += -weights[r] ** 2 * (Y.conj().T @ X)
    return K


if __name__ == '__main__':
    out = []
    base = sys.argv[1] if len(sys.argv) > 1 else '../../codex/work_C39'
    for n in (1, 2, 3):
        K = bond(n)
        ref = np.load(f'{base}/bond_n{n}.npz')['K2']
        N = 2 ** n
        hop_ids = ([(1 << j) * N for j in range(n)], [1 << j for j in range(n)])
        Th = K[np.ix_(*hop_ids)]; Tr = ref[np.ix_(*hop_ids)]
        row = dict(n=n, max_abs_diff=float(np.abs(K - ref).max()),
                   scalar_mine=float(K[0, 0].real), scalar_ref=float(ref[0, 0].real),
                   hop_rel_err=float(np.abs(Th - Tr).max() / np.abs(Tr).max()),
                   hop_mine=[[complex(z).__repr__() for z in rr] for rr in Th],
                   traceless_rel_err=float(np.linalg.norm((K - ref) - np.trace(K - ref) / (N * N) * np.eye(N * N)) / np.linalg.norm(ref - np.trace(ref) / (N * N) * np.eye(N * N))))
        print(json.dumps(row, indent=1)); out.append(row)
    json.dump(out, open('ref_bond_compare.json', 'w'), indent=1)
