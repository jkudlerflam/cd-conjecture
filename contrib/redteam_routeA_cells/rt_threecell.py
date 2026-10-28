"""Second order on an open chain of three cells (n = 1, 2): P V R V P with
V = V_12 + V_23 against the sum of the two-cell bonds.  A nonzero difference
would be a three-cell (next-nearest-neighbour) term at order lambda^2."""
import sys, json
import numpy as np
import scipy.sparse as sp
from rt_bond import P
from rt_cell import cell, code_basis, port_factors, port_V, second_order
out = []
for n in [int(x) for x in sys.argv[1:]] or [1]:
    p = P(td=1e-3)
    Hd, info = cell(n, p)
    cb = code_basis(n, Hd, info)
    G, E0 = cb['G'], cb['E0']
    w, U = np.linalg.eigh(Hd)
    D = len(w)
    FM, FX = port_factors(n, info)
    Vb = sum(c * sp.kron(sp.csr_matrix(A), sp.csr_matrix(B)) for (A, B, c) in port_V(p, FM, FX))
    I = sp.identity(D)
    V = (sp.kron(Vb, I) + sp.kron(I, Vb)).tocsr()
    N = G.shape[1]
    G3 = np.einsum('ia,jb,kc->ijkabc', G, G, G).reshape(D ** 3, N ** 3)
    X = V @ G3
    Xe = np.einsum('ai,bj,ck,abcK->ijkK', U.conj(), U.conj(), U.conj(), X.reshape(D, D, D, -1), optimize=True)
    den = w[:, None, None] + w[None, :, None] + w[None, None, :] - 3 * E0
    mask = np.abs(den) < 1e-9
    Xcode = np.abs(Xe[mask]).max()
    Ye = np.where(mask[..., None], 0.0, Xe / np.where(mask, 1.0, den)[..., None])
    H2_3 = -np.einsum('ijkK,ijkL->KL', Xe.conj(), Ye)
    H2 = second_order(n, p)['H2']
    pred = np.kron(H2, np.eye(N)) + np.kron(np.eye(N), H2)
    r = dict(n=n, dim=D ** 3, V_on_code_component=float(Xcode), max_diff=float(np.abs(H2_3 - pred).max()),
             H2_scale=float(np.abs(H2).max()))
    print(r, flush=True); out.append(r)
json.dump(out, open('threecell.json', 'w'), indent=1)
