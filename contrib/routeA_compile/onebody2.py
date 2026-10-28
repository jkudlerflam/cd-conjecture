"""Exact order-lambda^2 one-body logical bond for the lane layouts (model2)."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
import scipy.sparse.linalg as sla
from model import Chain, assemble, DEFAULT
from model2 import LaneChain

_C = {}


def prepared(layout, m):
    key = (layout, m)
    if key not in _C:
        ch = LaneChain(layout, m, 2, K=2, max_data=3, one_per_cell=True)
        _C[key] = (ch, ch.lane_parts())
    return _C[key]


def logical_h(layout, m, params=None):
    ch, P = prepared(layout, m)
    p = dict(DEFAULT)
    if params:
        p.update(params)
    H0, V = assemble(P, p)
    dm = ch.data_modes
    vac = [c for c in ch.code if all(((ch.states[c] >> q) & 1) == 0 for q in dm)][0]
    cols, signs = [vac], [1.0]
    for q in dm:
        sg, new = Chain.apply([(q, 1)], [], ch.states[vac])
        cols.append(ch.index[new]); signs.append(sg)
    signs = np.array(signs)
    mask = np.ones(ch.dim, bool); mask[ch.code] = False
    Q = np.where(mask)[0]
    VQP = V[Q][:, cols].toarray() * signs[None, :]
    HQQ = H0[Q][:, Q].tocsc()
    gap = float(sla.eigsh(HQQ, k=1, which='SA', return_eigenvectors=False)[0])
    X = sla.spsolve(HQQ, VQP)
    H2 = -VQP.conj().T @ X
    H2 = 0.5 * (H2 + H2.conj().T)
    E0 = H2[0, 0].real
    h = H2[1:, 1:] - E0 * np.eye(len(dm))
    return E0, h, gap


def logical_h_precise(layout, m, params=None):
    """Same as logical_h, but the conversion-free scalar is removed analytically:
    with H0 = H00 + C (C = conversion terms), R - R0 = -R C R0, so
    H2 - H2(tau=0) = P V R C R0 V P is evaluated as a product (no cancellation
    of the order-one scalar).  H2(tau=0) is proportional to the identity on the
    code (carriers never touch data without conversions)."""
    ch, P = prepared(layout, m)
    p = dict(DEFAULT)
    if params:
        p.update(params)
    H0, V = assemble(P, p)
    p0 = dict(p, tauA=0.0, tauB=0.0)
    H00, _ = assemble(P, p0)
    C = (H0 - H00).tocsr()
    dm = ch.data_modes
    vac = [c for c in ch.code if all(((ch.states[c] >> q) & 1) == 0 for q in dm)][0]
    cols, signs = [vac], [1.0]
    for q in dm:
        sg, new = Chain.apply([(q, 1)], [], ch.states[vac])
        cols.append(ch.index[new]); signs.append(sg)
    signs = np.array(signs)
    mask = np.ones(ch.dim, bool); mask[ch.code] = False
    Q = np.where(mask)[0]
    VQP = V[Q][:, cols].toarray() * signs[None, :]
    X0 = sla.spsolve(H00[Q][:, Q].tocsc(), VQP)
    Y = C[Q][:, Q] @ X0
    X1 = sla.spsolve(H0[Q][:, Q].tocsc(), Y)
    D = VQP.conj().T @ X1
    D = 0.5 * (D + D.conj().T)
    h = D[1:, 1:] - D[0, 0].real * np.eye(len(dm))
    return D[0, 0].real, h
