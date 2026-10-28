"""n = 2 logical chains of the four route-A cells, as one 16 x 16 nearest-neighbour bond.

Every bond acts on (cell j) x (cell j+1), each cell C^4 = two logical fermion modes, basis
|n_{j,0} n_{j,1} n_{j+1,0} n_{j+1,1}> with mode 0 the most significant bit and bit 1 = occupied.
Fermions are Jordan-Wigner across the whole chain, c_m = Z...Z sigma^-_m, sigma^- = |0><1|.
The on-site fields are split half and half onto the two bonds of a cell, so that
H = sum_j bond(j, j+1) on an infinite chain or a ring.  Scalars are removed.

The source code of each cell was copied unchanged into src/ (the source folders are read only).
Coupling parameter of each cell, the one that sets the non-Gaussian part relative to the
quadratic part (the overall lambda^2 scale drops out of the second-order chain):
  c28      g   (charge-conversion amplitude), bond in units W = lambda^2 g^2 / 9
  c29      g   (conversion amplitude),        bond in units J = lambda^2 alpha_2^4 g^2
  spatial  g   (gp = gm),                      bond in units of lambda^2
  compile  tau (tauA = tauB),                  bond in units of lambda^2
Each function returns (bond, info).  level='full' is the complete order-lambda^2 bond at that
coupling (all generated terms); level='free' replaces it by its quadratic part (normal ordered
with respect to the empty state, as in the sources' leading chains).
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28.
"""
import os, sys, functools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True


def _path(sub):
    """Put src/<sub> first on the path and drop clashing module names (c29 and compile both
    have a model.py)."""
    p = os.path.join(HERE, 'src', sub)
    for k in ('model', 'model2', 'checks', 'onebody', 'onebody2', 'effective_chain', 'switch_phase'):
        m = sys.modules.get(k)
        if m is not None and not os.path.dirname(os.path.abspath(m.__file__)) == p:
            del sys.modules[k]
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)
    return p


# ------------------------------------------------------------------ fermions (4 modes)
def jw(m):
    Z = np.diag([1.0, -1.0]); I = np.eye(2); sm = np.array([[0, 1.0], [0, 0]])
    out = []
    for j in range(m):
        M = np.ones((1, 1))
        for i in range(m):
            M = np.kron(M, Z if i < j else (sm if i == j else I))
        out.append(M.astype(complex))
    return out


CS = jw(4)
NUM = np.array([bin(i).count('1') for i in range(16)])


def traceless(A):
    return A - np.trace(A) / len(A) * np.eye(len(A))


def site_field(mu):
    """Diagonal one-cell operator mu_0 n_0 + mu_1 n_1 on C^4."""
    c = jw(2)
    return sum(m * (x.conj().T @ x) for m, x in zip(mu, c))


def split_field(F1):
    """One-cell operator F1 (4x4) split half on each adjacent bond."""
    return 0.5 * (np.kron(F1, np.eye(4)) + np.kron(np.eye(4), F1))


def normal_order(K):
    """K = const + sum h_ij c_i^+ c_j + sum_{i<j,k<l} U_{ij,kl} c_i^+ c_j^+ c_l c_k + R (R = 3,4-body).
    Returns const, h (4x4), U (6x6), the quadratic operator, the quartic operator, |R|."""
    vac = np.zeros(16); vac[0] = 1
    one = np.column_stack([c.conj().T @ vac for c in CS])
    const = K[0, 0].real
    h = one.conj().T @ K @ one - const * np.eye(4)
    Hq = const * np.eye(16, dtype=complex)
    for i in range(4):
        for j in range(4):
            Hq += h[i, j] * CS[i].conj().T @ CS[j]
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    two = np.column_stack([CS[i].conj().T @ CS[j].conj().T @ vac for i, j in pairs])
    U = two.conj().T @ (K - Hq) @ two
    H4 = np.zeros_like(K)
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            H4 += U[a, b] * CS[i].conj().T @ CS[j].conj().T @ CS[l] @ CS[k]
    return const, h, U, Hq, H4, float(np.linalg.norm(K - Hq - H4, 2))


def free_part(K):
    c, h, U, Hq, H4, r = normal_order(K)
    return traceless(Hq)


# ------------------------------------------------------------------ C28
@functools.lru_cache(None)
def _c28_fields():
    _path('c28')
    from verify import parameters
    values, K, h = parameters(2)
    return np.array(values)


def c28(g=0.2, level='full'):
    _path('c28')
    from bus_base import h2
    from ff_bus import PHASES
    values = _c28_fields()
    F = np.diag([0, values[1], values[0], values[0] + values[1]]).astype(complex)
    split = (np.kron(F, np.eye(4)) + np.kron(np.eye(4), F)) / 2
    if level == 'trunc':      # leading + g^2 C_4, as saved in logical_chain_n2.npz
        z = np.load(os.path.join(HERE, 'src', 'c28', 'logical_chain_n2.npz'))
        K = z['leading_bond'] + g ** 2 * z['g4_bond']
    else:
        exact, _ = h2(2, g=g, phases=PHASES)
        K = traceless(exact) / g ** 2 + traceless(split)
    K = 0.5 * (K + K.conj().T)
    if level == 'free':
        K = free_part(K)
    return traceless(K), dict(cell='c28', g=g, level=level, units='W = lambda^2 g^2/9', fields=values.tolist())


# ------------------------------------------------------------------ C29
def c29(g=0.2, level='full'):
    _path('c29')
    from model import combined_bond, target_data
    T, mu, W = target_data(2, g)
    if level == 'leading':
        from effective_chain import logical_bond
        K, W = logical_bond(g, 'leading')
    else:
        K = combined_bond(2, g) / g ** 2
    K = traceless(0.5 * (K + K.conj().T)) + split_field(W)
    if level == 'free':
        K = free_part(K)
    return traceless(K), dict(cell='c29', g=g, level=level, units='J = lambda^2 alpha_2^4 g^2', fields=mu.tolist())


# ------------------------------------------------------------------ routeA_spatial
SPATIAL_IR = dict(delta=0.0, Delta=0.0, Dx=4.0, DM=4.0, phd={'0': 0.9, '+': 0.0, '-': 0.0})


@functools.lru_cache(None)
def _spatial_H2(g, td):
    _path('spatial')
    from ycell import params, build_cell, code, second_order
    P = params(td=td, gp=g, gm=g, **SPATIAL_IR)
    cell = build_cell(2, P)
    E0, G, gap = code(cell)
    H2, k = second_order(cell, G, E0)
    return H2, gap


def spatial(g=0.05, level='full', td=1e-3, eps_frac=0.25):
    """Bond lambda^2 H2 + uniform field eps (n_0 + n_1) per cell, eps = eps_frac * 2 * min|eig R_s|
    of the g -> 0 hopping, as in logical_chain.py (a fixed field, independent of g)."""
    H2, gap = _spatial_H2(g, td)
    H2 = 0.5 * (H2 + H2.conj().T)
    H2r, _ = _spatial_H2(0.05, td)
    _, h, _, _, _, _ = normal_order(0.5 * (H2r + H2r.conj().T))
    T = h[:2, 2:]
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    rs = np.linalg.eigvalsh((R + R.conj().T) / 2)
    eps = eps_frac * 2 * np.abs(rs).min() * (g / 0.05) ** 2
    K = traceless(H2) + split_field(site_field([eps, eps]))
    K = 0.5 * (K + K.conj().T)
    if level == 'free':
        K = free_part(K)
    return traceless(K), dict(cell='spatial', g=g, td=td, level=level, units='lambda^2', eps=float(eps), cell_gap=float(gap))


# ------------------------------------------------------------------ routeA_compile (central switch)
COMPILE_BREAK = dict(wD=-0.3, xL=1.05, phiT=1e-3)   # 5% arm asymmetry and a phase on T hops


@functools.lru_cache(None)
def _compile_prep(m):
    _path('compile')
    from model2 import LaneChain, segment_parts
    ch = LaneChain('switch', m, 2, K=2, max_data=None, one_per_cell=True)
    return ch, segment_parts(ch)


def compile_(tau=0.003, level='full', muc=0.0):
    """Full order-lambda^2 bond on the 16 two-cell code states, with the conversion-free scalar
    removed analytically (R - R0 = -R C R0), as in src/compile/onebody2.logical_h_precise."""
    _path('compile')
    import scipy.sparse.linalg as sla
    from model import Chain, DEFAULT
    from model2 import assemble_seg
    ch, P = _compile_prep(2)
    p = dict(DEFAULT); p.update(COMPILE_BREAK); p.update(tauA=tau, tauB=tau)
    H0, V = assemble_seg(P, p)
    H00, _ = assemble_seg(P, dict(p, tauA=0.0, tauB=0.0))
    C = (H0 - H00).tocsr()
    dm = ch.data_modes                      # [cell0 d0, cell0 d1, cell1 d0, cell1 d1]
    assert len(dm) == 4
    vac = [c for c in ch.code if all(((ch.states[c] >> q) & 1) == 0 for q in dm)][0]
    cols, signs = [], []
    for b in range(16):
        bits = [(b >> (3 - i)) & 1 for i in range(4)]
        ops = [(dm[i], 1) for i in range(4) if bits[i]]   # c^+_{m0} c^+_{m1} ... (applied right to left)
        r = Chain.apply(ops, [], ch.states[vac])
        cols.append(ch.index[r[1]]); signs.append(r[0])
    signs = np.array(signs, float)
    mask = np.ones(ch.dim, bool); mask[ch.code] = False
    Q = np.where(mask)[0]
    VQP = V[Q][:, cols].toarray() * signs[None, :]
    X0 = sla.spsolve(H00[Q][:, Q].tocsc(), VQP)
    X1 = sla.spsolve(H0[Q][:, Q].tocsc(), C[Q][:, Q] @ X0)
    D = VQP.conj().T @ X1
    K = traceless(0.5 * (D + D.conj().T)) + split_field(site_field([muc, muc]))
    if level == 'free':
        K = free_part(K)
    return traceless(K), dict(cell='compile', tau=tau, level=level, units='lambda^2', muc=muc, **COMPILE_BREAK)


CELLS = dict(c28=c28, c29=c29, spatial=spatial, compile=compile_)


def bond(cell, x, level='full'):
    return CELLS[cell](x, level)[0]


if __name__ == '__main__':
    for cell, xs in [('c28', [0.2]), ('c29', [0.2]), ('spatial', [0.05]), ('compile', [0.003])]:
        for x in xs:
            K, info = CELLS[cell](x)
            c, h, U, Hq, H4, r = normal_order(K)
            Q = np.diag(NUM)
            print(cell, x, 'norm', np.linalg.norm(K, 2), 'U(1) defect', np.linalg.norm(K @ Q - Q @ K),
                  'herm', np.linalg.norm(K - K.conj().T), 'quartic', np.linalg.norm(U), 'R', r, info)
