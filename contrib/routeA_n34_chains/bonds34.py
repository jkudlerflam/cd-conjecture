"""Logical bonds of the route-A cells at n = 3, 4 flavours (also any n), as 4^n x 4^n matrices.

Basis of a bond: |n_{j,0} .. n_{j,n-1} n_{j+1,0} .. n_{j+1,n-1}>, mode 0 of cell j the most
significant bit, bit 1 = occupied, fermions Jordan-Wigner ordered cell j then cell j+1.  Every
bond conserves the logical fermion number, so the spin chain sum_j K_{j,j+1} is the fermion chain.

compile (contrib/routeA_compile, central-switch cell, layout P L^m M R^m P, m = n):
    K2 = exact order-lambda^2 bond (traceless), computed charge sector by charge sector with the
         conversion-free scalar removed analytically (R - R0 = -R C R0), so no order-one scalar
         is cancelled numerically;
    K3 = exact order-lambda^3 bond (traceless) = P V R V R V P.  First order vanishes, so this is
         the whole third-order Schrieffer-Wolff term; V^3 on one bond is the only contribution
         (a cell touched once by V keeps a carrier), and it stays in the one-carrier-per-cell
         sector, which is the basis used here.
spatial (contrib/routeA_spatial, 23-letter folded bus):
    K2 = exact order-lambda^2 bond from ycell.second_order (no odd orders exist).

Source files are copied unchanged into src/ (the source folders are not modified).
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28.
"""
import os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True


def _path(sub):
    p = os.path.join(HERE, 'src', sub)
    if p not in sys.path:
        sys.path.insert(0, p)


def traceless(A):
    return A - np.trace(A) / len(A) * np.eye(len(A))


# ---------------------------------------------------------------- compile (central switch)
COMPILE_BREAK = dict(wD=-0.3, xL=1.05, phiT=1e-3)   # as in routeA_compile/switch_break.py (n = 2)


def compile_prepare(m):
    _path('compile')
    from model2 import LaneChain, segment_parts
    t = time.time()
    ch = LaneChain('switch', m, 2, K=2, max_data=None, one_per_cell=True)
    P = segment_parts(ch)
    print('compile prepare m=%d dim=%d (%.0fs)' % (m, ch.dim, time.time() - t), flush=True)
    return ch, P


def compile_bonds(ch, P, tau, third=True, extra=None):
    """Return (K2, K3, info).  K2, K3 traceless, Hermitian, 4^m x 4^m."""
    _path('compile')
    from model import Chain, DEFAULT
    from model2 import assemble_seg
    m = ch.n
    p = dict(DEFAULT); p.update(COMPILE_BREAK); p.update(tauA=tau, tauB=tau)
    if extra:
        p.update(extra)
    H0, V = assemble_seg(P, p)
    H00, _ = assemble_seg(P, dict(p, tauA=0.0, tauB=0.0))
    C = (H0 - H00).tocsr()
    dm = ch.data_modes
    nm = len(dm)
    assert nm == 2 * m
    vac = [c for c in ch.code if all(((ch.states[c] >> q) & 1) == 0 for q in dm)][0]
    N = 2 ** nm
    cols, signs = [], []
    for b in range(N):
        bits = [(b >> (nm - 1 - i)) & 1 for i in range(nm)]
        ops = [(dm[i], 1) for i in range(nm) if bits[i]]
        r = Chain.apply(ops, [], ch.states[vac])
        cols.append(ch.index[r[1]]); signs.append(r[0])
    cols = np.array(cols); signs = np.array(signs, float)
    ncode = np.array([bin(b).count('1') for b in range(N)])
    num = ch.number
    mask = np.ones(ch.dim, bool); mask[ch.code] = False
    K2 = np.zeros((N, N), complex)
    K3 = np.zeros((N, N), complex)
    Vc = V.tocsc()
    for q in sorted(set(num[ch.code])):
        sel = np.where(num[cols] == q)[0]           # code columns with this physical number
        if len(sel) == 0:
            continue
        Q = np.where(mask & (num == q))[0]
        VQP = Vc[Q][:, cols[sel]].toarray() * signs[sel][None, :]
        A0 = sla.splu(H00[Q][:, Q].tocsc())
        A1 = sla.splu(H0[Q][:, Q].tocsc())
        Y0 = A0.solve(VQP)                          # H00^-1 V P
        Dl = -A1.solve(C[Q][:, Q] @ Y0)             # H0^-1 - H00^-1 applied to V P
        K2s = -VQP.conj().T @ Dl                    # H2 - H2(tau=0); the latter is a scalar
        K2[np.ix_(sel, sel)] = K2s
        if third:
            VQQ = V[Q][:, Q]
            Y = Y0 + Dl
            K3s = Dl.conj().T @ (VQQ @ Y) + Y0.conj().T @ (VQQ @ Dl)   # H3 - H3(tau=0)
            K3[np.ix_(sel, sel)] = K3s
        del A0, A1
    K2 = traceless(0.5 * (K2 + K2.conj().T))
    K3 = traceless(0.5 * (K3 + K3.conj().T))
    return K2, K3, dict(cell='compile', m=m, tau=tau, **{k: p[k] for k in ('wD', 'xL', 'phiT')})


# ---------------------------------------------------------------- spatial (folded bus)
SPATIAL_IR = dict(delta=0.0, Delta=0.0, Dx=4.0, DM=4.0, phd={'0': 0.9, '+': 0.0, '-': 0.0})


def spatial_bond(n, g, td, Dx=4.0):
    _path('spatial')
    from ycell import params, build_cell, code, second_order
    P = params(td=td, gp=g, gm=g, **dict(SPATIAL_IR, Dx=Dx, DM=Dx))
    cell = build_cell(n, P)
    E0, G, gap = code(cell)
    H2, k = second_order(cell, G, E0)
    K2 = traceless(0.5 * (H2 + H2.conj().T))
    return K2, dict(cell='spatial', n=n, g=g, td=td, Dx=Dx, cell_gap=float(gap), poles=int(k))


if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'compile':
        m = int(sys.argv[2]); taus = [float(x) for x in sys.argv[3:]]
        ch, P = compile_prepare(m)
        for tau in taus:
            t = time.time()
            K2, K3, info = compile_bonds(ch, P, tau)
            fn = os.path.join(HERE, 'data', 'bond_compile_n%d_tau%g.npz' % (m, tau))
            np.savez_compressed(fn, K2=K2, K3=K3, **{k: v for k, v in info.items()})
            print(fn, 'norm K2 %.3e K3 %.3e (%.0fs)' % (np.linalg.norm(K2, 2), np.linalg.norm(K3, 2),
                                                      time.time() - t), flush=True)
    elif which == 'spatial':
        n = int(sys.argv[2]); td = float(sys.argv[3]); Dx = float(sys.argv[4])
        for g in [float(x) for x in sys.argv[5:]]:
            t = time.time()
            K2, info = spatial_bond(n, g, td, Dx)
            fn = os.path.join(HERE, 'data', 'bond_spatial_n%d_g%g_td%g_Dx%g.npz' % (n, g, td, Dx))
            np.savez_compressed(fn, K2=K2, **{k: v for k, v in info.items()})
            print(fn, 'norm K2 %.3e (%.0fs)' % (np.linalg.norm(K2, 2), time.time() - t), info, flush=True)
