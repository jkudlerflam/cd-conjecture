"""C166 library (Claude Opus 5.5 agent, 2026-11-08).

Bonds of the XZ neighbourhood and the C163 path passes, a DMRG wrapper with arbitrary edge fields, and
observables on the open-chain MPS returned by dmrg.dmrg (left-canonical, orthogonality centre on the last site).

Conventions (C155/C163, est.py): v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz) normalized to |v| = 1,
J = est.J_f(v), H = sum_j [ s_j^T J s_{j+1} + f . s_j ], Pauli matrices, left factor = left site.
Open chains carry f on every site (f/2 from each adjacent bond plus f/2 at the two ends, dmrg_c.chain_bond).
VUMPS (fesx.mpo_cell) takes Js = [J]*n and Bs = [f]*n, one copy of f per site.
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import est as E
import dmrg as D
import dmrg_c as DC

PA = E.PAULI
I2 = np.eye(2, dtype=complex)

# bond vectors (copied from the JSON files named alongside, so the job does not depend on other folders)
VEC = {
    # codex/work_C163/typhon_out/big/hi_T_40_256.json ['v']
    'hi_T_40': [-0.4026997696385143, -0.4234362659483195, 0.6523832103833125, -0.4304897319791202,
                -0.024769178208296865, -0.21678526124768083, 0.0, 0.0, 0.0],
    # codex/work_C163/typhon_out/big/dist_T_50_256.json ['v']
    'dist_T_50': [-0.5810455731239875, 0.05369144058159372, 0.6336508921426355, 0.008284265794989032,
                  -0.4768573651647684, -0.1747233306587447, 0.0, 0.0, 0.0],
}
PATH_LAM = {'KI22': 0.4355555555555556, 'KI26': 0.2577777777777778, 'N34': 0.3022222222222222}  # C163 summary.txt


def load_vec(name):
    if name in VEC:
        v = np.array(VEC[name], float)
    elif name.startswith('cand_'):
        fn = os.path.join(HERE, 'inputs', name + '.json')
        v = np.array(json.load(open(fn))['v'], float)
    elif name in PATH_LAM:
        import c163lib as C
        spec = [s for s in C.make_specs() if s['tag'] == name][0]
        J, f = C.path_bond(spec, PATH_LAM[name])
        v = C.J_to_v(J, f)
    else:
        raise KeyError(name)
    return v / np.linalg.norm(v)


def JB(name):
    J, f = E.J_f(load_vec(name))
    return J, f


# ------------------------------------------------------------------ exact diagonalization (open and periodic)
def site_op(O, j, L):
    import scipy.sparse as sp
    return sp.kron(sp.kron(sp.identity(2 ** j, format='csr'), sp.csr_matrix(O)), sp.identity(2 ** (L - j - 1), format='csr'), format='csr')


def H_sparse(J, f, L, periodic=False):
    import scipy.sparse as sp
    H = sp.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    S = [[site_op(PA[a], j, L) for a in range(3)] for j in range(L)]
    nb = L if periodic else L - 1
    for j in range(nb):
        k = (j + 1) % L
        for a in range(3):
            for b in range(3):
                if J[a, b] != 0:
                    H = H + J[a, b] * (S[j][a] @ S[k][b])
    for j in range(L):
        for a in range(3):
            if f[a] != 0:
                H = H + f[a] * S[j][a]
    return H


def ed_e0(J, f, L, periodic=False, k=1):
    import scipy.sparse.linalg as spla
    H = H_sparse(J, f, L, periodic)
    w = spla.eigsh(H, k=k, which='SA', tol=1e-12)[0]
    return np.sort(w.real)


# ------------------------------------------------------------------ DMRG with arbitrary edge fields
def dmrg_edges(name, L, chis, left=None, right=None, seed=0, cutoff=1e-10, verbose=False):
    """Open chain of bond `name`; extra on-site fields left (3-vector, on site 0) and right (on site L-1),
    on top of the uniform field of the bond.  Implemented through dmrg.dmrg's edge_field slot by
    monkey-patching bond_mpo (dmrg.py itself is unchanged)."""
    v = load_vec(name)
    h, edge = DC.chain_bond(v)
    eL = edge + (0 if left is None else sum(left[a] * PA[a] for a in range(3)))
    eR = edge + (0 if right is None else sum(right[a] * PA[a] for a in range(3)))
    orig = D.bond_mpo

    def patched(hh, LL):
        Ws = orig(hh, LL)
        Ws[0][0, -1] = Ws[0][0, -1] + eL
        Ws[LL - 1][0, 0] = Ws[LL - 1][0, 0] + eR
        return Ws
    D.bond_mpo = patched
    try:
        r = D.dmrg(h, L, chis=chis, seed=seed, edge_field=None, cutoff=cutoff, verbose=verbose)
    finally:
        D.bond_mpo = orig
    return r, h, v


# ------------------------------------------------------------------ observables on a left-canonical open MPS
def _right_ids(M):
    L = len(M)
    R = [None] * (L + 1)
    R[L] = np.ones((1, 1), complex)
    for j in range(L - 1, -1, -1):
        A = M[j]
        R[j] = np.einsum('asc,cd,bsd->ab', np.conj(A), R[j + 1], A, optimize=True)
    return R


def onesite(M, ops, R=None):
    """<O_j> for all j; returns (len(ops), L)."""
    L = len(M)
    R = R or _right_ids(M)
    out = np.zeros((len(ops), L))
    Lf = np.ones((1, 1), complex)
    for j in range(L):
        A = M[j]
        for k, O in enumerate(ops):
            T = np.einsum('ab,asc,st,btd->cd', Lf, np.conj(A), O, A, optimize=True)
            out[k, j] = np.einsum('cd,cd->', T, R[j + 1], optimize=True).real
        Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A, optimize=True)
    return out


def bond_energies(M, h):
    """<h_{j,j+1}> (the 4x4 chain bond, including its f/2 + f/2 field part) for j = 0..L-2."""
    L = len(M)
    R = _right_ids(M)
    h4 = h.reshape(2, 2, 2, 2)  # (s1 s2, t1 t2)
    out = np.zeros(L - 1)
    Lf = np.ones((1, 1), complex)
    for j in range(L - 1):
        A, B = M[j], M[j + 1]
        K = np.einsum('asc,ctd->astd', A, B, optimize=True)                 # ket (a, s1, s2, d)
        Hk = np.einsum('uvst,astd->auvd', h4, K, optimize=True)
        val = np.einsum('ab,auvd,buve,de->', Lf, np.conj(K), Hk, R[j + 2], optimize=True)
        out[j] = val.real
        Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A, optimize=True)
    return out


def correlator_pairs(M, i, js, ops1, ops2):
    """<A_i B_j> for A in ops1, B in ops2, j in js (all j > i). Returns array (len(ops1), len(ops2), len(js))."""
    L = len(M)
    R = _right_ids(M)
    Lf = np.ones((1, 1), complex)
    for k in range(i):
        A = M[k]
        Lf = np.einsum('ab,asc,bsd->cd', Lf, np.conj(A), A, optimize=True)
    A = M[i]
    T = [np.einsum('ab,asc,st,btd->cd', Lf, np.conj(A), O, A, optimize=True) for O in ops1]
    js = sorted(js)
    out = np.zeros((len(ops1), len(ops2), len(js)), complex)
    cur = i + 1
    for n, j in enumerate(js):
        while cur < j:
            A = M[cur]
            T = [np.einsum('ab,asc,bsd->cd', t, np.conj(A), A, optimize=True) for t in T]
            cur += 1
        A = M[j]
        for p, t in enumerate(T):
            for q, O in enumerate(ops2):
                Y = np.einsum('ab,asc,st,btd->cd', t, np.conj(A), O, A, optimize=True)
                out[p, q, n] = np.einsum('cd,cd->', Y, R[j + 1], optimize=True)
    return out.real


def measure(M, h, L):
    """Profiles and the long-distance correlation matrix used by C166."""
    mag = onesite(M, PA)
    e = bond_energies(M, h)
    i = L // 4
    js = sorted(set([L // 2, 3 * L // 4 - 1, 3 * L // 4]))
    C = correlator_pairs(M, i, js, PA, PA)
    c = L // 2
    w = L // 8
    emid = float(np.mean(e[c - w:c + w]))
    return dict(mag=mag.tolist(), ebond=e.tolist(), emid=emid, corr_i=i, corr_js=js, corr=C.tolist())
