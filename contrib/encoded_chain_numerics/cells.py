#!/usr/bin/env python3
"""Route-B ("dressed Bell pair") cell of contrib/programmable_cell, chained by its ports.

One cell = one big site of dimension D = q^4 2^n:
    pL (q) x sL (q) x logical (2^n) x sR (q) x pR (q).
    H_L = Delta (1 - |Om><Om|)_{pL sL} + g sum_a (|Om_a><Om| + h.c.)_{pL sL} x gamma_a
    H_R = Delta (1 - |Om><Om|)_{pR sR} + g sum_a (|Om_a><Om| + h.c.)_{pR sR} x eta_a,  eta_a = i gamma_a P
    |Om> maximally entangled on (port, partner), |Om_a> = (F_a x 1)|Om>, F_a the r = 2 q+ q- graded
    (odd) port channels with tr(F_a F_b) = q delta_ab (graded_channels of majorana_bank_check.py).
Chain:  H = sum_j H_cell(j) + lam sum_j sum_a F_a(pR_j) F_a(pL_{j+1}).

For tensor networks the big site is factorized EXACTLY into three sub-sites
    [pL (q), mid = sL x logical x sR (q^2 2^n), pR (q)],
so H_L, H_R and the inter-cell bond are all nearest-neighbour two-sub-site terms.  Because the
outer sub-sites are q-dimensional, an MPS with bond chi at the cell boundary and q*chi on the two
inner bonds is exactly a big-site MPS of bond dimension chi (no extra truncation).  Entropies are
quoted at the cell boundary (bond 0 of the cell).  bigsite_* builds the literal D x D site.

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import os, sys
import numpy as np
from scipy import linalg

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'programmable_cell'))
from majorana_bank_check import majoranas, pair_states, graded_channels  # noqa: E402

kron = np.kron


def swap_pair(q):
    """Permutation matrix (p x s) -> (s x p)."""
    Sw = np.zeros((q * q, q * q))
    for i in range(q):
        for j in range(q):
            Sw[j * q + i, i * q + j] = 1
    return Sw


def cell(q, qp, n, Delta=1.0, g=0.5):
    """All ingredients of one cell.  Returns a dict.
    HL acts on pairL (port x partner) x logical; HR on logical x pairR (partner x port)."""
    F = graded_channels(q, qp)
    r = len(F)
    gam, P = majoranas(n)
    assert len(gam) >= r, 'logical register too small'
    gam = gam[:r]
    eta = [1j * x @ P for x in gam]
    N = 2 ** n
    S = pair_states(q, F)                    # (port x partner) columns |Om>, |Om_a>
    Sw = swap_pair(q)
    SR = Sw @ S                              # same states in (partner x port) order
    Iq2, IN = np.eye(q * q), np.eye(N)
    PomL = np.outer(S[:, 0], S[:, 0].conj())
    PomR = np.outer(SR[:, 0], SR[:, 0].conj())
    HL = Delta * kron(Iq2 - PomL, IN)
    HR = Delta * kron(IN, Iq2 - PomR)
    for a in range(r):
        K = np.outer(S[:, a + 1], S[:, 0].conj())
        HL = HL + g * kron(K + K.conj().T, gam[a])
        K = np.outer(SR[:, a + 1], SR[:, 0].conj())
        HR = HR + g * kron(eta[a], K + K.conj().T)
    Iq = np.eye(q)
    FportL = [kron(f, Iq) for f in F]        # port operators on pairL (port first)
    FportR = [kron(Iq, f) for f in F]        # port operators on pairR (port last)
    Zq = np.diag([1.0] * qp + [-1.0] * (q - qp))
    out = dict(q=q, qp=qp, n=n, r=r, N=N, Delta=Delta, g=g, F=F, gam=gam, eta=eta, P=P,
               HL=HL, HR=HR, FportL=FportL, FportR=FportR, PomL=PomL,
               # three-sub-site factorization [pL, mid = sL x log x sR, pR]
               dims=[q, q * q * N, q], HLs=kron(HL, Iq), HRs=kron(Iq, HR))
    return out


def bigsite_H(c):
    """D x D cell Hamiltonian on pairL x logical x pairR."""
    q2, N = c['q'] ** 2, c['N']
    return kron(c['HL'], np.eye(q2)) + kron(np.eye(q2), c['HR'])


def bigsite_ports(c):
    q2, N = c['q'] ** 2, c['N']
    FL = [kron(f, np.eye(N * q2)) for f in c['FportL']]
    FR = [kron(np.eye(q2 * N), f) for f in c['FportR']]
    return FL, FR


# ----------------------------------------------------------------------------- MPOs
def term_fragments(H, tdims, tol=1e-13):
    """Split an operator on kron(tdims) (leftmost most significant) into MPO pieces
    ops[i] of shape (k_i, d_i, d_i, k_{i+1}) with indices (left, out, in, right)."""
    m = len(tdims)
    T = np.asarray(H, complex).reshape(list(tdims) + list(tdims))
    perm = [x for i in range(m) for x in (i, m + i)]
    T = T.transpose(perm)
    ops, left = [], 1
    M = T.reshape(1, -1)
    for i in range(m - 1):
        M = M.reshape(left * tdims[i] ** 2, -1)
        U, S, Vh = linalg.svd(M, full_matrices=False)
        keep = max(1, int(np.sum(S > tol * S[0])))
        ops.append(U[:, :keep].reshape(left, tdims[i], tdims[i], keep))
        M = S[:keep, None] * Vh[:keep]
        left = keep
    ops.append(M.reshape(left, tdims[-1], tdims[-1], 1))
    return ops


def mpo_from_terms(dims, terms, tol=1e-13):
    """Periodic MPO (unit cell of len(dims) sites) for H = sum over cells of the given terms.
    terms: list of (start_site, nsites, matrix).  Channel 0 = nothing placed, last = done.
    Returns W[k] of shape (D_k, D_{k+1}, d_k, d_k) with bond k to the left of site k."""
    n = len(dims)
    nch = [0] * n
    frags, offs = [], {}
    for t, (s, m, H) in enumerate(terms):
        tdims = [dims[(s + i) % n] for i in range(m)]
        ops = term_fragments(H, tdims, tol)
        frags.append((s, ops))
        for i in range(1, m):
            b = (s + i) % n
            offs[(t, i)] = 1 + nch[b]
            nch[b] += ops[i].shape[0]
    D = [2 + x for x in nch]
    W = [np.zeros((D[k], D[(k + 1) % n], dims[k], dims[k]), complex) for k in range(n)]
    for k in range(n):
        W[k][0, 0] = np.eye(dims[k])
        W[k][-1, -1] = np.eye(dims[k])
    for t, (s, ops) in enumerate(frags):
        m = len(ops)
        for i, op in enumerate(ops):
            k = (s + i) % n
            kl, _, _, kr = op.shape
            for a in range(kl):
                wl = 0 if i == 0 else offs[(t, i)] + a
                for b in range(kr):
                    wr = D[(k + 1) % n] - 1 if i == m - 1 else offs[(t, i + 1)] + b
                    W[k][wl, wr] += op[a, :, :, b]
    return W


def chain_terms(c, lam, ncells=1):
    """Terms for a unit cell of `ncells` cells, three sub-sites each."""
    dims = c['dims'] * ncells
    bond = sum(kron(f, f) for f in c['F'])          # F_a(pR_j) F_a(pL_{j+1})
    terms = []
    for j in range(ncells):
        terms.append((3 * j, 2, c['HLs']))
        terms.append((3 * j + 1, 2, c['HRs']))
        terms.append((3 * j + 2, 2, lam * bond))
    return dims, terms


def chain_mpo(c, lam, ncells=1):
    dims, terms = chain_terms(c, lam, ncells)
    return dims, mpo_from_terms(dims, terms)


def bigsite_mpo(c, lam, ncells=1):
    """The literal big-site chain (site dimension D = q^4 2^n), for cross-checks."""
    H = bigsite_H(c)
    FL, FR = bigsite_ports(c)
    bond = sum(kron(fr, fl) for fr, fl in zip(FR, FL))
    D = H.shape[0]
    terms = []
    for j in range(ncells):
        terms.append((j, 1, H))
        terms.append((j, 2, lam * bond))
    return [D] * ncells, mpo_from_terms([D] * ncells, terms)


# ----------------------------------------------------------------------------- logical chain
def logical_chain_terms(c, t, ncells=1):
    """First-order logical chain t sum_a eta_a(j) gamma_a(j+1) on the 2^n-dim register."""
    h = sum(kron(e, g) for e, g in zip(c['eta'], c['gam']))
    return [c['N']] * ncells, [(j, 2, t * h) for j in range(ncells)]


def ising_terms(J=1.0, hx=1.0):
    X = np.array([[0, 1], [1, 0]], complex); Z = np.diag([1.0, -1.0]).astype(complex)
    return [2], [(0, 2, -J * kron(Z, Z)), (0, 1, -hx * X)]


# ----------------------------------------------------------------------------- cell checks
def cell_spectrum_and_alpha(c):
    """Exact cell spectrum and the compressions of the ports onto the 2^n ground space."""
    H = bigsite_H(c)
    E, V = linalg.eigh(H)
    N = c['N']
    G = V[:, :N]
    # fix the logical frame through the |Om>_L |Om>_R component (proportional to psi)
    q2 = c['q'] ** 2
    blk = G.reshape(q2, N, q2, N)
    i0 = 0  # |Om> is column 0 of pair_states, but here the pair basis is computational
    S = pair_states(c['q'], c['F'])
    SR = swap_pair(c['q']) @ S
    comp = np.einsum('a,aibk,b->ik', S[:, 0].conj(), blk, SR[:, 0].conj())
    U, s, Wh = linalg.svd(comp)
    G = G @ (U @ Wh).conj().T
    FL, FR = bigsite_ports(c)
    aL, aR, res = [], [], 0.0
    for b in range(c['r']):
        CL = G.conj().T @ FL[b] @ G
        CR = G.conj().T @ FR[b] @ G
        al = np.trace(c['gam'][b].conj().T @ CL).real / N
        ar = np.trace(c['eta'][b].conj().T @ CR).real / N
        aL.append(al); aR.append(ar)
        res = max(res, np.linalg.norm(CL - al * c['gam'][b]), np.linalg.norm(CR - ar * c['eta'][b]))
    return dict(E=E, split=E[N - 1] - E[0], gap=E[N] - E[N - 1], alphaL=aL, alphaR=aR, resid=res, G=G)


if __name__ == '__main__':
    np.set_printoptions(precision=6, suppress=True, linewidth=140)
    for (q, qp, n) in [(2, 1, 1), (3, 1, 2)]:
        for g in (0.2, 0.5, 1.0):
            c = cell(q, qp, n, 1.0, g)
            o = cell_spectrum_and_alpha(c)
            D = (q ** 4) * 2 ** n
            Eu = np.unique(np.round(o['E'], 8))
            print(f"q={q} q+={qp} r={c['r']} n={n} D={D} g={g}: split={o['split']:.1e} gap={o['gap']:.4f} "
                  f"alphaL={np.round(o['alphaL'], 6)} alphaR={np.round(o['alphaR'], 6)} resid={o['resid']:.1e}")
            print('   lowest distinct cell levels:', Eu[:8])
        dims, W = chain_mpo(c, 0.1, 1)
        print('   sub-site dims', dims, 'MPO bond dims', [w.shape[0] for w in W])


def u1_charge(c):
    """Big-site U(1) charge Q = (Pi0_pL - Pi0_sL) + Q_log + (Pi0_pR - Pi0_sR), Pi0 = |0><0| on the
    q+ = 1 port state, Q_log = s (i/2) sum_k gamma_{2k} gamma_{2k+1}; the sign s is fixed by
    [Q, H_cell] = 0.  Exists for q+ = 1 (the channels E_0k carry charge 1).  Also returns the check
    that Q commutes with the inter-cell bond."""
    q, N = c['q'], c['N']
    assert c['qp'] == 1
    Pi0 = np.zeros((q, q)); Pi0[0, 0] = 1
    Iq, IN = np.eye(q), np.eye(N)
    gam = c['gam']
    Ql = 0.5j * sum(gam[2 * k] @ gam[2 * k + 1] for k in range(len(gam) // 2))
    H = bigsite_H(c)
    k5 = lambda a, b, cc, d, e: kron(kron(kron(kron(a, b), cc), d), e)
    best = None
    for s in (1, -1):
        Q = (k5(Pi0, Iq, IN, Iq, Iq) - k5(Iq, Pi0, IN, Iq, Iq) + s * k5(Iq, Iq, Ql, Iq, Iq)
             - k5(Iq, Iq, IN, Pi0, Iq) + k5(Iq, Iq, IN, Iq, Pi0))
        err = np.linalg.norm(Q @ H - H @ Q)
        if best is None or err < best[1]:
            best = (Q, err)
    Q = best[0]
    # bond check without forming D^2 x D^2 matrices: [Q_pR + Q_pL', F_a(pR) F_a(pL')] summed over a
    Fs = c['F']
    bond = sum(kron(f, f) for f in Fs)
    qq = kron(Pi0, Iq) + kron(Iq, Pi0)
    berr = np.linalg.norm(qq @ bond - bond @ qq)
    return Q, best[1], berr


def truncated_cell(c, M):
    """Cell restricted to its M lowest eigenstates (M must end at a gap of the cell spectrum).
    Returns (E_M, FL_M, FR_M, Q_M, gapcut): H_cell = diag(E_M), port operators and U(1) charge
    projected on the kept states, and the energy gap at the cut."""
    H = bigsite_H(c)
    E, V = linalg.eigh(H)
    VM = V[:, :M]
    FL, FR = bigsite_ports(c)
    FLm = [VM.conj().T @ f @ VM for f in FL]
    FRm = [VM.conj().T @ f @ VM for f in FR]
    Qm = VM.conj().T @ u1_charge(c)[0] @ VM if c['qp'] == 1 else None
    return E[:M], FLm, FRm, Qm, (E[M] - E[M - 1]) if M < len(E) else np.inf
