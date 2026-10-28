#!/usr/bin/env python3
"""Task 2: a route-B cell whose ports carry the species-torus charges, and its exact
second-order effective chain.

Cell (as in ../encoded_chain_numerics/cells.py, Delta = 1): one big site
    pL (q) x sL (q) x logical (2^n) x sR (q) x pR (q),  q = 4, q+ = 2, n = 4, D = 4096.
The graded port C^{2|2} is the torus data register of contrib/chat_followup/torus_encoding.md
with b = 2r = 4: upper states |1,alpha> (port states 0, 1), lower states |0,beta> (2, 3).
Its r_M = 2 q+ q- = 8 odd Hermitian channels are the real and imaginary parts of
F_{alpha beta} = |1,alpha><0,beta|, so the logical register carries N = r^2 = 4 complex
flavours a = 2 alpha + beta (qubit a of the register, Majoranas gamma_{2a}, gamma_{2a+1}).
  H_L = Delta (1 - |Om><Om|)_L + sum_A gL_A (|Om_A><Om| + h.c.)_L x gamma_A
  H_R = same on the right pair with eta_A = i gamma_A P and gR_A
  bond  lam sum_A w_A F_A(pR_j) F_A(pL_{j+1})
Channel A = 2a, 2a+1 belongs to flavour a; gL, gR and w depend on the flavour only, which keeps
every species phase rotation (the torus) and breaks the U(2) x U(2) of the uniform cell.
A common chemical potential mu0 Q, Q the total charge (a torus element, conserved by every term),
is added exactly after Schrieffer-Wolff.

Second order (Schrieffer-Wolff, exact through lam^2, same formulas as effective.py, generalized
to channel-dependent compressions and bond weights):
  T1 = sum_A w_A CR_A x CL_A,       CR_A = G^ F^R_A G, CL_A = G^ F^L_A G
  T2 = sum_{(m,n) != code} <ik|V|mn><mn|V|jl> / (2 e_c - E_m - E_n)          (two cells)
  T3 = sum_AB w_A w_B CR_A x G^(F^R_B R F^L_A + F^L_A R F^R_B)G x CL_B,  R = sum_exc |m><m|/(e_c - E_m)
Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'encoded_chain_numerics'))
sys.path.insert(0, os.path.join(HERE, '..', 'redteam_construction'))
import numpy as np
from scipy import linalg
import cells as C                                   # noqa: E402
from majorana_bank_check import majoranas, pair_states, graded_channels   # noqa: E402
import route_b_second_order as RB                   # noqa: E402

kron = np.kron


def flavour_of_channel(r):
    return [A // 2 for A in range(r)]


def build_cell(q=4, qp=2, n=4, Delta=1.0, gL=0.5, gR=0.5):
    """gL, gR: scalars or per-flavour arrays (length r/2)."""
    F = graded_channels(q, qp)
    r = len(F)
    nf = r // 2
    gL = np.broadcast_to(np.asarray(gL, float), (nf,))
    gR = np.broadcast_to(np.asarray(gR, float), (nf,))
    gam, P = majoranas(n)
    gam = gam[:r]
    eta = [1j * x @ P for x in gam]
    N = 2 ** n
    S = pair_states(q, F)
    SR = C.swap_pair(q) @ S
    Iq2, IN = np.eye(q * q), np.eye(N)
    PomL = np.outer(S[:, 0], S[:, 0].conj())
    PomR = np.outer(SR[:, 0], SR[:, 0].conj())
    HL = Delta * kron(Iq2 - PomL, IN)
    HR = Delta * kron(IN, Iq2 - PomR)
    for A in range(r):
        K = np.outer(S[:, A + 1], S[:, 0].conj())
        HL = HL + gL[A // 2] * kron(K + K.conj().T, gam[A])
        K = np.outer(SR[:, A + 1], SR[:, 0].conj())
        HR = HR + gR[A // 2] * kron(eta[A], K + K.conj().T)
    Iq = np.eye(q)
    return dict(q=q, qp=qp, n=n, r=r, N=N, Delta=Delta, gL=gL, gR=gR, F=F, gam=gam, eta=eta, P=P,
                HL=HL, HR=HR, FportL=[kron(f, Iq) for f in F], FportR=[kron(Iq, f) for f in F], PomL=PomL)


# ----------------------------------------------------------------------------- torus generators
def logical_numbers(c):
    """n_a = f_a^dag f_a with f_a = (gamma_2a + i gamma_2a+1)/2 (qubit a, occupied = |1>)."""
    g = c['gam']
    N = c['N']
    out = []
    for a in range(c['r'] // 2):
        fa = (g[2 * a] + 1j * g[2 * a + 1]) / 2
        out.append(fa.conj().T @ fa)
    return out


def torus_charge_matrix(qp=2, q=4):
    """W[a, s]: torus charge of f_a^dag on species s (s = port basis state)."""
    nf = qp * (q - qp)
    W = np.zeros((nf, q), int)
    for al in range(qp):
        for be in range(q - qp):
            a = al * (q - qp) + be
            W[a, al] = 1; W[a, qp + be] = -1
    return W


def bigsite_torus(c, sign):
    """Species charges on the big site: T_s = Pi_s(pL) - Pi_s(sL) + sign*Qlog_s + Pi_s(pR) - Pi_s(sR)."""
    q, N = c['q'], c['N']
    W = torus_charge_matrix(c['qp'], q)
    nl = logical_numbers(c)
    Iq, IN = np.eye(q), np.eye(N)
    k5 = lambda a, b, cc, d, e: kron(kron(kron(kron(a, b), cc), d), e)
    Ts = []
    for s in range(q):
        Pi = np.zeros((q, q)); Pi[s, s] = 1
        Ql = sum(W[a, s] * nl[a] for a in range(len(nl)))
        # big-site order: pL, sL, logical, sR, pR
        Ts.append(k5(Pi, Iq, IN, Iq, Iq) - k5(Iq, Pi, IN, Iq, Iq) + sign * k5(Iq, Iq, Ql, Iq, Iq)
                  - k5(Iq, Iq, IN, Pi, Iq) + k5(Iq, Iq, IN, Iq, Pi))
    return Ts


def check_cell_torus(c, w=None):
    """Find the sign making [T_s, H_cell] = 0; check the bond (two-port operator) too."""
    H = C.bigsite_H(c)
    best = None
    for sign in (1, -1):
        Ts = bigsite_torus(c, sign)
        err = max(np.linalg.norm(T @ H - H @ T) for T in Ts)
        if best is None or err < best[1]:
            best = (sign, err, Ts)
    q = c['q']
    w = np.ones(c['r'] // 2) if w is None else w
    bond = sum(w[A // 2] * kron(f, f) for A, f in enumerate(c['F']))
    berr = 0.0
    for s in range(q):
        Pi = np.zeros((q, q)); Pi[s, s] = 1
        tt = kron(Pi, np.eye(q)) + kron(np.eye(q), Pi)
        berr = max(berr, np.linalg.norm(tt @ bond - bond @ tt))
    return best[0], best[1], berr, best[2]


# ----------------------------------------------------------------------------- spectrum, frame
def code_and_frame(c):
    H = C.bigsite_H(c)
    E, V = linalg.eigh(H, driver='evr')
    N = c['N']
    G = V[:, :N]
    q2 = c['q'] ** 2
    blk = G.reshape(q2, N, q2, N)
    S = pair_states(c['q'], c['F'])
    SR = C.swap_pair(c['q']) @ S
    comp = np.einsum('a,aibk,b->ik', S[:, 0].conj(), blk, SR[:, 0].conj())
    U, s, Wh = linalg.svd(comp)
    G = G @ (U @ Wh).conj().T
    return E, V, G


def second_order(c, w=None, E=None, V=None, G=None):
    """T1, T2, T3 and the compressions (bond weights folded in: CRw_A = w_A CR_A)."""
    if E is None:
        E, V, G = code_and_frame(c)
    N, r = c['N'], c['r']
    w = np.ones(r // 2) if w is None else np.asarray(w, float)
    wA = np.array([w[A // 2] for A in range(r)])
    FL, FR = C.bigsite_ports(c)
    D = V.shape[0]
    B = np.hstack([G, V[:, N:]])
    ec = E[0]
    Eb = np.concatenate([np.full(N, ec), E[N:]])
    Lm = [G.conj().T @ (FL[A] @ B) for A in range(r)]       # (N, D)
    Rm = [G.conj().T @ (FR[A] @ B) for A in range(r)]
    CL = [x[:, :N] for x in Lm]
    CR = [x[:, :N] for x in Rm]
    # compressions: CL_A = aL_A gamma_A, CR_A = aR_A eta_A
    aL = np.array([np.trace(c['gam'][A].conj().T @ CL[A]).real / N for A in range(r)])
    aR = np.array([np.trace(c['eta'][A].conj().T @ CR[A]).real / N for A in range(r)])
    res = max(max(np.linalg.norm(CL[A] - aL[A] * c['gam'][A]) for A in range(r)),
              max(np.linalg.norm(CR[A] - aR[A] * c['eta'][A]) for A in range(r)))
    # T2 = sum_AB wA wB sum_mn Rm_A[i,m] conj(Rm_B[j,m]) W[m,n] Lm_A[k,n] conj(Lm_B[l,n])
    Wd = 1.0 / (2 * ec - Eb[:, None] - Eb[None, :])
    Wd[:N, :N] = 0.0
    T2 = np.zeros((N, N, N, N), complex)     # [i, j, k, l]
    for A in range(r):
        XA = np.stack([np.einsum('im,jm->ijm', Rm[A], Rm[Bb].conj()).reshape(N * N, D) * (wA[A] * wA[Bb])
                       for Bb in range(r)])                                  # (r, N^2, D)
        YA = np.stack([np.einsum('kn,ln->kln', Lm[A], Lm[Bb].conj()).reshape(N * N, D) for Bb in range(r)])
        for Bb in range(r):
            Xr = XA[Bb]
            XW = Xr.real @ Wd + 1j * (Xr.imag @ Wd)                         # (N^2, D)
            T2 += (XW @ YA[Bb].T).reshape(N, N, N, N)
    T2 = T2.transpose(0, 2, 1, 3).reshape(N * N, N * N)                      # [(i,k), (j,l)]
    T2 = (T2 + T2.conj().T) / 2
    # T3
    exc = np.arange(N, D)
    rinv = 1.0 / (ec - E[N:])
    T3 = np.zeros((N ** 3, N ** 3), complex)
    YRL = np.zeros((r, r, N, N), complex)    # red-team convention: [b, a] = G^ F^R_b R+ F^L_a G, R+ = sum 1/eps
    YLR = np.zeros((r, r, N, N), complex)
    for A in range(r):
        for Bb in range(r):
            M1 = (Rm[Bb][:, exc] * rinv) @ Lm[A][:, exc].conj().T
            M2 = (Lm[A][:, exc] * rinv) @ Rm[Bb][:, exc].conj().T
            YRL[Bb, A] = -M1
            YLR[A, Bb] = -M2
            T3 += wA[A] * wA[Bb] * kron(kron(CR[A], M1 + M2), CL[Bb])
    T3 = (T3 + T3.conj().T) / 2
    T1 = sum(wA[A] * kron(CR[A], CL[A]) for A in range(r))
    return dict(T1=T1, T2=T2, T3=T3, CR=CR, CL=CL, aL=aL, aR=aR, resid=res, wA=wA, ec=ec, E=E,
                YRL=YRL, YLR=YLR, split=E[N - 1] - E[0], gap=E[N] - E[N - 1])


# ----------------------------------------------------------------------------- logical torus checks
def logical_torus_site(c):
    W = torus_charge_matrix(c['qp'], c['q'])
    nl = logical_numbers(c)
    return [sum(W[a, s] * nl[a] for a in range(len(nl))) for s in range(c['q'])]


def check_logical_torus(c, T, m):
    Ts = logical_torus_site(c)
    I = np.eye(c['N'])
    err = 0.0
    for t in Ts:
        tm = sum(C.kron(C.kron(np.eye(c['N'] ** x), t), np.eye(c['N'] ** (m - x - 1))) for x in range(m))
        err = max(err, np.linalg.norm(T @ tm - tm @ T) / max(np.linalg.norm(T), 1e-300))
    return err


# ----------------------------------------------------------------------------- Majorana expansion
class RBView:
    """Duck-typed object for route_b_second_order.coefficient / unit_cell_terms."""
    def __init__(self, c, so, lam):
        self.n, self.N, self.r = c['n'], c['N'], c['r']
        self.P, self.gam = c['P'], c['gam']
        self.S2 = so['T2']
        self.CR = [so['wA'][A] * so['CR'][A] for A in range(self.r)]
        self.CL = [so['wA'][A] * so['CL'][A] for A in range(self.r)]
        self.YRL, self.YLR = so['YRL'], so['YLR']


def quadratic_terms(c, so, lam, mu0, kappa=1.0):
    """All Majorana bilinears of lam T1 + kappa lam^2 (T2 + T3) + mu0 Q per unit cell, as a dict
    ((j,a),(j',b)) -> coefficient, j in {0,1,2} with the window logic of route_b_second_order."""
    view = RBView(c, so, lam)
    q2 = RB.unit_cell_terms(view, 2)
    out = {k: kappa * lam ** 2 * v for k, v in q2.items()}
    N = c['N']
    # first order: T1 on cells (0,1)
    for a in range(c['r']):
        for b in range(c['r']):
            mono = ((0, a), (1, b))
            f0, f1, f2 = RB.spin_form(view, mono)
            v = np.trace(np.kron(f0, f1).conj().T @ so['T1']) / N ** 2
            if abs(v) > 1e-13:
                out[mono] = out.get(mono, 0) + lam * v
    # mu0 Q = mu0 sum_a n_a = mu0 sum_a (1 + i gamma_2a gamma_2a+1 ...)/2: expand on cell 1
    nl = logical_numbers(c)
    Q = sum(nl)
    for a in range(c['r']):
        for b in range(a + 1, c['r']):
            mono = ((1, a), (1, b))
            f0, f1, f2 = RB.spin_form(view, mono)
            v = np.trace(f1.conj().T @ Q) / N
            if abs(v) > 1e-13:
                out[mono] = out.get(mono, 0) + mu0 * v
    return out


def flavour_dispersion(quad, m, ks):
    """Bloch Hamiltonian h(k) (m x m, Majoranas) from bilinears kappa c_x c_y, and the complex
    single-particle energies per flavour: block (2a, 2a+1) of h(k).  Returns eps[k, a] (the
    eigenvalue of the flavour block that is odd under k -> -k paired)."""
    out = []
    for k in ks:
        h = np.zeros((m, m), complex)
        for ((j1, a), (j2, b)), kap in quad.items():
            A = -2j * kap
            d = j2 - j1
            h[a, b] += 1j * A * np.exp(1j * k * d)
            h[b, a] += -1j * A * np.exp(-1j * k * d)
        out.append(h)
    return np.array(out)
