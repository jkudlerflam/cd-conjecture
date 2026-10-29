"""Attack (iii)/(iv): does the vacancy fluid move the cell-length selector?

Projected zero-port K = infinity model (C88 eq. (8), C91 gap form):
    H = t H_gap + J sum_j c_j b_j,   c_j = 1(g_j = 0),
with b_j the vacuum-clock bond operator of material bond j (unit clock hops, the
marker selector mu_sel split half to each adjacent bond).  For J << t the energy is
    E = J alpha sum_j <b_j> - (J^2/t) sum_ij S_ij <b_i b_j> + O(J^3/t^2),
with S the static contact susceptibility of the free vacancy fluid in label space.

Step 1: S_ij on gap rings N = 6..12 (exact sparse ED), fitted to a local part
        (s0, s1) plus a uniform -A/N piece.
Step 2: validate the second-order formula against exact ED of the coupled
        gap-fluid x clock model on small rings (two cells of length 3).
Step 3: apply to b = 14, 15, 16 cells at the scheduled mu_sel and find the J/t at
        which the selected length moves off 15.
"""
import json
import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

sys.path.insert(0, '.')
from gapfluid import compositions, hop_matrix, ground, response_vector, free_fermion_alpha

out = {}

# ---------------- step 1: susceptibility ----------------
susc = {}
for N in (6, 8, 10, 12):
    t0 = time.time()
    G = compositions(N, N)
    H = hop_matrix(G)
    E, V = ground(H)
    psi = V[:, 0]
    c = [(G[:, j] == 0).astype(float) for j in range(N)]
    alpha = float(psi ** 2 @ c[0])
    x = response_vector(H, E[0], psi, c[0])
    S = [float((c[r] * psi) @ x) for r in range(N)]
    A = -S[N // 2] * N
    s0 = S[0] + A / N
    s1 = S[1] + A / N
    susc[N] = dict(dim=len(G), E0=float(E[0]), gap=float(E[1] - E[0]), alpha=alpha,
                   alpha_free_fermion=float(free_fermion_alpha(2 * N, N)),
                   S=S, A=A, s0=s0, s1=s1, S_q0=float(sum(S)),
                   tail_spread=float(max(S[2:N - 1]) - min(S[2:N - 1])),
                   seconds=round(time.time() - t0, 1))
    print(N, susc[N], flush=True)
out['susceptibility'] = susc
s0, s1, A = susc[12]['s0'], susc[12]['s1'], susc[12]['A']
alpha_inf = 0.5 - 2 / np.pi ** 2


# ---------------- clock operators ----------------
def clock_bonds(l, mu):
    """Bond operators b_0..b_{l-1} of one cell of l material letters (p=l-1 data + marker).

    Clock positions s = 0..p (dimension l). Bond k < p carries the hop s=k <-> k+1;
    bond p-1 and bond p each carry mu/2 from the marker.
    """
    p = l - 1
    dim = p + 1
    bs = []
    for k in range(l):
        b = np.zeros((dim, dim))
        if k < p:
            b[k, k + 1] = b[k + 1, k] = -1.0
        if k in (p - 1, p):
            b += 0.5 * mu * np.eye(dim)
        bs.append(b)
    return bs


def cell_ground(l, mu):
    bs = clock_bonds(l, mu)
    h = sum(bs)
    w, v = np.linalg.eigh(h)
    return w[0], v[:, 0], bs


# ---------------- step 2: validate against exact ED ----------------
def coupled_ed(l, ncell, Nv, mu, J, t):
    N = l * ncell
    G = compositions(N, Nv)
    Hg = hop_matrix(G)
    bs = clock_bonds(l, mu)
    dim = l
    Ic = sp.identity(dim ** ncell, format='csr')
    Hm_terms = []
    for cell in range(ncell):
        for k in range(l):
            ops = [sp.identity(dim, format='csr')] * ncell
            ops = list(ops)
            ops[cell] = sp.csr_matrix(bs[k])
            M = ops[0]
            for o in ops[1:]:
                M = sp.kron(M, o, format='csr')
            j = cell * l + k
            cdiag = sp.diags((G[:, j] == 0).astype(float))
            Hm_terms.append(sp.kron(cdiag, M, format='csr'))
    H = t * sp.kron(Hg, Ic, format='csr') + J * sum(Hm_terms)
    w = spl.eigsh(H, k=1, which='SA', tol=1e-14)[0][0]
    return w, G, Hg


def second_order_finite(l, ncell, Nv, mu, J, t):
    """alpha J sum<b> - (J^2/t) sum S_ij <b_i b_j> using the exact finite-ring S."""
    N = l * ncell
    G = compositions(N, Nv)
    Hg = hop_matrix(G)
    E, V = ground(Hg)
    psi = V[:, 0]
    c = [(G[:, j] == 0).astype(float) for j in range(N)]
    al = float(psi ** 2 @ c[0])
    x = response_vector(Hg, E[0], psi, c[0])
    S = np.array([float((c[r] * psi) @ x) for r in range(N)])
    e, v, bs = cell_ground(l, mu)
    # material state: product of cell ground states (ground state of alpha*sum b)
    eb = np.array([v @ b @ v for b in bs])
    bb = np.array([[v @ bs[a] @ bs[b] @ v for b in range(l)] for a in range(l)])
    tot = 0.0
    for i in range(N):
        for j in range(N):
            ci, ki = divmod(i, l)
            cj, kj = divmod(j, l)
            corr = bb[ki, kj] if ci == cj else eb[ki] * eb[kj]
            tot += S[(j - i) % N] * corr
    return t * E[0] + J * al * ncell * e - J ** 2 / t * tot, t * E[0] + J * al * ncell * e


val = []
for (l, ncell, Nv, mu) in [(3, 2, 6, 1.3), (4, 2, 8, 1.5)]:
    for Jt in (0.02, 0.05, 0.1):
        J, t = Jt, 1.0
        t0 = time.time()
        Eex, _, _ = coupled_ed(l, ncell, Nv, mu, J, t)
        E2, E1 = second_order_finite(l, ncell, Nv, mu, J, t)
        rec = dict(l=l, ncell=ncell, Nv=Nv, mu=mu, J_over_t=Jt, E_exact=Eex, E_first=E1,
                   E_second=E2, second_order_coeff_exact=(Eex - E1) / J ** 2,
                   second_order_coeff_pred=(E2 - E1) / J ** 2,
                   residual_over_J3=(Eex - E2) / J ** 3, seconds=round(time.time() - t0, 1))
        print(rec, flush=True)
        val.append(rec)
out['validation'] = val


# ---------------- step 3: the scheduled selector ----------------
def a_(l):
    return 2 * (1 - np.cos(np.pi / (l + 1)))


def delta_(l):
    return (l + 1) * a_(l) - l * a_(l + 1)


mu_sel = 2 - (delta_(15) + delta_(14)) / 2


def per_site(l, s0, s1, A):
    """Return (first-order coefficient, second-order coefficient) per material site.

    e(l) = J * alpha * c1 - (J^2/t) * c2.
    """
    e, v, bs = cell_ground(l, mu_sel)
    eb = np.array([v @ b @ v for b in bs])
    bb = np.array([[v @ bs[a] @ bs[b] @ v for b in range(l)] for a in range(l)])
    c1 = e / l
    loc = 0.0
    for k in range(l):
        loc += s0 * bb[k, k]
        # neighbour k+1, in the same cell or the next one
        if k + 1 < l:
            loc += 2 * s1 * bb[k, k + 1]
        else:
            loc += 2 * s1 * eb[k] * eb[0]
    ebar = e / l
    c2 = loc / l - A * ebar ** 2
    return c1, c2, eb.tolist()


sel = {}
for l in range(10, 22):
    c1, c2, eb = per_site(l, s0, s1, A)
    sel[l] = dict(first=c1, second=c2)
res = {}
for l in range(10, 22):
    if l == 15:
        continue
    d1 = sel[l]['first'] - sel[15]['first']
    d2 = sel[l]['second'] - sel[15]['second']
    # e_l - e_15 = J alpha d1 - (J^2/t) d2 ; vanishes at J/t = alpha d1/d2 if d2 > 0
    thr = alpha_inf * d1 / d2 if d2 > 0 else None
    res[l] = dict(first_order_gap=d1, second_order_shift=d2, J_over_t_threshold=thr)
    print(l, res[l], flush=True)
out['scheduled_selector'] = dict(mu_sel=mu_sel, alpha_contact=alpha_inf, s0=s0, s1=s1, A=A,
                                 per_length=res,
                                 bond_energies_b15=per_site(15, s0, s1, A)[2])
thr = [r['J_over_t_threshold'] for r in res.values() if r['J_over_t_threshold']]
out['scheduled_selector']['J_over_t_critical'] = min(thr) if thr else None
print('critical J/t', out['scheduled_selector']['J_over_t_critical'])
json.dump(out, open('selector_shift.json', 'w'), indent=1, default=float)

# ---------------- fixed-marker-number comparison ----------------
# At fixed N and M (both conserved) wrong lengths come in pairs (15-k, 15+k).
pair = {}
for k in (1, 2, 3):
    lo, hi = 15 - k, 15 + k
    r = out['scheduled_selector']['per_length']
    d1 = lo * r[lo]['first_order_gap'] + hi * r[hi]['first_order_gap']
    d2 = lo * r[lo]['second_order_shift'] + hi * r[hi]['second_order_shift']
    pair[f'{lo}+{hi}'] = dict(first_order_cost_coeff=d1, second_order_coeff=d2,
                              note='cost = J alpha first - (J^2/t) second')
    print(k, pair[f'{lo}+{hi}'])
out['fixed_M_pairs'] = pair
json.dump(out, open('selector_shift.json', 'w'), indent=1, default=float)
