"""Quick reproduction of the numerical claims of RESULT.md and PROOF.md (contrib/cd_crossed_fock).
Run: OMP_NUM_THREADS=2 python3 run_checks.py      (a few minutes on a laptop; numpy and scipy only)
Each check prints PASS or FAIL with the numbers; results go to checks.json.
Heavier tables (m = 7, 8 spectra, m = 11 Lanczos, all-model insertion scans) are produced by the scripts named
in PROOF.md section 7; their stored outputs in data/ are spot-checked here from scratch."""
import os, sys, json, time
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cflib import Column, Transfer, sector_spectrum, direct_logZ, direct_logZ_dense, bond_xxz, bond_xxz_rot, realign, kron_power
from models import MODELS, get
from majorant import F, xx_logZt, vmin_needed
from insertion import Sectors, single_insertion, one_particle_levels, fock_counts
import xx_theorem
from eps0_lanczos import make_op
from scipy.sparse.linalg import eigsh

out = {}
ok_all = True


def report(name, ok, **vals):
    global ok_all
    ok_all &= bool(ok)
    out[name] = dict(ok=bool(ok), **{k: (float(v) if isinstance(v, (np.floating, float, int)) else v) for k, v in vals.items()})
    print(('PASS ' if ok else 'FAIL ') + name + ' ' + json.dumps(out[name], default=str), flush=True)


t0 = time.time()
# T1: crossed traces = direct Trotter traces, all models; T2: realigned gates PSD
worst = 0.0
for name in MODELS:
    h, mode, M, A, label = get(name)
    m, b = 3, 2.0
    col = Column(2, m); T = Transfer(h, 2, b, m, col); Qs = col.charge(mode)
    u = np.concatenate([sector_spectrum(T, np.where(Qs == q)[0], n) for q in np.unique(Qs) for n in range(m)])
    assert len(u) == col.N
    for L in (2, 4, 6):
        worst = max(worst, abs(logsumexp(L / 2 * np.log(np.clip(u, 1e-300, None))) - direct_logZ_dense(h, L, b, m)))
report('T1_crossed_equals_direct_trace', worst < 1e-10, max_abs_error=worst)

# X1: XX exact crossed product against exact diagonalization
err2, err4 = 0.0, 0.0
for b in (0.5, 2.0, 8.0):
    pc = quad(lambda k: np.logaddexp(b * np.cos(k) / 2, -b * np.cos(k) / 2), 0, 2 * np.pi, limit=200)[0] / (2 * np.pi)
    for L in (2, 4, 6, 8):
        e = abs(direct_logZ(bond_xxz(0.), L, b) - L * pc - xx_logZt([float(L)], b)[0])
        if L == 2: err2 = max(err2, e)
        else: err4 = max(err4, e)
report('X1_xx_product_vs_ED', err2 < 1e-6 and err4 < 1e-10, err_L2=err2, err_L_ge_4=err4)

# X2: XX theorem constants and the inequality on a reduced grid
res = {}
okx = True
for X, r, Kref in [(2.0, 1.5, 2.0384), (2.0, 2.0, 1.4413)]:
    K, xm, left, right = xx_theorem.K_of(X, r, grid=np.geomspace(1e-4, 1e3, 40001))
    V1 = X / np.arcsinh(X); Vp = r * V1
    worst = -np.inf
    for b in np.geomspace(0.05, 2048, 25):
        Ls = np.array([4, 6, 8, 16, 32, 64, 256, 1024])
        g = np.concatenate([xx_logZt([2.], b), xx_logZt(Ls, b)]); LL = np.concatenate([[2], Ls])
        worst = max(worst, np.max(g - K - 2 * F(LL / (Vp * b))))
    okx &= abs(K - Kref) < 2e-3 and worst < 0
    res['X%g_r%g' % (X, r)] = dict(Vprime_over_M=2 * Vp, K=K, max_violation=worst)
report('X2_xx_theorem', okx, **{k: str({a: round(c, 4) for a, c in v.items()}) for k, v in res.items()})

# X3: XX exact optimal majorant, S = 2, K = log 2
Ls = np.array([2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 128, 256, 512])
betas = np.geomspace(0.5, 4096, 40)
g = np.array([np.concatenate([xx_logZt([2.], b), xx_logZt(Ls[1:], b)]) for b in betas]).T
LL, BB = np.meshgrid(Ls, betas, indexing='ij')
v = vmin_needed(g, LL, BB, 2, np.log(2)).max() / 0.5
report('X3_xx_optimal_Vprime_S2_Kln2', abs(v - 4.08) < 0.03, Vprime_over_M=v)

# X4: finite-m XX crossed levels converge like m^-2 to the exact crossed spectrum (C125 eq. (7) and delta_beta)
from majorant import xx_delta_beta
b = 4.0; rows = []
for m in (6, 7):
    col = Column(2, m); T = Transfer(bond_xxz_rot(0.0), 2, b, m, col); Qs = col.charge('u1')
    u0 = sector_spectrum(T, np.where(Qs == 0)[0], 0)[0]
    e1 = -0.5 * np.log(sector_spectrum(T, np.where(Qs == 0)[0], 1)[0] / u0)
    q1 = -0.5 * np.log(sector_spectrum(T, np.where(Qs == 1)[0], 0)[0] / u0)
    rows.append((m, e1, q1))
ex1 = 2 * np.arcsinh(np.pi / b); exq = xx_delta_beta(b)
c1 = [(ex1 - e) * m ** 2 for m, e, _ in rows]; cq = [(exq - q) * m ** 2 for m, _, q in rows]
report('X4_xx_finite_m_levels_converge', abs(c1[0] - c1[1]) < 0.02 * abs(c1[0]) and abs(cq[0] - cq[1]) < 0.05 * abs(cq[0]),
       levels_n1=str([round(r[1], 5) for r in rows]), exact_n1=ex1, m2_times_diff_n1=str([round(c, 3) for c in c1]),
       levels_Q1=str([round(r[2], 5) for r in rows]), delta_beta=exq, m2_times_diff_Q1=str([round(c, 3) for c in cq]))

# Z1: zero-temperature deficits, D(L) maximal at L = 2 and decreasing (L <= 12 here)
from zeroT import E0_EXACT, ground
okz = True; dz = {}
for name in ['XX', 'XXZ05', 'HEIS', 'XXZ2', 'TFIM', 'XYI']:
    h, mode, M, A, label = get(name)
    Ls_ = [2, 4, 6, 8, 10, 12]
    E = [ground(h, L, False) for L in Ls_]
    D = np.array([L * (L * E0_EXACT[name] - e) for L, e in zip(Ls_, E)]) / M
    okz &= bool(np.all(np.diff(D) < 0))
    dz[name] = [round(x, 4) for x in D[[0, -1]]]
report('Z1_zeroT_deficit_decreasing', okz, D2_and_D12_over_M=str(dz))

# E1: eps0 Richardson extrapolation on XX against the exact p_c (beta = 4, m = 7..11)
b = 4.0
vals = []
for m in range(7, 12):
    A_, G = make_op(get('XX')[0], b, m)
    u = eigsh(A_, k=1, which='LA', v0=np.ones(4 ** m), tol=1e-13)[0][0]
    vals.append(-0.5 * np.log(u))
from infinite_m import richardson, xx_pc
e_inf, err, mm = richardson(list(range(7, 12)), vals, b)
report('E1_eps0_extrapolation_XX', abs(e_inf + xx_pc(b)) < 1e-7, extrapolated=e_inf, exact=-xx_pc(b), error=e_inf + xx_pc(b))

# I1: Heisenberg selection rule, composite singlet, level-by-level failure, Fock domination (m = 6, beta = 8)
S = Sectors('HEIS', 8.0, 6, h=bond_xxz_rot(1.0), Qmax=2)
W = single_insertion(S)
e00, f00 = S.eps((0, 0)), W[(0, 0)]['f']
iss = np.where((f00 < 1e-10) & (e00 > 1e-9))[0][0]
eps_s = e00[iss]
eps_t = min(l[0] for l in one_particle_levels(S, W))
partner = min(np.min(abs(S.eps((q, 0)) - eps_s)) for q in (-2, -1, 1, 2))
lev = one_particle_levels(S, W)
allE = np.sort(np.concatenate([S.eps(k) for k in S.data]))
FE, FQ, Fn = fock_counts(lev, 6, 3 * eps_t)
FEs = np.sort(FE)
kk = min(len(FEs), len(allE)); lam = np.min(allE[1:kk] / FEs[1:kk])
e1 = np.array([l[0] for l in lev])
dom = min(-np.log1p(-np.exp(-L * e1)).sum() - logsumexp(-L * allE) for L in (2, 4, 8))
okI = (f00[iss] < 1e-12 and partner > 1e-3 and eps_s < 2 * eps_t and dom > 0)
report('I1_heisenberg_composite', okI, eps_t=eps_t, eps_s=eps_s, ratio=eps_s / eps_t, overlap_singlet=f00[iss],
       nearest_partner_gap=partner, lambda_global=lam, eps_s_over_2eps_t=eps_s / (2 * eps_t), fock_margin_min_L248=dom)

# I2: Fock domination and lambda_global on all models at m = 5, beta = 4
from analysis_insertion import analyse
okI2 = True; li = {}
for name in MODELS:
    r = analyse(name, 4.0, 5)
    margin = min(a - b for a, b in zip(r['logZF'], r['logZt']))
    okI2 &= margin > -1e-9 and r['lambda_global'] >= 0.87
    li[name] = (round(r['lambda_global'], 3), round(margin, 4), round(r['fmin'], 3))
report('I2_fock_domination_all_models_m5_b4', okI2, lambda_global_margin_fmin=str(li))

# I3: SU(2) column representation V = (x)_j (conj u on even, u on odd slices) commutes with K_e, K_o
from scipy.linalg import expm
rng = np.random.default_rng(1)
Hm = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); Hm = Hm + Hm.conj().T; Hm -= np.trace(Hm) / 2 * np.eye(2)
u = expm(1j * Hm); m = 3
Wop = np.array([[1.]])
for j in range(2 * m):
    Wop = np.kron(Wop, u.conj() if j % 2 == 0 else u)
Wun = np.array([[1.]])
for j in range(2 * m):
    Wun = np.kron(Wun, u)
T = Transfer(bond_xxz_rot(1.0), 2, 2.0, m)
Ke, Ko = T.Ke.toarray(), T.Ko.toarray()
c1 = max(abs(Wop @ Ke - Ke @ Wop).max(), abs(Wop @ Ko - Ko @ Wop).max())
c2 = abs(Wun @ Ke - Ke @ Wun).max()
report('I3_su2_column_representation', c1 < 1e-12 and c2 > 1e-3, commutator_alternating=c1, commutator_uniform=c2)

# I4: two-insertion non-additivity (Heisenberg, m = 6, beta = 6): opposite-momentum triplet pair
import two_insertion as ti
S6 = Sectors('HEIS', 6.0, 6, h=bond_xxz_rot(1.0), Qmax=3)
S6, levels, pairs = ti.run(S6, 6)
opp = [p for p in pairs if p['sector'] == (0, 0) and p['a'][1] != 0]
best = max(opp, key=lambda p: p['weight_below_90pct'])
like = [p for p in pairs if p['sector'][0] == 2 and p['sector'][1] == 0][0]
report('I4_two_insertion_nonadditive', best['weight_below_90pct'] > 0.5 and like['dominant_eps'] > 1.5 * like['sum_eps'],
       opposite_pair=str((best['a'], best['b'], round(best['sum_eps'], 3), round(best['dominant_eps'], 3), round(best['weight_below_90pct'], 3))),
       like_pair=str((like['a'], like['b'], round(like['sum_eps'], 3), round(like['dominant_eps'], 3))))

# M1: m = infinity spot check of the stored small-ring table (Heisenberg, beta = 4, L = 4)
p = os.path.join(HERE, 'data', 'infinite_m.json')
if os.path.exists(p):
    stored = json.load(open(p))['HEIS']['4']
    from infinite_m import ring_logZ
    g4 = ring_logZ(bond_xxz_rot(1.0), 4, [4.0], 'u1')[0] + 4 * stored['eps0_inf']
    report('M1_infinite_m_spot_check', abs(g4 - stored['logZt']['4']) < 1e-9, recomputed=g4, stored=stored['logZt']['4'])

out['runtime_seconds'] = time.time() - t0
out['all_pass'] = bool(ok_all)
json.dump(out, open(os.path.join(HERE, 'checks.json'), 'w'), indent=1, default=str)
print('ALL PASS' if ok_all else 'SOME CHECKS FAILED', 'in %.0f s' % out['runtime_seconds'])
