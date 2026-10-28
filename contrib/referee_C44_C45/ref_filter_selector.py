"""Independent checks of C44 PROOFS Sec. 3: occupation filter and selector.

Filter (untwisted frame): r qubits, ferromagnetic singlet projectors on every
bond, left pin 1-|+><+| on site 1, right pin 1-|u><u| on site r with
|u> = (|0> + e^{i alpha}|1>)/sqrt2, alpha = (r-k) 2pi/(k+2).
Claims tested: left-pinned gap >= kappa_r = 1/[r(2r-1)];
F_r >= kappa_r beta_r/[(1+sqrt kappa_r)^2+beta_r] for rejected r;
J_f F_r >= 2U with J_f = 2U/eps_*; zero energy iff r = k (2 <= r <= k+2).

Selector: a_l = 2(1-cos(pi/l)) = 4 sin^2(pi/(2l)). Tests mu_c>0, eta_b>0,
e<0, mu_sel>0, the secant inequality (1) for every b up to 4470 (n up to
4467, which covers n = 4463, the first n beating d0-1 = 4463), and
reproduces the four quoted values of t_sel eta_b - U.
"""
import json, math
import numpy as np

out = {}

# ---------------- occupation filter ----------------
def filt(r, k, right_pin=True):
    I2 = np.eye(2)
    swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float)
    sing = (np.eye(4) - swap) / 2
    plus = np.array([1, 1]) / math.sqrt(2)
    H = np.zeros((2 ** r, 2 ** r), complex)
    for j in range(r - 1):
        H += np.kron(np.kron(np.eye(2 ** j), sing), np.eye(2 ** (r - j - 2)))
    H += np.kron(I2 - np.outer(plus, plus), np.eye(2 ** (r - 1)))
    if right_pin:
        al = (r - k) * 2 * math.pi / (k + 2)
        u = np.array([1, np.exp(1j * al)]) / math.sqrt(2)
        H += np.kron(np.eye(2 ** (r - 1)), I2 - np.outer(u, u.conj()))
    return np.linalg.eigvalsh(H)

rows, ok = [], True
for k in range(3, 9):                       # n = 1..6
    kap_s = 1 / ((k + 2) * (2 * k + 3))
    bet_s = math.sin(math.pi / (k + 2)) ** 2
    eps_s = kap_s * bet_s / ((1 + math.sqrt(kap_s)) ** 2 + bet_s)
    for r in range(2, k + 3):
        kap = 1 / (r * (2 * r - 1))
        gl = filt(r, k, right_pin=False)
        left_gap = gl[1] - gl[0]
        F = filt(r, k)[0]
        beta = math.sin((r - k) * math.pi / (k + 2)) ** 2
        bnd = kap * beta / ((1 + math.sqrt(kap)) ** 2 + beta)
        row = dict(k=k, r=r, left_gap=left_gap, kappa=kap, F=F, bound=bnd,
                   JfF_over_2U=F / eps_s)
        good = left_gap >= kap - 1e-12 and F >= bnd - 1e-12
        good = good and ((r == k and abs(F) < 1e-10) or (r != k and F / eps_s >= 1 - 1e-12))
        row['ok'] = bool(good)
        ok &= good
        rows.append(row)
out['filter_ok'] = bool(ok)
out['filter_min_left_gap_over_kappa'] = min(r['left_gap'] / r['kappa'] for r in rows)
out['filter_min_F_over_bound'] = min(r['F'] / r['bound'] for r in rows if r['r'] != r['k'])
out['filter_min_JfF_over_2U'] = min(r['JfF_over_2U'] for r in rows if r['r'] != r['k'])

# ---------------- selector ----------------
def a(l):
    return 4 * np.sin(np.pi / (2 * np.asarray(l, float))) ** 2

def selector(b, U=1.0, e1=None, e2=None):
    ab, ab1 = a(b), a(b + 1)
    delta_b = (b + 1) * ab - b * ab1
    mu_c = 2 - delta_b
    e_c = ab1 - ab
    eta = min((a(b - 1) - 2 * ab + ab1) / 2, (ab - 2 * ab1 + a(b + 2)) / 3)
    t_sel = 2 * U * (1 / eta + b / mu_c + 1 / abs(e_c))
    return dict(mu_c=mu_c, e_c=e_c, eta=eta, t_sel=t_sel, coeff=t_sel * eta - U)

quoted = [1.69357336, 1.39840043, 1.30289697, 1.25125262]
out['selector_quoted_check'] = [dict(n=n, ours=float(selector(n + 3)['coeff']), quoted=q)
                                for n, q in zip(range(1, 5), quoted)]
worst = dict(min_mu_c=1e9, min_eta_b4=1e9, min_secant_slack=1e9)
for b in range(4, 4471):
    s = selector(b)
    worst['min_mu_c'] = min(worst['min_mu_c'], float(s['mu_c']))
    worst['min_eta_b4'] = min(worst['min_eta_b4'], float(s['eta'] * b ** 4))
    # secant inequality (1) on l in [2, 6b], selected lengths excluded
    l = np.arange(2, 6 * b + 1)
    l = l[(l != b) & (l != b + 1)]
    lhs = a(l) - a(b) - (l - b) * (a(b + 1) - a(b))
    rhs = s['eta'] * (1 + np.abs(l - b))
    slack = float(np.min((lhs - rhs) / np.maximum(rhs, 1e-300)))
    worst['min_secant_slack'] = min(worst['min_secant_slack'], slack)
    # e<0 and mu_sel>0 hold for any 0<=e_1,e_2<=U given t_sel |e_c| >= 2U, t_sel mu_c >= 2Ub
    assert s['t_sel'] * abs(s['e_c']) >= 2 - 1e-9 and s['t_sel'] * s['mu_c'] >= 2 * b - 1e-9
out['selector_b_range'] = [4, 4470]
out['selector_worst'] = worst
out['selector_b4466'] = {k: float(v) for k, v in selector(4466).items()}
json.dump(dict(out, filter_rows=rows), open('ref_filter_selector.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
