# C193b diagnostics (unregistered, written after the candidate outputs were opened): test whether the qubit
# chain's open ends pin the charge phase theta_c (Dirichlet), which the effective model's ends do not.
import json, sys, numpy as np
import ana193 as A

def mixed_fit(dat, lo, hi):
    """ln|h| = c - eta ln tan u - 2 eta_c ln cos u on symmetric pairs, u = pi r / (2(L+1)); both valleys."""
    L = dat['L']
    rs, g, F = A.onebody_seq(dat['G'], L)
    out = []
    for k, amp in A.peaks(rs, g, 2)[:2]:
        hh = A.sm(g * np.exp(-1j * k * rs)); s = (rs >= 6) & (rs <= rs[-1] - A.NSM)
        k = k + np.polyfit(rs[s], np.unwrap(np.angle(hh[s])), 1)[0]
        h = np.abs(A.sm(g * np.exp(-1j * k * rs)))
        u = np.pi * rs / (2 * (L + 1)); s = (rs >= lo * L) & (rs <= min(hi * L, rs[-1] - A.NSM))
        X = np.stack([np.ones(s.sum()), -np.log(np.tan(u[s])), -2 * np.log(np.cos(u[s]))], 1)
        c, eta, etac = np.linalg.lstsq(X, np.log(h[s]), rcond=None)[0]
        res = np.log(h[s]) - X @ np.array([c, eta, etac])
        out.append((eta, etac, float(np.std(res))))
    return out

def pp_modes(fn):
    m = json.load(open(fn))['meas']; L = m['L']; res = {}
    for d in ('2', '3'):
        pp = m['pp'][d]; js = np.array(sorted(int(j) for j in pp)); v = np.array([complex(*pp[str(j)]) for j in js])
        s = (js >= L // 2) & (js <= 3 * L // 4)            # one quarter, away from the right edge
        a0 = abs(np.mean(v[s])); api = abs(np.mean(v[s] * (-1.0) ** js[s]))
        res[d] = (a0, api)
    return L, res

if __name__ == '__main__':
    print('mixed boundary fit of the one-body function: eta (bulk), eta_c (charge part with Dirichlet theta_c), rms')
    for t in sys.argv[1:]:
        dat = A.from_eff(t) if not t.endswith('rand1') and '_P' not in t else A.from_qubit(f'typhon_out/out/{t}_meas.json')
        rr = [(lo, hi, mixed_fit(dat, lo, hi)) for lo, hi in ((0.05, 0.40), (0.05, 0.60), (0.10, 0.80), (0.05, 0.85))]
        print(f'{t:26s} ' + '  '.join(f'[{lo},{hi}] ' + ' '.join(f'{e:.4f}/{ec:+.4f}({r:.3f})' for e, ec, r in o) for lo, hi, o in rr))


def table():
    """Unregistered: eta from the mixed-boundary fit on [0.10, 0.80] L (mean of valleys, error = half spread over
    [0.05,0.60], [0.10,0.80], [0.05,0.85]), 1/K_c from the frozen pair exponent, ground branch."""
    an = json.load(open('analysis.json'))
    cand = an['cand']
    print('key          L   eta_mixed(err)  eta_c   pair3   Ks_mixed(err)')
    for key in A.KEYS:
        Ks = []
        for L in (64, 128, 256):
            g = cand[f'{key}_L{L}']; o = g['sec'][str(g['ground'])]
            dat = A.from_qubit('typhon_out/out/' + o['file'])
            vals = [mixed_fit(dat, lo, hi) for lo, hi in ((0.05, 0.60), (0.10, 0.80), (0.05, 0.85))]
            ets = [np.mean([e for e, _, _ in v]) for v in vals]
            etc = np.mean([ec for _, ec, _ in vals[1]])
            eta, de = ets[1], (max(ets) - min(ets)) / 2 + 0.005
            x = 4 * eta - o['eta_p']; K = 1 / x; dK = np.hypot(4 * de, o['deta_p']) / x ** 2
            Ks.append((K, dK))
            print(f'{key:11s} {L:4d}  {eta:.4f}({de:.4f})  {etc:+.3f}  {o["eta_p"]:.4f}  {K:.3f}({dK:.3f})')
        b1 = A.bound(Ks[0][0], Ks[1][0], Ks[0][1], Ks[1][1]); b2 = A.bound(Ks[1][0], Ks[2][0], Ks[1][1], Ks[2][1])
        print(f'   d(64->128)={b1["dK"]:+.3f}({b1["sig"]:.3f})  d(128->256)={b2["dK"]:+.3f}({b2["sig"]:.3f})  |y1| {b2["y1c"]:.3f}/{b2["y1u"]:.3f}  y_par {b2["ypc"]:.3f}/{b2["ypl"]:.3f}')
