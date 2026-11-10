"""C189 Part 1 on the candidate series (registered estimators, RESULT.md amendment A2).
Usage: python3 kest_run.py > kest_run.out"""
import numpy as np, warnings, json
from kest import *
warnings.simplefilter('ignore')
D = '../work_C176/typhon_out/out/'
SER = {'sym316': (3, 16, (64, 96, 128, 192, 256)), 'sym532': (5, 32, (64, 128, 192, 256)), 'sym18': (1, 8, (64, 96, 128, 192, 256))}


def Kc_est(run, ffr):
    v = [K_from_Sq(run, lo, hi, 2)[0] for lo, hi in ((0.05, 0.3), (0.05, 0.2), (0.1, 0.3))]
    bias = abs(K_from_Sq(ffr, 0.05, 0.3, 2)[0] - 1.0)
    return v[0], np.hypot(0.5 * (max(v) - min(v)), bias), v


def F_est(run, ffr):
    L = run['L']; nb = n_mid(run)
    ks = [np.pi * nb, 2 * K0]
    v = [friedel_fit(run, lo, hi, ks)['eta'] for lo, hi in ((8, L // 4), (8, L // 2), (12, L // 4))]
    fb = abs(friedel_fit(ffr, 8, L // 4, [np.pi * n_mid(ffr), 2 * K0])['eta'] - 1.0)
    return v[0], np.hypot(0.5 * (max(v) - min(v)), fb), v


def D_est(run, ffr):
    L = run['L']
    v = [eta(run, 'D', lo, hi) for lo, hi in ((16, L // 4), (24, L // 4), (16, 3 * L // 8))]
    b = eta(ffr, 'D', 16, L // 4) - 2.0
    return v[0], np.hypot(0.5 * (max(v) - min(v)), abs(b)), v, b


def branch(run):
    """Fourier weight of the wall oscillation at 2kF against 2(k0 +- kF) and 2k0 (first quarter of the chain)."""
    L = run['L']; dn = np.asarray(run['dens']); nb = n_mid(run)
    prof = 0.5 * (dn + dn[::-1]); x = np.arange(L)
    osc = (prof - nb) * np.sin(np.pi * (x + 1) / (L + 1))
    sel = (x >= 4) & (x < L // 4)
    def ft(k):
        return float(abs(np.sum(np.exp(1j * k * x[sel]) * osc[sel])))
    kk = {'2kF': np.pi * nb, '2(k0+kF)': ks_two_valley(nb)[0], '2(k0-kF)': ks_two_valley(nb)[1], '2k0': 2 * K0}
    return {k: ft(v) for k, v in kk.items()}


print('Part 1 estimators, candidate series (phi = pi/2, |w| = 1, bare couplings), C176 final sweeps')
res = {}
for s, (p, q, Ls) in SER.items():
    print(f'\n=== {s} (N/L = {p}/{q})')
    for L in Ls:
        try:
            run = load_run(D + f'{s}_L{L}.jsonl')
        except FileNotFoundError:
            continue
        N = run['N']; ffr = ff_run(L, N)
        nb = n_mid(run)
        Kc, sKc, vKc = Kc_est(run, ffr)
        eF, seF, veF = F_est(run, ffr)
        KsF = 2 * eF - Kc; sKsF = np.hypot(2 * seF, sKc)
        line = (f'L={L} N={N} n_mid={nb:.4f} trunc={run["trunc"]:.1e}  K_c={Kc:.3f}+-{sKc:.3f} {np.round(vKc, 3)}  '
                f'eta_F={eF:.3f}+-{seF:.3f} {np.round(veF, 3)}  K_s^F={KsF:.3f}+-{sKsF:.3f}')
        rec = dict(L=L, N=N, n_mid=nb, Kc=Kc, sKc=sKc, etaF=eF, setaF=seF, KsF=KsF, sKsF=sKsF)
        if L >= 128:
            eD, seD, veD, bD = D_est(run, ffr)
            KsD = eD - 1 / Kc; sKsD = np.hypot(seD, sKc / Kc ** 2)
            KsDc = (eD - bD) - 1 / Kc
            eG = eta(run, 'G', 16, 64)
            ssum = 4 * eG - Kc - 1 / Kc
            line += (f'\n      eta_D={eD:.3f}+-{seD:.3f} {np.round(veD, 3)} (ff2v bias {bD:+.3f})  K_s^D={KsD:.3f}+-{sKsD:.3f}  '
                     f'bias-corrected K_s^D={KsDc:.3f}  [uncalibrated: eta_G={eG:.3f} -> K_s+1/K_s={ssum:.2f}]')
            rec.update(etaD=eD, setaD=seD, KsD=KsD, sKsD=sKsD, KsDc=KsDc, biasD=bD, etaG=eG)
        br = branch(run)
        line += '\n      wall Fourier weights: ' + ', '.join(f'{k} {v:.4f}' for k, v in br.items())
        rec['branch'] = br
        print(line)
        res[f'{s}_L{L}'] = rec
json.dump(res, open('kest_run.json', 'w'), indent=1)
