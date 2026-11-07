"""Collect the numerical results into the tables quoted in RESULT.md and PROOF.md (prints markdown)."""
import os, json, glob
import numpy as np
from analysis_majorant import table, LS
from majorant import vmin_needed
from models import MODELS, get

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ['XX', 'XXZ05', 'HEIS', 'XXZ2', 'TFIM', 'XYI', 'XYZ', 'XXZSTAG']


def zero_T():
    z = json.load(open(os.path.join(HERE, 'data', 'zeroT.json')))
    print('| model | M | A/M | D(2)/M | D(18)/M | pi A/(6M) | V0/M, S=1 | V0/M, S=2 | V0/M, S=6 |')
    print('| :-- | --: | --: | --: | --: | --: | --: | --: | --: |')
    for n in NAMES:
        r = z[n]; M = r['M']; D2 = r['D'][0] / M
        v = [12 * D2 / (np.pi * S) for S in (1, 2, 6)]
        print('| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %.3f | %.3f |' % (n, M, r['A'] / M, D2, r['D'][-1] / M, r['cft_value_over_M'], *v))


def finite_m(dmax=1.0, Ks=(np.log(2), 2.0), S=2):
    print('| model | m | beta range | ' + ' | '.join("V'/M (S=%d, K=%.2f)" % (S, K) for K in Ks) + ' | kappa(K=ln2) | max log Z~ |')
    print('| :-- | --: | :-- | ' + ' | '.join('--:' for _ in Ks) + ' | --: | --: |')
    for n in NAMES:
        M = get(n)[2]
        for m in range(2, 10):
            t = table(n, m)
            if t is None:
                continue
            betas, g = t
            sel = betas / m <= dmax
            if not sel.any():
                continue
            LL, BB = np.meshgrid(LS, betas[sel], indexing='ij')
            gg = g[:, sel]
            vals = [vmin_needed(gg, LL, BB, S, K).max() / M for K in Ks]
            kap = ((LL / BB) * (gg - np.log(2))).max() / M
            print('| %s | %d | %g..%g | %s | %.3f | %.3f |' % (n, m, betas[sel].min(), betas[sel].max(), ' | '.join('%.3f' % v for v in vals), kap, gg.max()))


def insertion_summary():
    files = sorted(glob.glob(os.path.join(HERE, 'data', 'insertion_m*_b*.json')))
    print('| model | m | beta | fmin | log Z~(2) | log(1+Z1)(2) | log Z~(4) | log(1+Z1)(4) | log Z~(8) | log(1+Z1)(8) | max_L [log Z~ - log Z_F] | max_L [log Z~ - log Z_F3] | first count violation E/eps1min |')
    print('| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |')
    rows = []
    for f in files:
        for r in json.load(open(f)):
            rows.append(r)
    rows.sort(key=lambda r: (NAMES.index(r['Delta']) if r['Delta'] in NAMES else 99, r['m'], r['beta']))
    for r in rows:
        i2, i4, i8 = r['L'].index(2), r['L'].index(4), r['L'].index(8)
        dF = max(a - b for a, b in zip(r['logZt'], r['logZF']))
        dF3 = max(a - b for a, b in zip(r['logZt'], r['logZF_one_per_sector'])) if 'logZF_one_per_sector' in r else float('nan')
        fv = r['first_violation_E'] / r['lowest_e1'] if r['first_violation_E'] else float('nan')
        print('| %s | %d | %g | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f |' % (
            r['Delta'], r['m'], r['beta'], r['fmin'], r['logZt'][i2], r['logZ1'][i2], r['logZt'][i4], r['logZ1'][i4],
            r['logZt'][i8], r['logZ1'][i8], dF, dF3, fv))


if __name__ == '__main__':
    print('## zero temperature'); zero_T()
    print('\n## finite m, beta/m <= 1'); finite_m()
    print('\n## insertion'); insertion_summary()


def insertion_compact(models=NAMES):
    """one row per (model, largest m, beta): fmin, single-insertion share of log Z~, Fock margins, lambda"""
    files = sorted(glob.glob(os.path.join(HERE, 'data', 'insertion_m*_b*.json')))
    best = {}
    for f in files:
        for r in json.load(open(f)):
            if 'lambda_global' not in r:
                continue
            k = (r['Delta'], r['beta'])
            if k not in best or r['m'] > best[k]['m']:
                best[k] = r
    print('| model | m | beta | fmin | share L=2 | share L=4 | share L=8 | min_L(log Z_F - log Z~) over L=2..8 | same, Z_F3 | lambda_global | lambda_sector |')
    print('| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |')
    for n in models:
        for b in (2.0, 4.0, 6.0, 8.0):
            r = best.get((n, b))
            if r is None:
                continue
            L = r['L']; sh = [r['logZ1'][L.index(x)] / r['logZt'][L.index(x)] for x in (2, 4, 8)]
            idx = [L.index(x) for x in (2, 4, 6, 8)]
            mF = min(r['logZF'][i] - r['logZt'][i] for i in idx)
            mF3 = min(r['logZF_one_per_sector'][i] - r['logZt'][i] for i in idx)
            print('| %s | %d | %g | %.3f | %.2f | %.2f | %.2f | %.3f | %.3f | %.3f | %.3f |' % (
                n, r['m'], b, r['fmin'], *sh, mF, mF3, r['lambda_global'], r['lambda_sector']))
