"""Second-order defect coefficients on rings: python3 -B run.py N V [procs]. Writes out/N<N>_V<V>.json."""
import sys, json, time, os
from pathlib import Path
import numpy as np
from kernel import kernel
from sources import selector, frozen, B_coeff
HERE = Path(__file__).resolve().parent
b = 15


def patterns(k):
    g = [b] * k
    P = {'good': (g, [False] * k)}
    P['pair_14_16_adjacent'] = ([14, 16] + g[2:], [False] * k)
    if k >= 4:
        h = k // 2
        P['pair_14_16_far'] = ([14] + g[1:h] + [16] + g[h + 1:], [False] * k)
    P['pair_13_17_adjacent'] = ([13, 17] + g[2:], [False] * k)
    P['long_30'] = ([30] + g[2:], [False] * (k - 1))
    P['short_7_8'] = ([7, 8] + g[1:], [False] * (k + 1))
    P['pin'] = (g, [True] + [False] * (k - 1))
    return P


def W_of(lengths, pins):
    return sum(l for l in lengths if l != b) + sum(1 for l, p in zip(lengths, pins) if l == b and p)


def main(N, V, procs):
    t0 = time.time()
    chi, diag = kernel(N, V, procs=procs)
    diag['seconds'] = time.time() - t0
    mu = selector(b); a = diag['a']
    rows = {}
    Bg = None
    for name, (L, p) in patterns(N // b).items():
        B = B_coeff(chi, L, p, mu)
        if name == 'good':
            Bg = B; Eg = frozen(L, p, mu); continue
        dE = frozen(L, p, mu) - Eg
        W = W_of(L, p)
        dB = B - Bg
        rows[name] = dict(lengths=L, pins=p, W=W, dE_material=dE, first_order_over_J=a * dE,
                          delta_B=dB, delta_B_over_W=dB / W,
                          relative_coeff=-dB / (a * dE))   # cost = J a dE [1 + (J/t) relative_coeff + ...]
    out = dict(N=N, V=V, L=N + V, density=N / (N + V), b=b, mu_sel=mu, B_good=Bg, diagnostics=diag,
               chi=chi.tolist(), defects=rows)
    (HERE / 'out').mkdir(exist_ok=True)
    (HERE / 'out' / f'N{N}_V{V}.json').write_text(json.dumps(out, indent=1) + '\n')
    print(N, V, 'sec', round(diag['seconds']), 'sum_rule', diag['sum_rule_error'], 'imag', diag['max_imag'],
          'tail', diag['tail_integrand'], flush=True)
    for k, r in rows.items():
        print('  ', k, 'W', r['W'], 'aDE', '%.4e' % r['first_order_over_J'], 'dB', '%.6e' % r['delta_B'],
              'dB/W', '%.4e' % r['delta_B_over_W'], 'rel', '%.4e' % r['relative_coeff'], flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1)
