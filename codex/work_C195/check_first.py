"""C195 first-sweep check (Claude Opus 5.5 agent, 2026-11-10).  numpy only.

Usage: python3 check_first.py <dir with the L = 512 run logs (out512/logs)> [<C190 outG folder>]
For every logged sweep of the L = 512 runs: residual against its target, energy, density, and the comparison with
the L = 256 and L = 128 runs of p153c (rand:1) at the same sweep index (the same stage of the schedule up to
chi = 512): the bulk energy per site (E(512) - E(256))/256 against (E(256) - E(128))/128 where L = 128 is at the
same chi, E/L at both sizes, and n_mid against the L = 256 value.
"""
import os
import sys
import json

d512 = sys.argv[1]
dG = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'work_C190',
                                                         'gpu', 'typhon_out', 'outG')
for P in (0, 1):
    recs = []
    for line in open(os.path.join(d512, f'p153c_L512_P{P}_rand1.log')):
        if line.startswith('{'):
            recs.append(json.loads(line))
    h256 = json.load(open(os.path.join(dG, f'p153c_L256_P{P}_rand1.json')))['hist']
    h128 = json.load(open(os.path.join(dG, f'p153c_L128_P{P}_rand1.json')))['hist']
    print(f'== P = {P}: {len(recs)} sweeps logged')
    for i, r in enumerate(recs):
        ok = r['res_max'] <= r['res_tol'] * (1 + 1e-9)
        s = (f"{i:2d} chi {r['chi']:4d} a {r['alpha']:g} res {r['res_max']:.2e}/{r['res_tol']:.0e} {'met' if ok else 'MISSED'} "
             f"mv {r['mv_mean']:.1f}/{r['mv_max']} trunc {r['trunc']:.1e} t {r['t']:.0f}s E {r['E_phys']:.8f} E/L {r['E_phys'] / 512:.9f} "
             f"n_mid {r['n_mid']:.4f} N {r['N_tot']:.2f}")
        if r.get('dev_used_GB') is not None:
            s += f" dev {r['dev_used_GB']:.1f}GB"
        if r.get('host_rss_GB') is not None:
            s += f" rss {r['host_rss_GB']:.1f}GB"
        if i < len(h256) and h256[i]['chi'] == r['chi']:
            a = h256[i]
            s += f" | L256: E/L {a['E_phys'] / 256:.9f} n_mid {a['n_mid']:.4f} (dn {r['n_mid'] - a['n_mid']:+.4f})"
            bulk512 = (r['E_phys'] - a['E_phys']) / 256
            s += f" bulk(512-256) {bulk512:.9f}"
            if i < len(h128) and h128[i]['chi'] == r['chi']:
                b = h128[i]
                s += f" bulk(256-128) {(a['E_phys'] - b['E_phys']) / 128:.9f} diff {bulk512 - (a['E_phys'] - b['E_phys']) / 128:+.1e}"
        if 'Sbar6' in r:
            s += f" S6 {r['Sbar6']:.4f}"
        print(s)
