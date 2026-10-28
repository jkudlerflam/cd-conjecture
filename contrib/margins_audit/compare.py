"""Comparison with Astra C33 and C37 where they overlap.  Reads their files and
imports their float code read-only (no bytecode is written).
  A. code level: leading.py vs C33 and C37 structured.Leading at n = 2, 3, 4
     (D = 10, a = 4, t_d = 0.01, their respective loop phases);
  B. positivity brackets: log10 lambda_min, lambda_max of F = Herm(-T)/(2 beta^2)
     at D = 1000 vs C33 certificate_n*.json and C37 positivity_n*.json;
  C. C33's certified n = 23 finite-loop point vs results/n23.json;
  D. C37: margins at C37's loop phases, carried to its literal t_d, field.
Usage: python3 -B compare.py   Writes compare.json."""
import sys, os, json, glob, importlib.util
sys.dont_write_bytecode = True
import numpy as np
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
C33 = os.path.join(ROOT, 'codex', 'work_C33')
C37 = os.path.join(ROOT, 'codex', 'work_C37')
import hpmat                                     # noqa: E402
from leading import leading, default_params, alpha4   # noqa: E402


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def part_a():
    s33 = load(os.path.join(C33, 'structured.py'), 's33')
    s37 = load(os.path.join(C37, 'structured.py'), 's37')
    out = []
    for n in (2, 3, 4):
        hpmat.setprec(60)
        for label, phd, ref in [('C33', {'0': 'pi', '+': 'pi/2', '-': 'pi/2'},
                                 lambda: s33.Leading(n, D=10., a=4., td=0.01, loopflux=np.pi,
                                                     chargedloopflux=np.pi / 2).calculate()),
                                ('C37', {'0': 'pi', '+': '0', '-': '0'},
                                 lambda: s37.Leading(n, D=10., a=4., td=0.01, loopflux=np.pi).calculate())]:
            P = default_params(n)
            P.update(Dx=Decimal(10), DM=Decimal(10), td=Decimal('0.01'), phd=phd)
            L = leading(n, P)
            a4 = float(alpha4(n))
            T, mu = a4 * L['T'].tocomplex(), a4 * L['mu'].tocomplex()
            Tr, Mr, _ = ref()
            sc = np.abs(T).max()
            rec = dict(n=n, code=label, T_rel=float(np.abs(T - Tr).max() / sc),
                       mu_rel_to_T=float(np.abs(mu - Mr).max() / sc), mu_over_T=float(np.abs(mu).max() / sc))
            print(json.dumps(rec), flush=True)
            out.append(rec)
    return out


def part_b():
    out = []
    for n in range(1, 24):
        f = os.path.join(HERE, 'results', 'n%d.json' % n)
        if not os.path.exists(f):
            continue
        R = json.load(open(f))
        c = os.path.join(C33, 'certificate_n%d.json' % n)
        src = 'C33'
        if not os.path.exists(c):
            c = os.path.join(C37, 'positivity_n%d.json' % n)
            src = 'C37'
        if not os.path.exists(c):
            continue
        C = json.load(open(c))
        lmin, lmax = R['probe']['log10_lambda_min_F'], R['probe']['log10_lambda_max_F']
        rec = dict(n=n, source=src, mine_lmin=lmin, cert_lmin=[C['log10_lambda_min_lower'], C['log10_lambda_min_upper']],
                   inside_min=C['log10_lambda_min_lower'] <= lmin + 1e-9 and lmin <= C['log10_lambda_min_upper'] + 1e-9,
                   mine_lmax=lmax, cert_lmax=[C['log10_lambda_max_lower'], C['log10_lambda_max_upper']],
                   inside_max=C['log10_lambda_max_lower'] <= lmax + 1e-9 and lmax <= C['log10_lambda_max_upper'] + 1e-9,
                   mine_log10_beta=R['log10_beta'], cert_log10_beta=C['log10_beta'])
        print(json.dumps(rec), flush=True)
        out.append(rec)
    return out


def part_c():
    f = os.path.join(HERE, 'results', 'n23.json')
    if not os.path.exists(f):
        return None
    R = json.load(open(f))
    r = R['runs'].get('C33_exact_n23')
    if r is None:
        return None
    B = json.load(open(os.path.join(C33, 'full_loop_bands_n23.json')))
    mine = sorted(float(x) for x in r['roots_k'])
    theirs = sorted(B['roots'])
    mirror = sorted(-x for x in theirs)
    d_same = max(abs(a - b) for a, b in zip(mine, theirs))
    d_mirror = max(abs(a - b) for a, b in zip(mine, mirror))
    mc = sorted(float(b['centre_k']) for b in r['bands'])
    tc = sorted(B['centres'])
    rec = dict(max_root_diff=d_same, max_root_diff_if_mirrored=d_mirror,
               mine_delta2=float(r['delta2']), C33_delta2=B['delta2'], C33_cert_delta2_lower=0.0140193,
               mine_delta4=float(r['delta4']), C33_delta4=B['delta4'], C33_cert_delta4_lower=3.1558e-7,
               mine_witness4=r['witness4'], C33_witness4=B['witness'],
               max_centre_diff_sorted=max(abs(a - b) for a, b in zip(mc, tc)),
               mine_log10_vmin_2beta2=r['log10_vmin_over_2beta2'], C33_log10_min_speed=B['log10_min_speed_normalized'])
    print(json.dumps(rec), flush=True)
    return rec


def part_d():
    out = []
    for n in range(2, 24):
        f = os.path.join(HERE, 'results', 'n%d.json' % n)
        if not os.path.exists(f):
            continue
        r = json.load(open(f))['runs'].get('C37_phases_linear')
        if r is None:
            continue
        out.append(dict(n=n, loop_whitened_run=r['loop_whitened'], field_whitened_run=r['field_whitened'],
                        dPT_run=float(r['dPT']), dU2_run=float(r['dU2']), dU1_run=float(r['dU1']),
                        delta2_run=float(r['delta2']), **r['C37_point']))
    for x in out:
        print(json.dumps(x), flush=True)
    return out


if __name__ == '__main__':
    res = dict(A=part_a(), B=part_b(), C=part_c(), D=part_d())
    json.dump(res, open(os.path.join(HERE, 'compare.json'), 'w'), indent=1)
