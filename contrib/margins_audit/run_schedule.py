"""Margins of the leading logical chain for n = 1..23 under an explicit schedule.

Fixed for every n (C33 / C37 family, all inside the original 23 letters):
  t = t_x = t_sw = 1, stub energy a = 4 above the vacuum (Delta = 4 - 2cos(pi/(n+2))),
  D_x = D_M = 1000, delta = 0, phi^x_0 = pi/(2n), phi^x_+- = 0, port phase pi/2,
  stub and switch phases 0, leading order in g (coefficient of lambda^2 g^2 alpha^4).
Loop and field (units: beta = 1/det C_express, T in units alpha^4 g^2 lambda^2):
  t_d = beta 2^-p_n,  p_n = round(log2(kappa1_n / kappa_star)),
      kappa1_n = whitened loop norm ||W mu W^dag|| per unit t_d/beta (linear probe),
  eps = 2 beta^2 2^-f_n,  f_n = ceil(-log2 lambda_min(F) + 2),  F = Herm(-T)/(2 beta^2),
so the whitened field is between 1/8 and 1/4 of the narrowest band.  kappa_star
= 1e-3 reproduces C33's n = 23 choice p = 136.
Loop phases: 'C33' = (phi^d_0, phi^d_+, phi^d_-) = (pi, pi/2, pi/2),
             'C37' = (pi, 0, 0).
Usage: python3 -B run_schedule.py [n ...]   (default 1..23, parallel)
Writes results/n{n}.json."""
import sys, os, json, time, math
from decimal import Decimal
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hpmat                                    # noqa: E402
from leading import leading, default_params, beta, alpha4   # noqa: E402
from bands import analyse                       # noqa: E402

KAPPA_STAR = 1e-3
KAPPA_BIG = 2.0 ** -3
PROBE = Decimal(10) ** -40


def prec_for(n):
    return 6 * (n + 1) + 3 * n + 100


def S(x, d=30):
    """Decimal -> string with d significant digits"""
    if x is None:
        return None
    if isinstance(x, (int, float, str)):
        return x
    return format(x, '.%de' % (d - 1))


def summarize(o, b):
    """JSON-ready summary of an analyse() result.  Velocities in units of
    alpha^4 g^2 lambda^2 t (one cell = one lattice unit) and of 2 beta^2 that."""
    two_b2 = 2 * b * b
    bands = [dict(centre_k=S(o['centres_k'][m], 25), kF=S(o['kF'][m], 25),
                  v_left=S(o['vel'][o['bands'][m]['kL']], 12), v_right=S(o['vel'][o['bands'][m]['kR']], 12),
                  vw_left=o['vel_whitened'][o['bands'][m]['kL']], vw_right=o['vel_whitened'][o['bands'][m]['kR']])
             for m in range(o['n'])]
    rec = dict(Phi=o['Phi'], unit_circle_dev=o['unit_circle_dev'], newton_shift=o['newton_shift'],
               newton_iters=o['newton_iters'], newton_imag_ratio=o['newton_imag_ratio'],
               inertia_steps_ok=o['inertia_steps_ok'], unimodal=o['unimodal'],
               min_abs_whitened_eig_on_arcs=o['min_abs_eig_on_arcs'],
               lambda_min_Rs=S(o['lmin_Rs'], 15), lambda_max_Rs=S(o['lmax_Rs'], 15),
               log10_lambda_min_F=float((o['lmin_Rs'] / two_b2).log10()),
               log10_lambda_max_F=float((o['lmax_Rs'] / two_b2).log10()),
               chol_min_pivot=S(o['chol_min_pivot'], 10),
               norm_A_whitened=o['normA'], norm_B_whitened=o['normB'],
               loop_whitened=o['loop_whitened'], field_whitened=o['field_whitened'],
               roots_k=[S(x, 25) for x in o['roots_k']], bands=bands,
               delta2=S(o['delta2'], 12), witness2=o['w2'],
               delta4=S(o['delta4'], 12), witness4=o['w4'],
               dPT=S(o['dPT'], 12), witnessPT=o['wPT'], dU1=S(o['dU1'], 12),
               dU2=S(o['dU2'], 12), dF=S(o['dF'], 12),
               small_gaps=[(S(g, 8), w1, w2) for g, w1, w2 in o['small_gaps']],
               log10_vmin=float(o['vmin'].log10()), log10_vmax=float(o['vmax'].log10()),
               log10_vmin_over_2beta2=float((o['vmin'] / two_b2).log10()),
               vmin_whitened=o['vw_min'])
    return rec


def run_n(n, prec=None, do_c37=True, do_big=True):
    t0 = time.time()
    prec = prec or prec_for(n)
    hpmat.setprec(prec)
    P = default_params(n)
    b = beta(n, P['a'], P['Dx'])
    out = dict(n=n, prec=prec, log10_beta=float(b.log10()), log10_alpha4=float(alpha4(n).log10()),
               kappa_star=KAPPA_STAR)
    # linear loop probe, also gives lambda_min of the (essentially tree) hopping
    Pp = dict(P, td=b * PROBE)
    L = leading(n, Pp)
    o, wh = analyse(L['T'], L['mu'], Decimal(0))
    lminF = o['lmin_Rs'] / (2 * b * b)
    kappa1 = o['loop_whitened'] / float(PROBE)
    out['probe'] = dict(td_over_beta='1e-40', kappa1=kappa1, log10_lambda_min_F=float(lminF.log10()),
                        log10_lambda_max_F=float((o['lmax_Rs'] / (2 * b * b)).log10()),
                        herm_err=float(L['herm_err']), eig=L['eig'])
    f = math.ceil(-float(lminF.ln() / Decimal(2).ln()) + 2)
    eps = 2 * b * b * Decimal(2) ** (-f)
    out['f'] = f
    runs = {}
    plist = [('C33_kappa1e-3', 'C33', KAPPA_STAR)]
    if do_big:
        plist.append(('C33_kappa1/8', 'C33', KAPPA_BIG))
    for label, phases, kst in plist:
        if n == 1 or kappa1 == 0:
            p = None
            td = Decimal(0)
        else:
            p = int(round(math.log2(kappa1 / kst)))
            td = b * Decimal(2) ** (-p)
        Pr = dict(P, td=td)
        L = leading(n, Pr)
        o, wh = analyse(L['T'], L['mu'], eps)
        rec = summarize(o, b)
        rec.update(p=p, f=f, herm_err=float(L['herm_err']))
        runs[label] = rec
        print(json.dumps(dict(n=n, run=label, p=p, f=f, d2=rec['delta2'], d4=rec['delta4'], dPT=rec['dPT'],
                              dU1=rec['dU1'], dU2=rec['dU2'], dF=rec['dF'], unimodal=rec['unimodal'],
                              loopW=rec['loop_whitened'], t=round(time.time() - t0))), flush=True)
    if do_c37 and n > 1:
        # C37 loop phases (pi, 0, 0) in the C37 hierarchy, loop << field, both small
        # but resolvable: whitened loop about 1e-25 (scaled with the C33 kappa1),
        # whitened field 1e-4 of the schedule value.  Margins split by the loop scale
        # linearly in t_d, those split by the field linearly in eps, so both can be
        # carried to C37's literal t_d = 1e-2000, e = 1e-1000.
        tauc = Decimal(1e-25) / Decimal(kappa1)
        Pc = dict(P, td=b * tauc, phd={'0': 'pi', '+': '0', '-': '0'})
        L = leading(n, Pc)
        epsc = eps * Decimal(10) ** -4
        o, wh = analyse(L['T'], L['mu'], epsc, rdigits=70)
        rec = summarize(o, b)
        rec.update(td_over_beta=S(tauc, 6), eps=S(epsc, 6), f=f,
                   small_gaps=[(S(g, 8), w1, w2) for g, w1, w2 in o['small_gaps']])
        # extrapolation to C37's literal point, in leading.py units
        lmin = o['lmin_Rs']
        loopw_c37 = Decimal(rec['loop_whitened']) / tauc * (Decimal(10) ** -2000 / b)
        fieldw_c37 = Decimal(10) ** -1000 / lmin
        rec['C37_point'] = dict(loop_whitened=S(loopw_c37, 6), field_whitened=S(fieldw_c37, 6),
                                log10_dPT=float((o['dPT'] * loopw_c37 / Decimal(rec['loop_whitened'])).log10()),
                                log10_dU2=float((o['dU2'] * fieldw_c37 / Decimal(rec['field_whitened'])).log10()),
                                log10_dU1=float((o['dU1'] * fieldw_c37 / Decimal(rec['field_whitened'])).log10())
                                if o['dU1'] < Decimal('1e-2') else None)
        runs['C37_phases_linear'] = rec
        print(json.dumps(dict(n=n, run='C37_phases_linear', d2=rec['delta2'], d4=rec['delta4'], dPT=rec['dPT'],
                              dU1=rec['dU1'], dU2=rec['dU2'], loopW=rec['loop_whitened'],
                              fieldW=rec['field_whitened'], C37=rec['C37_point'],
                              t=round(time.time() - t0))), flush=True)
    if n == 23:
        # C33's certified point exactly: t_d = beta 2^-136, eps = 2 beta^2 2^-214
        Pr = dict(P, td=b * Decimal(2) ** -136)
        L = leading(n, Pr)
        o, wh = analyse(L['T'], L['mu'], 2 * b * b * Decimal(2) ** -214)
        rec = summarize(o, b)
        rec.update(p=136, f=214)
        runs['C33_exact_n23'] = rec
        print(json.dumps(dict(n=n, run='C33_exact_n23', d2=rec['delta2'], d4=rec['delta4'], t=round(time.time() - t0))), flush=True)
    out['runs'] = runs
    out['seconds'] = time.time() - t0
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'results', 'n%d.json' % n), 'w'), indent=1)
    return n, out['seconds']


if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1:]] or list(range(1, 24))
    ns = sorted(ns, reverse=True)
    with Pool(min(9, len(ns))) as pool:
        for n, sec in pool.imap_unordered(run_n, ns):
            print('done', n, round(sec), flush=True)
