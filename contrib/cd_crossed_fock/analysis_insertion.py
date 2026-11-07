"""Single-insertion sector and Fock tests at finite Trotter number.
Usage: python3 analysis_insertion.py MLIST BETALIST   (comma separated), writes data/insertion_m*_b*.json

For each model (Delta = 0, 0.5, 1, 2), m and beta it records
 - fmin: smallest overlap |P_W e|^2 of the chosen one-particle eigenstates with the single-insertion span W,
 - weight1: total spectral weight of single insertions carried by the chosen one-particle levels,
 - logZt(L), log(1 + Z1(L)) with Z1 = sum_j f_j exp(-L eps_j) (single-insertion part), log Z_Fock(L),
 - count violations: sectors (Q,n) and energies E with N_true(E) > N_Fock(E), Fock built from the
   one-particle levels (bosonic multisets, additive energies, charge and momentum mod m conserved).
"""
import sys, os, json, time
import numpy as np
from scipy.special import logsumexp
from insertion import Sectors, single_insertion, one_particle_levels, fock_counts

HERE = os.path.dirname(os.path.abspath(__file__))
LS = [2, 4, 6, 8, 12, 16, 24, 32]


def analyse(name, beta, m, Ecut=None):
    from models import get
    h, mode, M, A, label = get(name)
    S = Sectors(name, beta, m, Qmax=4, h=h, mode=mode)
    Delta = name
    W = single_insertion(S)
    lev = one_particle_levels(S, W)
    e1 = np.array([l[0] for l in lev])
    if Ecut is None:
        Ecut = min(4.0, 3.0 * np.sort(e1)[min(3, len(e1) - 1)] + 0.5)
    allE = np.concatenate([S.eps(k) for k in S.data])
    res = dict(Delta=Delta, beta=beta, m=m, n_levels=len(lev), fmin=float(min(l[3] for l in lev)))
    w1 = []
    for key, w in W.items():
        top = np.argsort(-w['f'])[:w['dim']]
        w1.append(w['mu'][top].sum())
    res['weight1_min'] = float(min(w1))
    res['L'] = LS
    res['logZt'] = [float(logsumexp(-L * allE)) for L in LS]
    res['logZ1'] = [float(np.log1p(sum((w['f'] * np.exp(-L * w['eps'])).sum() for w in W.values()))) for L in LS]
    res['logZF'] = [float(-np.log1p(-np.exp(-L * e1)).sum()) for L in LS]
    # reduced Fock model: one level per (Q, n) sector (the largest-overlap one), i.e. d^2 - 1 = 3 per momentum
    e1r = []
    for key, w in W.items():
        e1r.append(w['eps'][int(np.argmax(w['f']))])
    e1r = np.array(e1r)
    res['logZF_one_per_sector'] = [float(-np.log1p(-np.exp(-L * e1r)).sum()) for L in LS]
    res['one_particle_levels'] = [[float(l[0]), int(l[1]), int(l[2]), float(l[3])] for l in lev]
    # singlet-type composites: lowest level in each sector with negligible single-insertion overlap
    comp = []
    for key, w in W.items():
        low = np.where(w['f'] < 1e-8)[0]
        if low.size:
            comp.append([int(key[0]), int(key[1]), float(w['eps'][low[0]])])
    res['lowest_zero_overlap_levels'] = comp
    # counting test
    while True:
        try:
            FE, FQ, Fn = fock_counts(lev, m, Ecut, qmod=2 if mode == 'z2' else None)
            break
        except RuntimeError:
            Ecut *= 0.85
    res['Ecut'] = float(Ecut)
    viol = []
    grid = np.linspace(0, Ecut, 81)[1:]
    for key in S.data:
        et = np.sort(S.eps(key))
        ef = np.sort(FE[(FQ == key[0]) & (Fn == key[1])])
        nt = np.searchsorted(et, grid, side='right'); nf = np.searchsorted(ef, grid, side='right')
        bad = np.where(nt > nf)[0]
        if bad.size:
            i = bad[0]
            viol.append(dict(Q=int(key[0]), n=int(key[1]), E=float(grid[i]), true=int(nt[i]), fock=int(nf[i]),
                             max_excess=int((nt - nf).max())))
    res['violations'] = viol
    # level-by-level domination with scaled one-particle energies: largest lambda with e_k >= lambda f_k
    # for every k (k-th level counted with multiplicity), globally and within each (Q, n) sector
    allE_s = np.sort(allE)
    FE_s = np.sort(FE)
    kmax = min(len(FE_s), len(allE_s))
    ratios = allE_s[1:kmax] / FE_s[1:kmax]
    kbad = int(np.argmin(ratios)) + 1
    res['lambda_global'] = float(ratios.min())
    res['lambda_global_at'] = [kbad, float(allE_s[kbad]), float(FE_s[kbad])]
    lam_sec = []
    for key in S.data:
        et = np.sort(S.eps(key))
        ef = np.sort(FE[(FQ == key[0]) & (Fn == key[1])])
        k = min(len(et), len(ef))
        st = 1 if key == (0, 0) else 0
        if k > st:
            r = et[st:k] / ef[st:k]
            j = int(np.argmin(r))
            lam_sec.append((float(r[j]), int(key[0]), int(key[1]), float(et[st + j]), float(ef[st + j])))
        elif len(et) > len(ef) and et[len(ef)] < Ecut:
            lam_sec.append((0.0, int(key[0]), int(key[1]), float(et[len(ef)]), float('inf')))
    lam_sec.sort()
    res['lambda_sector'] = lam_sec[0][0] if lam_sec else None
    res['lambda_sector_at'] = list(lam_sec[0][1:]) if lam_sec else None
    res['first_violation_E'] = float(min([v['E'] for v in viol])) if viol else None
    res['lowest_e1'] = float(np.sort(e1)[0])
    return res


if __name__ == '__main__':
    mlist = [int(x) for x in sys.argv[1].split(',')]
    blist = [float(x) for x in sys.argv[2].split(',')]
    for m in mlist:
        for b in blist:
            out = []
            for Delta in (sys.argv[3].split(',') if len(sys.argv) > 3 else ['XX', 'XXZ05', 'HEIS', 'XXZ2', 'TFIM', 'XYI', 'XYZ', 'XXZSTAG']):
                t0 = time.time()
                r = analyse(Delta, b, m)
                r['seconds'] = time.time() - t0
                out.append(r)
                print('m=%d beta=%g %s: fmin=%.3f w1min=%.3f  logZt=%s logZ1=%s logZF=%s  #viol=%d first_viol_E=%s lam_glob=%.3f lam_sec=%s (%.0fs)' % (
                    m, b, Delta, r['fmin'], r['weight1_min'], np.round(r['logZt'][:4], 3), np.round(r['logZ1'][:4], 3),
                    np.round(r['logZF'][:4], 3), len(r['violations']), r['first_violation_E'], r['lambda_global'], r['lambda_sector'], r['seconds']), flush=True)
            with open(os.path.join(HERE, 'data', 'insertion_m%d_b%g.json' % (m, b)), 'w') as f:
                json.dump(out, f, indent=1)
