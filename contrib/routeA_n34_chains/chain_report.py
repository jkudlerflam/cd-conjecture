"""Weak-coupling report on a logical chain H = sum_j K_{j,j+1}:
Hartree-Fock Fermi points, velocities and centres; the closest four-fermion resonances with their
amplitudes; the dimensionless forward-scattering (Luttinger) matrix.  K is normalized so that
max|T| = 1, the field mu adds mu * (number of cell j + number of cell j+1)/2 per bond.
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28.
"""
import itertools, json
import numpy as np
from analyze34 import nbody, free_chain, fermi, band_data, hartree_fock, vertex, resonances, wrap, ops


def number_bond(n):
    cs = ops(2 * n)
    N = sum(c.conj().T @ c for c in cs)
    return 0.5 * N


def with_field(K, n, mu):
    return K + mu * number_bond(n)


def report(K, n, mu=0.0, hf=True, top=10, verbose=True, W=None):
    Kf = with_field(K, n, mu)
    Ops, Ws = nbody(Kf) if W is None else W
    T0, M0, _ = free_chain(Kf, n, Ws[1])
    if hf:
        T, M, hist = hartree_fock(Kf, n, iters=60, mix=0.5, tol=1e-11, W=Ws)
    else:
        T, M, hist = T0, M0, []
    fp = fermi(T, M)
    bd = band_data(T, M)
    out = dict(n=n, mu=mu, hf=hf, hf_last_change=hist[-1] if hist else None,
               body_norms=[float(np.linalg.norm(O, 2)) for O in Ops],
               n_fermi_points=len(fp))
    out['fermi_points'] = [dict(k=f['k'], band=f['band'], v=f['v']) for f in fp]
    # centres from the two Fermi points of each sorted band
    cen = []
    for b in range(n):
        ks = sorted([f['k'] for f in fp if f['band'] == b])
        if len(ks) == 2:
            # occupied arc: the band is negative between the points on one side
            kL, kR = ks
            mid = 0.5 * (kL + kR)
            from analyze34 import hk
            e_mid = np.linalg.eigvalsh(hk(T, M, mid))[b]
            c = mid if e_mid < 0 else wrap(mid + np.pi)
            width = (kR - kL) if e_mid < 0 else 2 * np.pi - (kR - kL)
            cen.append(dict(band=b, centre=float(c), filling=float(width / (2 * np.pi))))
        else:
            cen.append(dict(band=b, centre=None, n_points=len(ks)))
    out['bands'] = cen
    out['band_widths'] = [x['width'] for x in bd]
    if len(fp) >= 2:
        vs = np.array([abs(f['v']) for f in fp])
        out['v_max'], out['v_min'] = float(vs.max()), float(vs.min())
        closest, allr = resonances(fp, top=10 ** 6)
        Op2 = Ops[2]
        ph = [(f['k'], f['u']) for f in fp]
        # forward-scattering matrix (dimensionless)
        P = len(fp)
        A = np.eye(P)
        for p in range(P):
            for q in range(P):
                if p != q:
                    g = vertex(Op2, n, ph[p], ph[q], ph[p], ph[q]).real
                    A[p, q] = g / (2 * np.pi * np.sqrt(vs[p] * vs[q]))
        out['forward_matrix_min_eig'] = float(np.linalg.eigvalsh(0.5 * (A + A.T)).min())
        out['forward_max_offdiag'] = float(np.max(np.abs(A - np.eye(P))))
        # Luttinger parameters per sorted band, from the R-L forward coupling
        Ks = []
        for b in range(n):
            idx = [i for i, f in enumerate(fp) if f['band'] == b]
            if len(idx) == 2:
                Ks.append(float(np.sqrt((1 - A[idx[0], idx[1]]) / (1 + A[idx[0], idx[1]]))) if abs(A[idx[0], idx[1]]) < 1 else None)
        out['K_bands_estimate'] = Ks
        rows = []
        for r in allr:
            a, b, c, d = r['points']
            if a == b or c == d:
                continue                        # contact amplitude vanishes (Pauli); dimension >= 3
            g = vertex(Op2, n, ph[a], ph[b], ph[c], ph[d])
            gam = abs(g) / (2 * np.pi * (vs[a] * vs[b] * vs[c] * vs[d]) ** 0.25)
            rows.append(dict(points=[a, b, c, d], kind=r['kind'], mismatch=r['mismatch'],
                             gamma=float(gam), flow=float(gam * np.log(1 / max(r['mismatch'], 1e-12)))))
        rows.sort(key=lambda r: r['mismatch'])
        out['closest_resonances'] = rows[:top]
        out['max_flow_indicator'] = max((r['flow'] for r in rows), default=None)
        out['worst_channel'] = max(rows, key=lambda r: r['flow']) if rows else None
        reps = [r for r in allr if (r['points'][0] == r['points'][1] or r['points'][2] == r['points'][3])]
        out['closest_repeated_point_channel'] = reps[0] if reps else None
    out['T_hf'] = [[complex(x) for x in row] for row in T]
    out['M_hf'] = [[complex(x) for x in row] for row in M]
    if verbose:
        print('n=%d mu=%.3e HF=%s  #FP=%d  body norms %s' % (n, mu, hf, len(fp), ['%.1e' % x for x in out['body_norms']]))
        for f in fp:
            print('   k=%+.6f band %d v=%+.3e' % (f['k'], f['band'], f['v']))
        print('   bands', [(c['band'], None if c['centre'] is None else round(c['centre'], 5), round(c.get('filling', -1), 4)) for c in cen])
        if len(fp) >= 2:
            print('   forward: min eig %.4f max offdiag %.3e  K %s' % (out['forward_matrix_min_eig'], out['forward_max_offdiag'], out['K_bands_estimate']))
            for r in out['closest_resonances'][:6]:
                print('   res %s %-40s mismatch %.2e gamma %.2e flow %.2e' % (r['points'], r['kind'], r['mismatch'], r['gamma'], r['flow']))
            print('   worst', out['worst_channel'])
            print('   repeated-point closest', out['closest_repeated_point_channel'])
    return out


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o
