"""Analysis of spec2.py output: twist-relaxed, valley-labelled spectra.

For every level (theta_j, m, E, C1): M_est = C1 / cos(s) with valleys at k_a = -s, k_b = pi + s, rounded to the
nearest integer of the parity of N. Twist relaxation: for each (m, M) the minimum over the theta grid, refined by
a parabola through the grid minimum and its two neighbours (periodic in theta). The relaxed ground state fixes m_0,
and kappa = (m - m_0) mod N is the pseudospin chain momentum (q = 2 pi kappa / N) of an M = 0 level.
Prints, per run: relaxed gaps by (M, kappa), Delta_ps, the Neel/triplet-0 level at kappa = N/2, R, and the charge
twist stiffness v_J = pi L E''(theta) of the ground-state branch.
Usage: python3 analyze2.py file.jsonl [more files]"""
import sys, json, math, numpy as np
from collections import defaultdict


def valley_cos(phi, wa):
    ks = np.linspace(-np.pi, np.pi, 200001)
    e = 2 * wa * np.cos(ks - phi) - 2 * np.cos(2 * ks)
    # local minima
    i = np.nonzero((e < np.roll(e, 1)) & (e < np.roll(e, -1)))[0]
    i = i[np.argsort(e[i])][:2]
    return ks[i], e[i], np.abs(np.cos(ks[i])).mean()


def relax(levels, nth, N, cs, model, L):
    """Twist relaxation on the unfolded flux line. Gauge invariance gives (m, theta + 2 pi) = (m + N, theta), so a
    branch is labelled by kappa = m mod N and its energy is minimised over all m = kappa mod N and all theta_j.
    Returns dict (m_mod_N, M) -> (E_relaxed, (m, j) of the grid minimum, profile dict {(m, j): E})."""
    prof = defaultdict(dict)
    for (j, m, E, c) in levels:
        if model in ('rods', 'legs'):
            M = 0
        else:
            x = c / cs
            M = int(np.round((x - N % 2) / 2) * 2 + N % 2)
        key = (m, M)
        if j not in prof[key] or E < prof[key][j]:
            prof[key][j] = E
    G = math.gcd(N, L)
    rel = {}
    for (m, M), d in prof.items():
        for j, E in d.items():
            # neighbours on the unfolded line
            jm, mm = (j - 1, m) if j > 0 else (nth - 1, (m - N) % L)
            jp, mp = (j + 1, m) if j < nth - 1 else (0, (m + N) % L)
            a = prof.get((mm, M), {}).get(jm)
            c = prof.get((mp, M), {}).get(jp)
            Eref = E
            if a is not None and c is not None and a >= E and c >= E:
                den = a - 2 * E + c
                if den > 1e-12:
                    x = 0.5 * (a - c) / den
                    Eref = E - 0.25 * (a - c) * x
            key = (m % G, M)
            if key not in rel or Eref < rel[key][0]:
                rel[key] = (Eref, (m, j), d)
    return rel


def summarize(r, verbose=True):
    L, N, phi, wa, model, nth = r['L'], r['N'], r['phi'], r['wabs'], r['model'], r['nth']
    ks, es, cs = valley_cos(phi, wa)
    rel = relax(r['levels'], nth, N, cs, model, L)
    Eg = min(v[0] for v in rel.values())
    (m0, M0) = min(rel, key=lambda k: rel[k][0])
    out = dict(L=L, N=N, n=N / L, phi=phi, wa=wa, model=model, E0=Eg, m0=m0, M0=M0)
    # ground-state branch stiffness from the grid minimum and its unfolded neighbours
    (mg, jg) = rel[(m0, M0)][1]
    P = {}
    for (j, m, E, c) in r['levels']:
        P[(m, j)] = min(E, P.get((m, j), 1e9))
    jm, mm = (jg - 1, mg) if jg > 0 else (nth - 1, (mg - N) % L)
    jp, mp = (jg + 1, mg) if jg < nth - 1 else (0, (mg + N) % L)
    if (mm, jm) in P and (mp, jp) in P:
        h = 2 * np.pi / nth
        out['vJ'] = np.pi * L * (P[(mm, jm)] - 2 * P[(mg, jg)] + P[(mp, jp)]) / h ** 2
    # classes
    byM = defaultdict(list)
    for (m, M), v in rel.items():
        byM[M].append((v[0] - Eg, m))
    for M in byM: byM[M].sort()
    out['byM'] = {M: byM[M][:6] for M in sorted(byM)}
    if model not in ('rods', 'legs'):
        if N % 2 == 0:
            z = byM.get(0, [])
            kap = defaultdict(list)
            for (dE, m) in z:
                kap[(m - m0) % N].append(dE)  # m here is already m mod gcd(N,L)
            out['M0_by_kappa'] = {k: sorted(v)[:3] for k, v in sorted(kap.items())}
            p2 = [x for M in (2, -2) for x in byM.get(M, [])]
            out['D_ps'] = min(p2)[0] if p2 else None
            p4 = [x for M in (4, -4) for x in byM.get(M, [])]
            out['D_4'] = min(p4)[0] if p4 else None
            neel = kap.get(N // 2, [])
            out['D_neel'] = min(neel) if neel else None
            if out['D_ps'] and out['D_neel'] is not None:
                out['R'] = out['D_neel'] / out['D_ps']
            k1 = sorted(kap.get(1, []) + kap.get(N - 1, []))
            out['D_k1'] = k1[0] if k1 else None
        else:
            p1 = [x for M in (1, -1) for x in byM.get(M, [])]
            p3 = [x for M in (3, -3) for x in byM.get(M, [])]
            out['D_ps'] = (min(p3)[0] - min(p1)[0]) if (p1 and p3) else None
    else:
        z = byM[0]
        kap = defaultdict(list)
        for (dE, m) in z:
            kap[(m - m0) % N].append(dE)
        out['rods_k1'] = sorted(kap.get(1, []) + kap.get(N - 1, []))[:2]
    if verbose:
        print(f"\n== {model} L={L} N={N} n={N/L:.4f} phi={phi:.4f} |w|={wa}: E0_rel={Eg:.6f} (E0/L={Eg/L:.6f}) m0={m0} M0={M0} vJ={out.get('vJ', float('nan')):.3f}")
        for M in sorted(byM, key=abs):
            print(f"   M={M:+d}: " + "  ".join(f"{dE:.4f}@{m}" for dE, m in byM[M][:6]))
        if 'M0_by_kappa' in out:
            print("   M=0 by kappa: " + "  ".join(f"k{k}:" + ",".join(f"{x:.4f}" for x in v) for k, v in out['M0_by_kappa'].items()))
            print(f"   D_ps={out['D_ps']}  D_neel={out['D_neel']}  R={out.get('R')}  D_4={out['D_4']}  D_kappa1={out['D_k1']}")
            if out['D_ps']:
                print(f"   L*D_ps={L*out['D_ps']:.3f}  N*D_ps={N*out['D_ps']:.3f}")
        elif model not in ('rods', 'legs'):
            print(f"   (odd N) D_ps(|M|=3 above |M|=1) = {out['D_ps']}")
        else:
            print(f"   rods: lowest at dm=+-1: {out['rods_k1']}  -> v = L dE/2pi = {[L*x/(2*np.pi) for x in out['rods_k1']]}")
    return out


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        for line in open(fn):
            summarize(json.loads(line))
