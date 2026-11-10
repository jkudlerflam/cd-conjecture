"""C176 estimators. reduce_meas() is called by run.py; the __main__ block summarizes jsonl outputs.
Usage: python3 analyze.py <dir with *.jsonl>
"""
import numpy as np, json, glob, os, sys

KG = np.linspace(-np.pi, np.pi, 256, endpoint=False)


def nk_from_G(G, ks=KG):
    W = G.shape[0]
    x = np.arange(W)
    ph = np.exp(1j * np.outer(ks, x))          # e^{ik i}
    return np.real(np.einsum('ki,ij,kj->k', ph, G, ph.conj())) / W


def reduce_meas(m, L, N, store_G=False):
    out = dict(norm=m['norm'], nn_adjacent_max=m['nn_adjacent_max'], dens=m['dens'], gwin=m['gwin'],
               P_exp_maxabs=float(np.max(np.abs(m['P_exp']))), P_exp=m['P_exp'], Cnn=m['Cnn'], Cpp=m['Cpp'], Dpair=m.get('Dpair'))
    for lab in ('f', 'b'):
        G = np.array(m[f'G{lab}_re']) + 1j * np.array(m[f'G{lab}_im'])
        nk = nk_from_G(G)
        out[f'nk_{lab}'] = nk.tolist()
        Wa = float(nk[np.cos(KG) > 0].sum()); Wb = float(nk[np.cos(KG) < 0].sum())
        out[f'Wa_{lab}'] = Wa; out[f'Wb_{lab}'] = Wb
        if store_G:
            out[f'G{lab}_re'] = m[f'G{lab}_re']; out[f'G{lab}_im'] = m[f'G{lab}_im']
    return out


# ---------------- summary estimators ----------------
def Sbar(S, L, W):
    S = np.asarray(S)
    l = np.arange(1, L)
    sel = (l >= L // 2 - W) & (l <= L // 2 + W)
    return float(S[sel].mean())


def c_prof(S, L, sm=12, lo=1 / 6, hi=5 / 6):
    S = np.asarray(S)
    l = np.arange(1, L)
    k = np.ones(sm) / sm
    Ss = np.convolve(S, k, mode='valid')
    ls = np.convolve(l.astype(float), k, mode='valid')
    x = np.log((2 * L / np.pi) * np.sin(np.pi * ls / L))
    sel = (ls >= lo * L) & (ls <= hi * L)
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    c6, a = np.linalg.lstsq(A, Ss[sel], rcond=None)[0]
    return 6 * c6


def c_prof_raw(S, L, lo=1 / 6, hi=5 / 6):
    S = np.asarray(S)
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    sel = (l >= lo * L) & (l <= hi * L)
    A = np.vstack([x[sel], np.ones(sel.sum())]).T
    return 6 * np.linalg.lstsq(A, S[sel], rcond=None)[0][0]


def cpp_at(meas, L, frac=0.25, half=4):
    vals = []
    for k, row in meas['Cpp'].items():
        i0 = int(k)
        r = np.arange(2, 2 + len(row))
        sel = (np.abs(r - frac * L) <= half)
        vals.append(np.mean(np.asarray(row)[sel]))
    return float(np.mean(vals))


def load(d):
    runs = {}
    for f in sorted(glob.glob(os.path.join(d, '*.jsonl'))):
        recs = [json.loads(l) for l in open(f) if l.strip()]
        if recs:
            runs[recs[0]['tag']] = recs
    return runs


def series_key(tag):
    """tags are <series>_L<L>[...]"""
    i = tag.rfind('_L')
    return tag[:i], int(tag[i + 2:].split('_')[0])


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else '.'
    runs = load(d)
    ser = {}
    for tag, recs in runs.items():
        s, _ = series_key(tag)
        ser.setdefault(s, {})[recs[0]['cfg']['L']] = recs
    for s in sorted(ser):
        print(f"=== {s}")
        last = {}
        for L in sorted(ser[s]):
            recs = ser[s][L]
            r = recs[-1]
            N = r['cfg']['N']
            Es = [(x['chi'], round(x['E'], 9), f"{x['trunc']:.1e}") for x in recs]
            meas = None
            for x in reversed(recs):
                if 'meas' in x:
                    meas = x['meas']; break
            line = f"L={L} N={N} sweeps={len(recs)}/{len(r['cfg']['chis'])} chi={r['chi']} trunc={r['trunc']:.1e} E={r['E']:.10f}"
            line += f" Sbar1={Sbar(r['S'], L, 1):.4f} Sbar6={Sbar(r['S'], L, 6):.4f} Sbar12={Sbar(r['S'], L, 12):.4f}"
            line += f" cprof={c_prof(r['S'], L):.3f} cprof_raw={c_prof_raw(r['S'], L):.3f}"
            if meas:
                Wa, Wb = meas['Wa_f'], meas['Wb_f']
                line += f" Wa/Wb={min(Wa, Wb) / max(Wa, Wb):.3f} CP(L/4)={cpp_at(meas, L):.2e} nnadj={meas['nn_adjacent_max']:.0e}"
                line += f" n_mid={np.mean(meas['dens'][L // 4: 3 * L // 4]):.4f}"
            print(line)
            print("    E(chi):", Es[-6:])
            last[L] = r
        Ls = sorted(last)
        for W in (1, 6, 12):
            cs = []
            for a in range(len(Ls)):
                for b in range(a + 1, len(Ls)):
                    L1, L2 = Ls[a], Ls[b]
                    c = 6 * (Sbar(last[L2]['S'], L2, W) - Sbar(last[L1]['S'], L1, W)) / np.log(L2 / L1)
                    cs.append(f"{L1}->{L2}:{c:.3f}")
            print(f"  c_two W={W}:", ' '.join(cs))
        # S at fixed x = l/L
        for x in (0.25, 0.5):
            print(f"  S(x={x}):", ' '.join(f"{L}:{last[L]['S'][int(x * L) - 1]:.4f}" for L in Ls))


# ---------------- registered criterion and extra diagnostics ----------------
def ppol(meas, L):
    n = np.mean(meas['dens'][L // 4: 3 * L // 4])
    return 2 * n  # cos(s) ~ 0.97 is dropped, which makes the threshold slightly stricter


def sq_peak(meas, L):
    """Density structure factor from the Cnn rows: S(q) = sum_r C(r) cos(q r) over 1 <= r <= L/4, averaged over refs."""
    qs = np.linspace(0, np.pi, 257)
    tot = np.zeros_like(qs)
    for k, row in meas['Cnn'].items():
        row = np.asarray(row)[: L // 4]
        r = np.arange(1, len(row) + 1)
        tot += 2 * np.cos(np.outer(qs, r)) @ row
    tot /= len(meas['Cnn'])
    i = np.argmax(tot)
    return float(qs[i]), float(tot[i])


def last_meas(recs):
    for x in reversed(recs):
        if 'meas' in x:
            return x['meas']
    return None


def criterion(ser, main='sym316', c1='c1', c2='c2'):
    out = []
    S = ser.get(main, {})
    if 128 not in S or 64 not in S:
        return ["criterion: main series lacks L = 64 or 128"]
    r64, r128 = S[64][-1], S[128][-1]
    m64, m128 = last_meas(S[64]), last_meas(S[128])
    ok1 = r128['trunc'] < 1e-7
    cts = {W: 6 * (Sbar(r128['S'], 128, W) - Sbar(r64['S'], 64, W)) / np.log(2) for W in (1, 6, 12)}
    ok2 = all(v > 1.2 for v in cts.values())
    cp = c_prof(r128['S'], 128); ok3 = cp > 1.2
    rat = min(m128['Wa_f'], m128['Wb_f']) / max(m128['Wa_f'], m128['Wb_f']); ok4 = rat >= 0.8
    c128, c64 = cpp_at(m128, 128), cpp_at(m64, 64)
    thr = 0.05 * ppol(m128, 128) ** 2
    ok5 = (abs(c128) < thr) and (abs(c128) < abs(c64))
    out.append(f"P1 trunc(128)={r128['trunc']:.1e} <1e-7: {ok1}")
    out.append(f"P2 c_two(64->128) W=1,6,12: {', '.join(f'{v:.3f}' for v in cts.values())} >1.2: {ok2}")
    out.append(f"P3 c_prof(128)={cp:.3f} >1.2: {ok3}")
    out.append(f"P4 Wmin/Wmax(128)={rat:.3f} >=0.8: {ok4}")
    out.append(f"P5 |C_P(L/4)| 64: {abs(c64):.2e}, 128: {abs(c128):.2e}, threshold {thr:.2e}: {ok5}")
    for nm, tgt, tol in ((c1, 1, 0.15), (c2, 2, 0.25)):
        C = ser.get(nm, {})
        if 64 in C and 128 in C:
            v = 6 * (Sbar(C[128][-1]['S'], 128, 6) - Sbar(C[64][-1]['S'], 64, 6)) / np.log(2)
            out.append(f"calibration {nm}: c_two(64->128, W=6) = {v:.3f} (target {tgt} +- {tol}): {abs(v - tgt) <= tol}")
    out.append(f"FIRST STAGE {'PASS' if all((ok1, ok2, ok3, ok4, ok5)) else 'FAIL'}")
    if 256 in S:
        r256 = S[256][-1]
        v = 6 * (Sbar(r256['S'], 256, 6) - Sbar(r128['S'], 128, 6)) / np.log(2)
        out.append(f"prediction check: c_two(128->256, W=6) = {v:.3f}, c_two(64->128, W=6) = {cts[6]:.3f}, trunc(256) = {r256['trunc']:.1e}")
        m256 = last_meas(S[256])
        if m256:
            out.append(f"   |C_P(L/4)| at 256: {abs(cpp_at(m256, 256)):.2e}")
    return out


def diagnostics(ser):
    out = []
    for s in sorted(ser):
        for L in sorted(ser[s]):
            recs = ser[s][L]
            m = last_meas(recs)
            # energy against truncation at the last sweep of each chi
            byc = {}
            for x in recs:
                byc[x['chi']] = (x['E'], x['trunc'])
            cs = sorted(byc)[-3:]
            if len(cs) >= 2:
                tr = np.array([byc[c][1] for c in cs]); E = np.array([byc[c][0] for c in cs])
                A = np.vstack([tr, np.ones_like(tr)]).T
                slope, E0 = np.linalg.lstsq(A, E, rcond=None)[0]
                dE = E[-1] - E0
            else:
                dE = np.nan
            line = f"{s} L={L}: E(last chis) {[(c, round(byc[c][0], 8)) for c in cs]}  E - E_extrap(trunc->0) = {dE:.1e}"
            if m:
                q, h = sq_peak(m, L)
                cp = [np.mean(np.asarray(row)[int(f * L) - 3:int(f * L) + 1]) for f in (1 / 16, 1 / 8, 1 / 4) for row in [m['Cpp'][sorted(m['Cpp'])[0]]]]
                line += f"  S_nn peak q={q:.3f} height={h:.3f}  C_P(L/16,L/8,L/4)={', '.join(f'{v:.1e}' for v in cp)}"
            out.append(line)
    return out


if __name__ == "__main__":
    print("\n=== registered first-stage criterion (main = sym316, n = 3/16) ===")
    for l in criterion(ser):
        print(l)
    print("\n=== diagnostics ===")
    for l in diagnostics(ser):
        print(l)


def e_at_chi(recs, chi):
    v = [x['E'] for x in recs if x['chi'] == chi]
    return v[-1] if v else None


def parity_table(ser):
    out = []
    for s in ('sym316', 'sym532', 'sym18'):
        if s not in ser:
            continue
        for L in sorted(ser[s]):
            P, M = ser.get(s + 'p1', {}).get(L), ser.get(s + 'm1', {}).get(L)
            if not P or not M:
                continue
            R = ser[s][L]
            common = sorted(set(x['chi'] for x in R) & set(x['chi'] for x in P) & set(x['chi'] for x in M))
            rows = []
            for chi in common[-2:]:
                E0, Ep, Em = e_at_chi(R, chi), e_at_chi(P, chi), e_at_chi(M, chi)
                d = Ep + Em - 2 * E0
                rows.append(f"chi={chi}: D_P={d:.6f} L*D_P={L * d:.4f}")
            out.append(f"{s} L={L} N={R[0]['cfg']['N']}: " + '; '.join(rows) +
                       f"  (complete sweeps {len(R)}/{len(R[0]['cfg']['chis'])}, {len(P)}/{len(P[0]['cfg']['chis'])}, {len(M)}/{len(M[0]['cfg']['chis'])})")
    return out


def envelope(y, w=8):
    y = np.abs(np.asarray(y))
    return np.array([y[max(0, i - w // 2): i + w // 2 + 1].max() for i in range(len(y))])


def pair_vs_single(ser, Lsel=None):
    out = []
    for s in sorted(ser):
        if s.endswith('p1') or s.endswith('m1') or s.startswith('c'):
            continue
        for L in sorted(ser[s]):
            if Lsel and L not in Lsel:
                continue
            recs = ser[s][L]
            m = None
            for x in reversed(recs):
                if 'meas' in x and 'Gf_re' in x['meas']:
                    m = x['meas']; break
            if m is None or not m.get('Dpair'):
                continue
            G = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
            a0, a1 = m['gwin']; W = a1 - a0
            i = W // 2 - 4
            Gr = np.array([np.mean([abs(G[i + d, i + d + r]) for d in range(4)]) for r in range(1, W - i - 4)])
            Ds = []
            for key, row in m['Dpair'].items():
                if key.startswith('2_'):
                    Ds.append(np.abs(np.array(row)[:, 0] + 1j * np.array(row)[:, 1]))
            nmin = min(len(d) for d in Ds)
            Dr = np.mean([d[:nmin] for d in Ds], axis=0)       # r = 3, 4, ...
            rD = np.arange(3, 3 + nmin); rG = np.arange(1, 1 + len(Gr))
            eG, eD = envelope(Gr), envelope(Dr)
            def slope(r, e, lo=8, hi=64):
                sel = (r >= lo) & (r <= min(hi, r.max()))
                return np.polyfit(np.log(r[sel]), np.log(e[sel]), 1)[0]
            def at(r, e, x):
                return e[np.argmin(np.abs(r - x))]
            hi = min(64, rG.max(), rD.max())
            if hi < 24:
                continue
            out.append(f"{s} L={L}: |G_f| env slope={slope(rG, eG, hi=hi):.2f} ratio(16/{hi})={at(rG, eG, 16) / at(rG, eG, hi):.1f};"
                       f" |D| env slope={slope(rD, eD, hi=hi):.2f} ratio(16/{hi})={at(rD, eD, 16) / at(rD, eD, hi):.1f};"
                       f" D/G at r=16, {hi}: {at(rD, eD, 16) / at(rG, eG, 16):.2e}, {at(rD, eD, hi) / at(rG, eG, hi):.2e}")
    return out


if __name__ == "__main__":
    print("\n=== parity gap ===")
    for l in parity_table(ser):
        print(l)
    print("\n=== pair vs single-particle correlations ===")
    for l in pair_vs_single(ser):
        print(l)
    print("\n=== secondary reading, n = 5/32 ===")
    for l in criterion(ser, main='sym532'):
        print(l)


def density_corrected(ser, ntargets=(0.15, 0.17, 0.19), series=('sym18', 'sym532', 'sym16', 'sym316'), W=6):
    """Secondary estimator (not registered): at each L, interpolate S_bar_W linearly in ln(n_mid) across the
    phi = pi/2 series, evaluate at a fixed bulk density, then form two-size c. Removes the O(1/L) drift of the
    bulk density at fixed N/L (open edges are depleted)."""
    out = []
    byL = {}
    for s in series:
        for L, recs in ser.get(s, {}).items():
            m = last_meas(recs)
            if m is None:
                continue
            nm = float(np.mean(m['dens'][L // 4: 3 * L // 4]))
            byL.setdefault(L, []).append((nm, Sbar(recs[-1]['S'], L, W), s))
    for nt in ntargets:
        vals = {}
        for L, pts in byL.items():
            pts = sorted(pts)
            xs = np.log([p[0] for p in pts]); ys = [p[1] for p in pts]
            if len(pts) >= 2 and xs[0] <= np.log(nt) <= xs[-1]:
                vals[L] = float(np.interp(np.log(nt), xs, ys))
        Ls = sorted(vals)
        cs = [f"{Ls[a]}->{Ls[b]}:{6 * (vals[Ls[b]] - vals[Ls[a]]) / np.log(Ls[b] / Ls[a]):.3f}" for a in range(len(Ls)) for b in range(a + 1, len(Ls))]
        out.append(f"n_bulk={nt}: S_bar_W{W} " + ' '.join(f"{L}:{vals[L]:.4f}" for L in Ls) + "  c_two " + ' '.join(cs))
    for L in sorted(byL):
        out.append(f"   L={L}: " + ', '.join(f"{s} n_mid={nm:.4f} S={sb:.4f}" for nm, sb, s in sorted(byL[L])))
    return out


if __name__ == "__main__":
    print("\n=== density-corrected two-size estimate (secondary, not registered) ===")
    for W in (6, 12):
        for l in density_corrected(ser, W=W):
            print(l)
