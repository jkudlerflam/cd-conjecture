"""C189 secondary (unregistered) estimator: pseudospin density correlations from the valley-polarization operator
P_j = c+_{j+1} c_j + h.c. (c frame), whose long-wavelength part is 2 cos(k_v) Z rho_s (PROOF 4.4).
Smooth part: <P_i P_j>_c = -Z^2 (2 cos k_v)^2 K_s / (pi^2 r^2). Z is calibrated on ff2v (exact Wick) at the same L
and N/L, so the output is K_s Z_cand^2 / Z_ff^2. Usage: python3 kest_cp.py > kest_cp.out"""
import numpy as np, json, warnings
from kest import *
warnings.simplefilter('ignore')
D = '../work_C176/typhon_out/out/'


def cp_rows_ff(L, N):
    h = ff_hamiltonian(L); e, U = np.linalg.eigh(h); occ = U[:, :N]
    G = np.conj(occ) @ occ.T                     # <c+_i c_j>
    Gh = np.eye(L) - G.T                         # <c_i c+_j>
    refs = [L // 2 - 3, L // 2 - 2, L // 2 - 1, L // 2]
    rows = {}
    for i0 in refs:
        row = []
        A = [(i0 + 1, i0), (i0, i0 + 1)]         # c+_a c_b terms of P_i0
        for j in range(i0 + 2, L - 1):
            B = [(j + 1, j), (j, j + 1)]
            v = 0
            for (a, b) in A:
                for (c, d) in B:
                    v += G[a, d] * Gh[b, c]       # connected Wick contraction <c+_a c_d><c_b c+_c>
            row.append(float(np.real(v)))
        rows[i0] = np.array(row)
    return rows


def cp_smooth(rows, L, frac, half=4):
    vals = []
    for i0, row in rows.items():
        r = np.arange(2, 2 + len(row))
        sel = np.abs(r - frac * L) <= half
        vals.append(np.mean(row[sel]))
    return float(np.mean(vals))


def cp_rows_run(path):
    recs = [json.loads(l) for l in open(path) if l.strip()]
    for x in reversed(recs):
        if 'meas' in x:
            return {int(k): np.array(v) for k, v in x['meas']['Cpp'].items()}


print('C_P smooth amplitude A = -C_P(r) r^2 averaged over |r - xL| <= 4; K_s Z^2/Z_ff^2 = A_cand / A_ff (ff2v: K_s = 1)')
for s, (p, q) in (('sym316', (3, 16)), ('sym532', (5, 32)), ('sym18', (1, 8))):
    for L in (128, 192, 256):
        N = p * L // q
        rc = cp_rows_run(D + f'{s}_L{L}.jsonl'); rf = cp_rows_ff(L, N)
        out = []
        for x in (1 / 8, 3 / 16, 1 / 4):
            ac = -cp_smooth(rc, L, x) * (x * L) ** 2; af = -cp_smooth(rf, L, x) * (x * L) ** 2
            out.append(f'x={x:.3f}: A_cand={ac:.4f} A_ff={af:.4f} ratio={ac / af:.3f}')
        print(f'{s} L={L}: ' + ' | '.join(out))


def Pmid_ff(L, N, phi):
    h = ff_hamiltonian(L, phi); e, U = np.linalg.eigh(h); occ = U[:, :N]
    G = np.conj(occ) @ occ.T
    Pj = np.array([2 * np.real(G[j + 1, j]) for j in range(L - 1)])   # <c+_{j+1} c_j + h.c.>
    r = ff_run(L, N, phi)
    p = (r['Wa'] - r['Wb']) / (r['Wa'] + r['Wb'])
    return float(np.mean(Pj[L // 4: 3 * L // 4])), p, n_mid(r)


print('\nCalibrated version: A = mean of -C_P(r) r^2 over L/8 <= r <= L/4 (many oscillation periods);'
      ' Z_P = <P>_mid / (p n_mid) from the phi = 1.54 run; K_s^P = pi^2 A / Z_P^2')
for L in (128, 256, 512):
    N = 3 * L // 16
    rf = cp_rows_ff(L, N)
    A = np.mean([np.mean(-row[(np.arange(2, 2 + len(row)) >= L // 8) & (np.arange(2, 2 + len(row)) <= L // 4)] *
                         np.arange(2, 2 + len(row))[(np.arange(2, 2 + len(row)) >= L // 8) & (np.arange(2, 2 + len(row)) <= L // 4)] ** 2) for row in rf.values()])
    P, p, nm = Pmid_ff(L, N, 1.54)
    Z = P / (p * nm)
    print(f'ff2v L={L}: A={A:.4f} <P>={P:+.4f} p={p:+.4f} Z_P={Z:.3f} K_s^P={np.pi ** 2 * A / Z ** 2:.3f} (exact 1)')


def A_window(rows, L, lo, hi):
    vals = []
    for row in rows.values():
        r = np.arange(2, 2 + len(row)); sel = (r >= lo) & (r <= hi)
        vals.append(np.mean(-row[sel] * r[sel] ** 2))
    return float(np.mean(vals))


def zp_from(path):
    recs = [json.loads(l) for l in open(path) if l.strip()]
    r = recs[-1]; m = [x for x in recs if 'meas' in x][-1]['meas']; L = r['cfg']['L']
    P = float(np.mean(m['P_exp'][L // 4: 3 * L // 4])); p = (m['Wa_f'] - m['Wb_f']) / (m['Wa_f'] + m['Wb_f'])
    nm = float(np.mean(m['dens'][L // 4: 3 * L // 4]))
    return P / (p * nm)


print('\nK_s^P = pi^2 A / Z_P^2 for the candidates, A over [L/8, L/4] (neighbours [L/16, L/4], [L/8, 3L/8])')
Z316 = [zp_from(D + 'p154_L128.jsonl'), zp_from(D + 'p154_L256.jsonl'), zp_from(D + 'p152_L256.jsonl'), zp_from(D + 'p150_L256.jsonl'),
        zp_from('../work_C187q/desktop_run3/out/w1n316_p1.5300_L192_N36.jsonl')]
Z18 = [zp_from('../work_C181/typhon_out/out/w1n8_p1.5500_L128_N16.jsonl')]
print('Z_P at 3/16 (phi = 1.54 L=128, 1.54 L=256, 1.52 L=256, 1.50 L=256, 1.53 L=192):', np.round(Z316, 3), ' at 1/8 (1.55, L=128):', np.round(Z18, 3))
Zs = {'sym316': np.mean(Z316), 'sym532': 0.5 * (np.mean(Z316) + Z18[0]), 'sym18': Z18[0]}
sZ = np.std(Z316)
for s, (p, q) in (('sym316', (3, 16)), ('sym532', (5, 32)), ('sym18', (1, 8))):
    for L in (128, 192, 256):
        rc = cp_rows_run(D + f'{s}_L{L}.jsonl')
        A = [A_window(rc, L, L // 8, L // 4), A_window(rc, L, L // 16, L // 4), A_window(rc, L, L // 8, 3 * L // 8)]
        Af = A_window(cp_rows_ff(L, p * L // q), L, L // 8, L // 4)
        Pf, pf, nf = Pmid_ff(L, p * L // q, 1.54); Zf = Pf / (pf * nf); Kf = np.pi ** 2 * Af / Zf ** 2
        K = [np.pi ** 2 * a / Zs[s] ** 2 for a in A]
        sig = np.sqrt((0.5 * (max(K) - min(K))) ** 2 + (K[0] * (1 - Kf)) ** 2 + (2 * K[0] * sZ / Zs[s]) ** 2)
        print(f'{s} L={L}: A={A[0]:.4f} K_s^P={K[0]:.3f} +- {sig:.3f}  [windows {K[1]:.3f}, {K[2]:.3f}; ff2v reads {Kf:.3f}; Z_P={Zs[s]:.3f}]')
