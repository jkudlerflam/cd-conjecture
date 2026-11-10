# KT drift bound (PROOF §5) from candidates.json, plus the unregistered C_P amplitude ratio.
import json, numpy as np, c193 as C
cand = json.load(open('candidates.json')); ln2 = np.log(2)
def bound(K1, K2, dK1, dK2, dl):
    dK = K2 - K1; sig = np.hypot(dK1, dK2); Kb = 0.5 * (K1 + K2)
    y2c = max(0.0, -2 * dK / (Kb ** 2 * dl)); y2u = max(0.0, -2 * (dK - sig) / (Kb ** 2 * dl))
    return dK, sig, np.sqrt(y2c), np.sqrt(y2u), 2 * (K2 - 1), 2 * (K2 - dK2 - 1)
print('key  branch  Ks128  Ks256  dK(sig)  |y1| central/upper   y_par central/lower   rdrift: deta(R3->R4) dK/dl')
for key in ['p153c', 'sym40_n316', 'rat_n316', 'sym40_n18', 'p153_n316', 'sym80_n316']:
    rows = []
    for lab, pick in (('ground', None), ('P0', '0'), ('P1', '1')):
        o = [cand[f'{key}_L{L}']['sec'][pick or str(cand[f'{key}_L{L}']['ground'])] for L in (128, 256)]
        dK, sig, yc, yu, yp, ypl = bound(o[0]['Ks'], o[1]['Ks'], o[0]['dKs'], o[1]['dKs'], ln2)
        # r-scale drift at L = 256: R3 [8,28] -> R4 [12,34], dl = ln(sqrt(12*34)/sqrt(8*28))
        de = np.mean([p['eta'][3] - p['eta'][2] for p in o[1]['b']]); dl = np.log(np.sqrt(12 * 34 / (8 * 28)))
        K = o[1]['Ks']; dKdl = -4 * K * K * de / dl; yr = np.sqrt(max(0, -2 * dKdl / K ** 2))
        print(f"{key:11s} {lab:6s} {o[0]['Ks']:.3f} {o[1]['Ks']:.3f} {dK:+.3f}({sig:.3f})  {yc:.3f}/{yu:.3f}   {yp:.3f}/{ypl:.3f}   "
              f"{de:+.4f} {dKdl:+.3f} |y1|_r={yr:.3f} remaining dK={dKdl/(2*(2*K-2)):+.3f}")
    m = []
    for L in (128, 256):
        g = cand[f'{key}_L{L}']; mm = C.load(C.GPU + g['sec'][str(g['ground'])]['file'])['meas']
        cp = np.array(mm['CP']); r = cp[:, 0]; v = cp[:, 1]; i0 = mm['CP_i0']; F = C.Fform(i0 + 1., i0 + 1. + r, L); h = C.sm(v)
        s = (r >= 16) & (r <= int(.45 * L)); m.append((np.exp(np.mean(np.log(np.abs(h[s])) + 2 * np.log(F[s]))), np.mean(np.sign(h[s]))))
    print(f"   C_P smooth amplitude at fixed exponent 2: A128={m[0][0]:.4g} A256={m[1][0]:.4g} ratio={m[1][0]/m[0][0]:.3f} sign={m[0][1]:+.2f},{m[1][1]:+.2f}")
# effective model cross-check, C192 numbers (full window): (eta_b128, eta_b256, pair128, pair256)
for lab, e1, e2, p1, p2 in (('eff 3/16 fixed 0.62', .325, .326, .62, .62), ('eff 3/16 measured pair', .325, .326, .635, .619), ('eff 1/8 fixed 0.62', .347, .351, .62, .62)):
    K1, K2 = 1 / (4 * e1 - p1), 1 / (4 * e2 - p2); d = 4 * 0.002 * K2 ** 2
    dK, sig, yc, yu, yp, ypl = bound(K1, K2, d, d, ln2)
    print(f"{lab:24s} Ks {K1:.3f} {K2:.3f} dK={dK:+.3f}  |y1| {yc:.3f}/{yu:.3f}  y_par {yp:.3f}")
