"""Open-chain DMRG on the uniform-gauge Kitaev chain in a [110] field:
H = K sum S^x_j S^y_{j+1} - (h/sqrt2) sum (S^x_j + S^y_j), S = sigma/2, K = 1.
Field split half/half into the bond so the two end sites carry half the field.
Usage: python3 kit110_dmrg.py h L chis(comma) out.jsonl"""
import sys, json, time
import numpy as np
from dmrg import dmrg, expect_bonds
sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
I2 = np.eye(2)
h = float(sys.argv[1]); L = int(sys.argv[2]); chis = tuple(int(x) for x in sys.argv[3].split(',')); out = sys.argv[4]
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
f = -(h / np.sqrt(2)) * (sx + sy)
hb = np.kron(sx, sy) + 0.5 * (np.kron(f, I2) + np.kron(I2, f))
t0 = time.time()
r = dmrg(hb, L, chis=chis, seed=seed, verbose=True)
S = r['S']
ops = expect_bonds(r['M'], [(sx, I2), (sy, I2), (sz, I2), (sx, sy)])
l = np.arange(1, L); x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
def fit(w):
    m = (l >= w * L) & (l <= (1 - w) * L)
    return [round(float(6 * np.polyfit(x[m & (l % 2 == p)], S[m & (l % 2 == p)], 1)[0]), 3) for p in (0, 1)]
rec = dict(h=h, L=L, chis=chis, seed=seed, E=r['E'], hist=r['hist'], S=S.tolist(), S_mid=float(S[L // 2 - 1]),
           c_fit={str(w): fit(w) for w in (0.05, 0.125, 0.25)},
           mx=ops[0].real.tolist(), my=ops[1].real.tolist(), mz=ops[2].real.tolist(), xy=ops[3].real.tolist(), time=time.time() - t0)
print(f"h={h} L={L} E={r['E']:.10f} S_mid={rec['S_mid']:.4f} c_fit={rec['c_fit']} time={rec['time']:.0f}s")
open(out, 'a').write(json.dumps(rec) + "\n")
