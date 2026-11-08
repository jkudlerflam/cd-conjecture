"""C154: DMRG along the spiral-frame path J(lam), B(lam) from the U(1) parent (lam = 0) to a bond (lam = 1)
(Claude Opus agent, 2026-11-07). Usage: python3 path_dmrg.py point L chis out.jsonl lam1,lam2,...
point = '1282' (C68 symmetry-free floating point) or a JSON file with keys J, B.
Records E, S(l), the window-averaged entropy Sbar = mean S(l) over l in [L/4, 3L/4] (C68/d2 convention:
for a CFT Sbar(2L) - Sbar(L) = (c/6) ln 2 up to corrections), and the energy gap is not computed."""
import sys, json, time
import numpy as np
from xzlib import bond
from dmrg import dmrg
from spiral import path, parent

pt, L, chis, fout = sys.argv[1], int(sys.argv[2]), tuple(int(c) for c in sys.argv[3].split(",")), sys.argv[4]
lams = [float(x) for x in sys.argv[5].split(",")]
if pt == '1282':
    d = json.load(open('../work_C68/checks.json'))['point_1282']
else:
    d = json.load(open(pt))
J1, B1 = np.array(d['J']), np.array(d['B'])
for lam in lams:
    J, B = path(J1, B1, lam)
    t0 = time.time()
    r = dmrg(bond(J, B), L, chis=chis)
    S = np.array(r['S']); l = np.arange(1, L)
    w = (l >= L // 4) & (l <= 3 * L // 4)
    rec = dict(point=pt, lam=lam, L=L, chis=list(chis), E_per_site=r['E'] / L, trunc=r['hist'][-1][2],
               dE_last=(r['hist'][-1][1] - r['hist'][-2][1]), Sbar=float(S[w].mean()), S=[round(float(s), 6) for s in S],
               sec=round(time.time() - t0, 1))
    with open(fout, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print(f"{pt} lam={lam} L={L}: E/L={rec['E_per_site']:.8f} trunc={rec['trunc']:.1e} dE={rec['dE_last']:.1e} Sbar={rec['Sbar']:.4f} S_mid={S[L//2-1]:.4f} ({rec['sec']}s)", flush=True)
