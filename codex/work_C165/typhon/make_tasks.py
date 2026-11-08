#!/usr/bin/env python3
"""C165 (Claude Opus agent, 2026-11-08): task lists for the packed Typhon job.
tasks_long.txt : DMRG (L = 64..256) and VUMPS to chi = 256 at calibrations and the Q2 line (start first)
tasks_A.txt    : two-site VUMPS screen, chi 8..64, families A1 (charge-2 coexistence family), A2 (plateau edge),
                 A3 (random symmetry-free harmonics)
Each line: python3 <script> ... (run from the scratch directory)."""
import numpy as np
L, A = [], []
def fes(tag, kw, chis, n=2):
    return f"python3 runq.py \"{kw}\" fes/{tag}.json {' '.join(map(str, chis))} --n {n} --npz npz/{tag}.npz"
# ---- long tasks
for Lc, chi in ((128, 128), (256, 160)):
    for Dz in (1.15, 1.2, 1.25, 1.3):
        L.append(f"python3 dmrgq.py 'fam:a=0.9,Dz={Dz},g=0.3,hx=0.2,hz=0.05,Dx=0.03' {Lc} {chi} dm/Q2_Dz{Dz}_L{Lc}.json")
    L.append(f"python3 dmrgq.py 'fam:nfloat=dict(eta=0.1,b=1.05)' {Lc} {chi} dm/nfloat_L{Lc}.json")
    L.append(f"python3 dmrgq.py xx {Lc} {chi} dm/xx_L{Lc}.json")
L.append("python3 dmrgq.py tfi 128 64 dm/tfi_L128.json")
L.append("python3 dmrgq.py isxx 64 128 dm/isxx_L64.json")
L.append("python3 dmrgq.py isxx 128 128 dm/isxx_L128.json")
big = [8, 12, 16, 24, 32, 48, 64, 96, 128, 160, 192, 256]
L.append("python3 fes_c153.py isxx --chis " + ' '.join(map(str, big)) + " --outdir fes_cal")
L.append("python3 fesx.py xx --n 2 --chis " + ' '.join(map(str, big)) + " --out fes_cal/xx.json")
L.append("python3 fesx.py tfi --n 2 --chis " + ' '.join(map(str, big[:-2])) + " --out fes_cal/tfi.json")
L.append(fes('nfloat', 'nfloat=dict(eta=0.1,b=1.05)', big))
for Dz in (1.2, 1.25, 1.3):
    L.append(fes(f'Q2big_Dz{Dz}', f'a=0.9,Dz={Dz},g=0.3,hx=0.2,hz=0.05,Dx=0.03', big))
# ---- screen
sc = [8, 12, 16, 24, 32, 48, 64]
for a in (0.7, 0.9, 1.1):
    for g in (0.2, 0.35):
        for hz in (0.05, 0.3):
            for Dz in np.round(np.arange(1.0, 1.81, 0.1), 2):
                A.append(fes(f'A1_a{a}_g{g}_hz{hz}_Dz{Dz}', f'a={a},Dz={Dz},g={g},hx=0.2,hz={hz},Dx=0.03', sc))
for Dz, hzs in ((1.6, (0.3, 0.6, 0.9, 1.2, 1.4, 1.6, 1.8, 2.0)), (2.0, (0.6, 1.0, 1.4, 1.8, 2.1, 2.3, 2.5, 2.7))):
    for j, (g, hx, G) in enumerate(((0.3, 0.2, 0.0), (0.3, 0.4, 0.2), (0.6, 0.3, 0.1))):
        for hz in hzs:
            A.append(fes(f'A2_Dz{Dz}_h{j}_hz{hz}', f'a=0.9,Dz={Dz},g={g},hx={hx},G={G},hz={hz},Dx=0.03', sc))
rng = np.random.default_rng(165)
for i in range(16):
    v = rng.standard_normal(8); v *= 0.35 / np.linalg.norm(v)
    g, gxy, G, Gp, Dx, Dy, hx, hy = np.round(v, 4)
    for Dz in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5):
        A.append(fes(f'A3_v{i}_Dz{Dz}', f'a=0.9,Dz={Dz},g={g},gxy={gxy},G={G},Gp={Gp},Dx={Dx},Dy={Dy},hx={hx},hy={hy},hz=0.05', sc))
open('tasks_long.txt', 'w').write('\n'.join(L) + '\n')
open('tasks_A.txt', 'w').write('\n'.join(A) + '\n')
print(len(L), len(A))
