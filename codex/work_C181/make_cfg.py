"""C181: Typhon task list. Fixed-N open chains of model (3.2) with C176's fibdmrg.py (run.py driver)."""
import json, numpy as np, os
PI2 = float(np.pi / 2)
lo = [16, 16, 16, 16, 32, 32, 32, 32, 64, 64, 64, 64, 128, 128, 128, 192, 192, 256, 256]
def chis(L, cmax):
    s = list(lo)
    for c in (384, 512):
        if c <= cmax: s += [c, c]
    return s
tasks64, eos, t128N, t128P = [], [], [], []
def add(tag, L, N, wa, phi, cmax, meas=True, **kw):
    c = dict(tag=tag, L=L, N=N, wa=wa, phi=phi, chis=chis(L, cmax), **kw)
    if meas: c['meas_chi'] = 10 ** 9
    else: c['nomeas'] = True
    json.dump([c], open(f"typhon/cfg/{tag}.json", "w"))
    line = f"timeout 2100m python3 run.py typhon/cfg/{tag}.json out > logs/{tag}.log 2>&1"
    return line
points = []  # (label, wa, phi, nfrac)
for p in (PI2, 1.560, 1.553, 1.545, 1.535, 1.52, 1.50, 1.48):
    points.append((f"w2n932_p{p:.4f}", 2.0, p, 9 / 32))
for p in (PI2, 1.545, 1.52):
    points.append((f"w2n38_p{p:.4f}", 2.0, p, 0.375))
for p in (1.565, 1.555, 1.54, 1.50):
    points.append((f"w1n316_p{p:.4f}", 1.0, p, 3 / 16))
for p in (1.565, 1.55):
    points.append((f"w1n8_p{p:.4f}", 1.0, p, 1 / 8))
for (lab, wa, p, nf) in points:
    for L in (64, 128):
        N = int(round(nf * L))
        cmax = 512 if (wa == 1.0 or L == 64) else 384
        tgt = tasks64 if L == 64 else t128N
        tgt.append(add(f"{lab}_L{L}_N{N}", L, N, wa, p, cmax))
        tgtp = tasks64 if L == 64 else t128P
        for dN in (-1, 1):
            tgtp.append(add(f"{lab}_L{L}_N{N+dN}", L, N + dN, wa, p, min(cmax, 384), meas=False))
# controls at n = 3/16 (c = 1: one valley, no V; c = 2: two decoupled XX legs) and the V2 = V3 = 0 two-valley control
for L in (64, 128):
    N = 3 * L // 16
    tgt = tasks64 if L == 64 else t128N
    tgt.append(add(f"ctl1_L{L}_N{N}", L, N, 1.0, 0.0, 512, V2=0.0, V3=0.0))
    tgt.append(add(f"ctl2_L{L}_N{N}", L, N, 0.0, 0.0, 512, V2=0.0, V3=0.0, blockade=False))
    tgt.append(add(f"ctlV0_w2_p1.5350_L{L}_N{9*L//32}", L, 9 * L // 32, 2.0, 1.535, 384, V2=0.0, V3=0.0))
# equation of state at |w| = 2, L = 64
for p in (PI2, 1.535):
    for N in range(4, 30, 2):
        eos.append(add(f"eos_w2_p{p:.4f}_L64_N{N}", 64, N, 2.0, p, 256, meas=False))
lines = tasks64 + eos + t128N + t128P
open("typhon/tasks.txt", "w").write("\n".join(lines) + "\n")
print(len(tasks64), len(eos), len(t128N), len(t128P), len(lines))
