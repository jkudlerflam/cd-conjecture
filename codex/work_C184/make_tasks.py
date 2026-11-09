"""Task list for the C184 Typhon job (sorted longest first). Slice: parity class, h_z = 0 (mu = V), t = 1."""
import numpy as np
from model import spin_from_fermion
def pt(L, chis, phi, Dl, V):
    p = spin_from_fermion(1.0, phi, Dl, V, V)
    lab = f"ph{phi:.2f}_D{Dl:.2f}_V{V:.2f}"
    return f"pt {L} {chis} {lab} spin:" + ",".join(f"{x:.10f}" for x in p)
C192 = "16,32,64,96,128,160,192,224,224"; C96 = "16,32,64,96,128,160,160"; C48 = "16,32,64,96,128,128"
T = []
Vf = [round(v, 2) for v in np.arange(-2.8, -0.95, 0.1)]
Vc = [-2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2]
for Dl in [0.4, 0.45, 0.5, 0.55, 0.6, 0.7]:
    for V in Vc: T.append((192, pt(192, C192, 1.0, Dl, V)))
for Dl in [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]:
    for V in Vf:
        T.append((96, pt(96, C96, 1.0, Dl, V))); T.append((48, pt(48, C48, 1.0, Dl, V)))
for Dl in [0.4, 0.6, 0.8]:
    for V in [-2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0]:
        T.append((96, pt(96, C96, 1.2, Dl, V))); T.append((48, pt(48, C48, 1.2, Dl, V)))
for L, ch in [(192, C192), (96, C96), (48, C48)]:
    for lab in ['xx', 'tfi']: T.append((L, f"ctrl {L} {ch} {lab}"))
T.sort(key=lambda x: -x[0])
open('typhon/tasks.txt', 'w').write("\n".join(t for _, t in T) + "\n")
print(len(T), sum(1 for L, _ in T if L == 192))
