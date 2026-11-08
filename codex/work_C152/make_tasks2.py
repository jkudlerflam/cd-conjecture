"""Task list for Typhon job 2 (finite-size scaling at the L = 64 anomalies). Fermion units t = 1, Delta = s sin(phi), mu = V + 2h.
Set 1: h = 0 near XZ (phi = 1.45 and the Ashkin-Teller plane phi = pi/2), L = 32, 64, 128.
Set 2: attractive transition at phi = 1.3, V = -0.4, s = 1.0, 1.1, fine h scan, L = 32, 64, 128.
Set 3: strong attraction V = -1.5, Delta = 0.5, phi = 0, 0.8, mu scan, L = 48, 64.
Set 4: Peotta et al. convergence region, phi = 1, h = 0, Delta = 0.25, 0.32, 0.4, V = -2.3 to -1.6, L = 48, 96.
Set 5: XZ-window points with h_z = 0.05, 0.1 (phi = 1.45, s = 1.0, 1.1, V = 0.4, 0.8), attractive transition at phi = 1.45,
s = 1.1, V = -0.4, h_z = 0.03 to 0.07 (L = 32, 64, 128), and L = 192 at the two largest h_z = 0 fits."""
import numpy as np
chis = {32: "16,32,64,96,96", 48: "16,32,64,96,96", 64: "16,32,64,96,128,128", 128: "16,32,64,96,128,160,160"}
T = []
def add(L, lab, phi, s, V, mu, Dl=None):
    D = s * np.sin(phi) if Dl is None else Dl
    T.append((L, f"{L} {chis[L]} {lab} 1.0 {phi:.6f} {D:.6f} {V} {mu:.6f}"))
pts = [(phi, s, V) for phi in (1.45, np.pi / 2) for s in (0.9, 1.0) for V in (-0.8, -0.4, 0.4, 0.8)]
pts += [(1.45, 1.1, V) for V in (-0.8, -0.4, 0.4, 0.8)]
for (phi, s, V) in pts:
    for L in (128, 64, 32):
        add(L, f"fss_phi{phi:.3f}_s{s}_V{V}_h0", phi, s, V, V)
for s in (1.0, 1.1):
    for h in np.arange(0.11, 0.195, 0.01):
        for L in (128, 64, 32):
            add(L, f"atr_phi1.3_s{s}_V-0.4_h{h:.2f}", 1.3, s, -0.4, -0.4 + 2 * h)
for phi in (0.0, 0.8):
    for mu in np.linspace(-2.5, 0.5, 16):
        for L in (64, 48):
            add(L, f"sat_phi{phi}_D0.5_V-1.5_mu{mu:.2f}", phi, None, -1.5, mu, Dl=0.5)
# Set 4: the region where the phases of Peotta et al. (J. Stat. Mech. 2014, P09005) converge: phi = 1, h_z = 0,
# t = J, Delta = J_Delta, V = J_z in their notation; L = 48, 96.
chis[96] = "16,32,64,96,128,128"
for Dl in (0.25, 0.32, 0.4):
    for V in np.round(np.arange(-2.3, -1.55, 0.1), 2):
        for L in (96, 48):
            add(L, f"pmz_phi1.0_D{Dl}_V{V}_h0", 1.0, None, V, V, Dl=Dl)
# Set 5 (added after the full L = 64 table): h_z > 0 points of the XZ window with large L = 64 fits, the attractive
# transition at phi = 1.45, and L = 192 at the two largest h_z = 0 anomalies.
chis[192] = "16,32,64,96,128,160,192,192"
for s in (1.0, 1.1):
    for V in (0.4, 0.8):
        for h in (0.05, 0.1):
            for L in (128, 64, 32):
                add(L, f"fss_phi1.450_s{s}_V{V}_h{h}", 1.45, s, V, V + 2 * h)
for h in (0.03, 0.04, 0.05, 0.06, 0.07):
    for L in (128, 64, 32):
        add(L, f"atr_phi1.45_s1.1_V-0.4_h{h:.2f}", 1.45, 1.1, -0.4, -0.4 + 2 * h)
add(192, "fss_phi1.450_s1.1_V0.4_h0", 1.45, 1.1, 0.4, 0.4)
add(192, "fss_phi1.450_s0.9_V-0.8_h0", 1.45, 0.9, -0.8, -0.8)
T.sort(key=lambda x: -x[0])
open("typhon/dtasks2.txt", "w").write("\n".join(x[1] for x in T) + "\n")
print(len(T), sum(1 for x in T if x[0] == 128))
