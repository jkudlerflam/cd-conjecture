"""Reproduces the numerical claims of RESULT.md (quick version, a few minutes on a laptop).

Run: python3 run_checks.py      (add --long for L = 14 to 20 and SU(3) L = 9)
"""
import sys
import numpy as np
import scipy.linalg as sla
import fr_lib as F

LONG = '--long' in sys.argv
fails = []


def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    if not ok:
        fails.append(name)


def S(r):
    w = sla.eigvalsh((r + r.conj().T) / 2, driver='evr'); w = w[w > 1e-14]
    return float(-np.sum(w * np.log(w)))


# 1. Box smoothing, XX: E_L >= -csc(pi/L), equality on even L
h = F.xxz_bond(0.0)
ok = True
for L in range(2, 13):
    E0 = F.ground_states(h, L)[2]
    b = -1 / np.sin(np.pi / L)
    ok &= E0 >= b - 1e-10
    if L % 2 == 0:
        ok &= abs(E0 - b) < 1e-9
check('1 XX box smoothing E_L >= -csc(pi/L), equality on even L = 2..12', ok)

# 2. Petz defect L^2 (e_trial - E_L/L)/M at l = L/2 for d = 2 critical chains
Ls = [6, 8, 10, 12] + ([14, 16, 18, 20] if LONG else [])
e0 = {'XX': -1 / np.pi, 'Heis': 0.25 - np.log(2)}
Avals = {'XX': 1.0, 'XXZ0.5': np.pi * np.sin(np.pi / 3) / (np.pi / 3) / 2, 'Heis': np.pi / 2}
for name, hb in [('XX', F.xxz_bond(0.0)), ('XXZ0.5', F.xxz_bond(0.5)), ('Heis', F.heis_bond()), ('XXZ2', F.xxz_bond(2.0))]:
    M = F.M_of_bond(hb)
    vals = []
    ok = True
    for L in Ls:
        w, v, E0 = F.ground_states(hb, L, k=2)
        et = F.petz_defect(hb, v, L, L // 2)
        vals.append(L ** 2 * (et - E0 / L) / M)
        if name in e0:
            ok &= et >= e0[name] - 1e-12      # trial is a Z state: e_trial >= e_0
    dec = all(vals[i + 1] <= vals[i] + 1e-9 for i in range(len(vals) - 1))
    extra = '' if name not in Avals else ', pi A/(3M) = %.4f' % (np.pi * Avals[name] / (3 * M))
    print('   %s M=%.4f  L^2 delta/M = %s%s' % (name, M, ' '.join('%.4f' % x for x in vals), extra))
    check('2 %s Petz defect bounded and nonincreasing in L, e_trial >= e_0 where known' % name, ok and dec and max(vals) < 5)

# 3. SU(3) swap on vacuum-compatible lengths
hb = F.swap_bond(3)
vals = []
for L in ([3, 6, 9] if LONG else [3, 6]):
    w, v, E0 = F.ground_states(hb, L, k=3)
    vals.append(min(L ** 2 * (F.petz_defect(hb, v, L, l) - E0 / L) for l in range(2, min(L, 6))))
print('   SU(3) L^2 delta/M =', ' '.join('%.4f' % x for x in vals))
check('3 SU(3) Petz defect decreasing on L = 3, 6 (and 9)', all(vals[i + 1] < vals[i] for i in range(len(vals) - 1)))

# 4. Ferromagnetic Heisenberg, full ground multiplet mixture (A = 0)
hb = -F.heis_bond()
vals = []
for L in [4, 6, 8, 10]:
    w, v, E0 = F.ground_states(hb, L, k=L + 3)
    assert v.shape[1] == L + 1
    vals.append(min(L ** 2 * (F.petz_defect(hb, v, L, l) - E0 / L) / 0.5 for l in range(2, L)))
print('   FM Heisenberg L^2 delta/M =', ' '.join('%.4f' % x for x in vals))
check('4 FM Heisenberg multiplet: defect bounded (< 1.5) and decreasing', max(vals) < 1.5 and all(vals[i + 1] < vals[i] for i in range(3)))

# 5. Massive modes: TFIM off criticality is exponentially suppressed
Z = np.diag([1., -1]); X = np.array([[0, 1.], [1, 0]]); I = np.eye(2)
for g in [1.0, 1.6]:
    hb = -np.kron(Z, Z) - g * (np.kron(X, I) + np.kron(I, X)) / 2
    M = F.M_of_bond(hb)
    vals = []
    for L in [8, 10, 12]:
        w, v, E0 = F.ground_states(hb, L, k=2)
        vals.append(L ** 2 * (F.petz_defect(hb, v, L, L // 2) - E0 / L) / M)
    print('   TFIM g=%.1f L^2 delta/M =' % g, ' '.join('%.4f' % x for x in vals))
    if g > 1:
        check('5 TFIM g=1.6 defect decays with L (massive modes suppressed)', vals[-1] < 0.5 * vals[0])

# 6. Heisenberg: recovery trace distance and CMI scale as 1/L^2
hb = F.heis_bond(); d = 2
for L in [8, 12]:
    w, v, E0 = F.ground_states(hb, L, k=2)
    l = L // 2
    rl = F.ti_marginal(v, d, L, l); rl1 = F.ti_marginal(v, d, L, l + 1)
    rB = F.ptrace_last(rl, d, l)
    KB = np.kron(np.eye(d), F.psd_power(rB, -0.5)); Kl = np.kron(np.eye(d), F.psd_power(rl, 0.5))
    Xr = Kl @ np.kron(KB @ rl @ KB, np.eye(d)) @ Kl
    td = float(np.abs(sla.eigvalsh(Xr - rl1, driver='evr')).sum())
    cmi = 2 * S(rl) - S(rB) - S(rl1)
    dl = F.petz_defect(hb, v, L, l) - E0 / L
    print('   Heis L=%d  L^2 ||X - rho||_1 = %.3f  L^2 CMI = %.3f  (c/3)pi^2 = %.3f  L^2 delta = %.3f <= M L^2 ||.||_1 = %.3f'
          % (L, L * L * td, L * L * cmi, np.pi ** 2 / 3, L * L * dl, 0.5 * L * L * td))
    check('6 Heis L=%d reduction (4.1): delta <= M ||X - rho_{l+1}||_1' % L, dl <= 0.5 * td + 1e-12)

# 7. Decoupled stack: defect is additive and scales with each bond (velocity weighting)
L = 6
h1 = F.xxz_bond(0.0); h2 = 0.1 * F.heis_bond()
# stack bond on (q1 q2)_x (q1 q2)_{x+1}: reorder kron(h1, h2) from (x1, x+1_1, x2, x+1_2)
P = np.zeros((16, 16))
for a in range(2):
    for b in range(2):
        for c in range(2):
            for e in range(2):
                P[((a * 2 + c) * 2 + b) * 2 + e, ((a * 2 + b) * 2 + c) * 2 + e] = 1
hstack = P @ (np.kron(h1, np.eye(4)) + np.kron(np.eye(4), h2)) @ P.T
w1, v1, E1 = F.ground_states(h1, L); w2, v2, E2 = F.ground_states(h2, L)
ws, vs, Es = F.ground_states(hstack, L, k=4)
d1 = F.petz_defect(h1, v1, L, 3) - E1 / L
d2 = F.petz_defect(h2, v2, L, 3) - E2 / L
ds = F.petz_defect(hstack, vs, L, 3) - Es / L
gH = F.ground_states(F.heis_bond(), L)
dH = F.petz_defect(F.heis_bond(), gH[1], L, 3) - gH[2] / L
print('   stack: delta = %.6e, delta_1 + delta_2 = %.6e, delta_2/delta_Heis = %.6f' % (ds, d1 + d2, d2 / dH))
check('7 decoupled stack: Petz defect additive, slow copy weighted by its bond', abs(ds - d1 - d2) < 1e-9 and abs(d2 / dH - 0.1) < 1e-8)

# 8. Lemma 3.1: every recorded bond marginal of the iterated construction equals sigma_2*
w, v, E0 = F.ground_states(F.heis_bond(), 8)
rl = F.ti_marginal(v, 2, 8, 4)
outs = F.petz_periodize(rl, 2, 4, steps=6, tol=-1)
check('8 Lemma 3.1: iterated Petz bond marginals all equal', max(np.abs(o - outs[0]).max() for o in outs) < 1e-10)

# 9. Long-crossover two-band chain (C106), small sizes
vals = []
for eps in [1.0, 0.1, 0.03]:
    hb = F.c106_bond(1.0, eps, 2.0); M = F.M_of_bond(hb)
    for L in [4, 5]:
        w, v, E0 = F.ground_states(hb, L, k=4)
        vals.append(min(L ** 2 * (F.petz_defect(hb, v, L, l) - E0 / L) / M for l in range(2, L)))
print('   C106 L^2 delta/M (eps = 1, 0.1, 0.03; L = 4, 5):', ' '.join('%.3f' % x for x in vals))
check('9 C106 chain: defect uniformly bounded in eps (< 4)', max(vals) < 4)

if LONG:
    # 10. SU(3) L = 12, fast/slow ladder, C106 at L = 6 (logged L = 8 values in logs/ took longer)
    hb = F.swap_bond(3)
    w, v, E0 = F.ground_states(hb, 12, k=3)
    x = 144 * (F.petz_defect(hb, v, 12, 6) - E0 / 12)
    print('   SU(3) L=12 l=6 L^2 delta/M = %.4f' % x)
    check('10 SU(3) L=12 below L=9 value 5.387', x < 5.387)
    vals = []
    for eps in [1.0, 0.03, 0.01]:
        hb = F.ladder_bond(eps, 0.5); M = F.M_of_bond(hb)
        for L in [4, 6]:
            w, v, E0 = F.ground_states(hb, L, k=4)
            vals.append(min(L ** 2 * (F.petz_defect(hb, v, L, l) - E0 / L) / M for l in range(2, min(L, 5))))
    print('   ladder L^2 delta/M (eps = 1, 0.03, 0.01; L = 4, 6):', ' '.join('%.3f' % x for x in vals))
    check('11 ladder bounded uniformly in eps', max(vals) < 3.5)
    hb = F.c106_bond(1.0, 0.03, 2.0); M = F.M_of_bond(hb)
    w, v, E0 = F.ground_states(hb, 6, k=4)
    x = min(36 * (F.petz_defect(hb, v, 6, l) - E0 / 6) / M for l in range(2, 5))
    print('   C106 eps=0.03 L=6 L^2 delta/M = %.3f' % x)
    check('12 C106 eps=0.03 L=6 bounded', x < 3)

print('\nALL PASS' if not fails else '\nFAILED: %s' % fails)
