"""C191 frame check: two hard-core bosons (= the qubit-chain particles) on a ring of L sites, total momentum K = pi.
Relative-coordinate ED, bosonic sector f(L-r) = e^{-iKr} f(r).  Compare levels near the n = 3/16 Fermi kinematics
with the reflectionless Bethe-Yang prediction p_1 L = 2 pi j + tau(E) using (i) the boson-frame transmission phase
tau = arg t and (ii) the fermion-frame alternative tau + pi."""
import numpy as np
from scat191 import Model, unlike_amplitudes
m = Model(w=1j)
ka, kb = m.valleys()[:2]
L = 600
K = np.pi
dim = L - 3
rs = np.arange(2, L - 1)
H = np.zeros((dim, dim), complex)
def Vr(r): return m.V(r) + m.V(L - r)
for i, r in enumerate(rs):
    H[i, i] += Vr(r)
    for d, ad in m.a.items():
        for rp, coef in ((r + d, ad * np.exp(-1j * K * d)), (r - d, ad)):
            if 2 <= rp <= L - 2:
                H[i, rp - 2] += coef
ev, vec = np.linalg.eigh(H)
# bosonic sector test
phase = np.exp(-1j * K * rs)
sym = []
for k in range(dim):
    f = vec[:, k]
    fr = f[(L - rs) - 2]
    sym.append(np.vdot(fr, phase * f).real)
sym = np.array(sym)
bos = ev[sym > 0.5]
print("hermiticity", np.abs(H - H.conj().T).max(), " #bosonic", len(bos), " #fermionic", (sym < -0.5).sum())
# Bethe-Yang predictions: unlike pair with the a-particle at ka + p (p>0 moving right) and b at kb - p.
# reflectionless quantization (p_a) L = 2 pi j + tau  (mod 2 pi), E = eps(ka+p) + eps(kb-p)
def tau_of(p):
    # tau for in_ab with aR at ka+p (any p>0): use unlike_amplitudes with pR = p and pL from energy match
    E = m.eps(ka + p) + m.eps(kb - p)
    # find pL<0 with eps(ka+pL) = eps(ka+p)
    lo, hi = -p - 0.2, -1e-9
    f = lambda x: m.eps(ka + x) - m.eps(ka + p)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.sign(f(mid)) == np.sign(f(lo)): lo = mid
        else: hi = mid
    U = unlike_amplitudes(m, ka, kb, 0.5 * (lo + hi), p)
    return np.angle(U['t']), abs(U['r'])
target = 0.3062  # pR at n = 3/16
Es = []
for frame, shift in (("boson", 0.0), ("fermion", np.pi)):
    preds = []
    for j in range(-3, 4):
        # solve (ka + p) L = 2 pi J + tau(p) + shift for p near target
        J0 = np.round(((ka + target) * L) / (2 * np.pi))
        J = J0 + j
        p = target
        for _ in range(60):
            tau, r = tau_of(p)
            p = (2 * np.pi * J + tau + shift) / L - ka
        preds.append(m.eps(ka + p) + m.eps(kb - p))
    print(frame, "frame predicted levels:", np.round(sorted(preds), 6))
E0 = m.eps(ka + target) + m.eps(kb - target)
near = bos[np.abs(bos - E0) < 0.12]
print("ED bosonic levels near E_F:", np.round(near, 6))
print("free spacing ~", 2*np.pi*m.deps(ka+target)/L)

print("\n--- both families, both signs of the phase, both frames ---")
def amps(p):
    """return (tau for in_ab with aR at ka+|p|, tau' for in_ba), for |p| the right Fermi offset"""
    lo, hi = -abs(p) - 0.2, -1e-9
    f = lambda x: m.eps(ka + x) - m.eps(ka + abs(p))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.sign(f(mid)) == np.sign(f(lo)): lo = mid
        else: hi = mid
    pL = 0.5 * (lo + hi)
    U = unlike_amplitudes(m, ka, kb, pL, abs(p))
    return np.angle(U['t']), np.angle(U['tp']), pL
def family(sgn_tau, shift, which):
    out = []
    for j in range(-8, 9):
        pref = target if which == 1 else -0.2829
        J = np.round(((ka + pref) * L) / (2 * np.pi)) + j
        p = pref
        for _ in range(80):
            t1, t2, pL = amps(p if which == 1 else 0.3)
            if which == 1:
                tau = t1
            else:
                # family 2: a at ka + pL'; need pL' consistent with energy of some pR: iterate on pR
                tau = t2
            pn = (2 * np.pi * J + sgn_tau * tau + shift) / L - ka
            p = 0.5 * p + 0.5 * pn
        out.append(m.eps(ka + p) + m.eps(kb - p))
    return np.array(out)
def match(pred, ed):
    sel = pred[np.abs(pred - E0) < 0.1]
    return np.mean([np.min(np.abs(ed - x)) for x in sel])
EDn = bos[np.abs(bos - E0) < 0.15]
for which in (1,):
    for frame, shift in (("boson", 0.0), ("fermion", np.pi)):
        for sg in (+1, -1):
            pr = family(sg, shift, which)
            print(f"family {which} {frame:7s} sign {sg:+d}: mean |ED - pred| = {match(pr, EDn):.2e}")
