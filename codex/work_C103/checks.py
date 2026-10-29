"""Small C103 normalization checks. Run only inside work_C103."""
from pathlib import Path
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

out = ["# C103 checks", "", "[Numerical] Formula and normalization checks, not a general proof.", ""]
t1, t3 = 1.0, 0.1
A = 2 * (t1 - 3 * t3) / 3
pred = 7 * math.pi**2 * (t1 - 27 * t3) / 810
out += ["## Blocked free chain", "", f"A = {A:.12f}; predicted L^2 correction = {pred:.12f}.", "",
        "| L | A_ring(L) | L^2(A_ring-A) |", "| :-- | --: | --: |"]
err = 0.0
values = []
for L in [8, 12, 16, 24, 32, 48, 64, 96]:
    N = 3 * L
    k = 2 * math.pi * (np.arange(N) + 0.5) / N
    eps = -2 * t1 * np.cos(k) - 2 * t3 * np.cos(3*k)
    energy = eps[eps < 0].sum()
    exact = -2*t1 / np.sin(math.pi/N) + 2*t3 / np.sin(3*math.pi/N)
    err = max(err, abs(energy-exact))
    kp = 2*math.pi*np.arange(N)/N
    ep = -2*t1*np.cos(kp)-2*t3*np.cos(3*kp)
    periodic = ep[ep < 0].sum()
    periodic_formula = -2*t1/np.tan(math.pi/N)+2*t3/np.tan(3*math.pi/N)
    err = max(err, abs(periodic-periodic_formula))
    assert exact < periodic
    einfty = -2*t1/math.pi + 2*t3/(3*math.pi)
    aring = 6*L/math.pi * (einfty*N-exact)
    values.append(aring)
    out.append(f"| {L} | {aring:.12f} | {L*L*(aring-A):.12f} |")
assert np.all(np.diff(values) > 0) and values[-1] < A
assert abs(96**2*(values[-1]-A)-pred) < 1e-4
assert err < 1e-11
A0 = np.diag([-t1, -t1], 1) + np.diag([-t1, -t1], -1)
B = -t3*np.eye(3)
B[2, 0] -= t1
band_error = 0.0
for K in np.linspace(-math.pi, math.pi, 101):
    bands = np.linalg.eigvalsh(A0 + B*np.exp(1j*K) + B.T*np.exp(-1j*K))
    ks = (K+2*math.pi*np.arange(3))/3
    unfolded = np.sort(-2*t1*np.cos(ks)-2*t3*np.cos(3*ks))
    band_error = max(band_error, np.max(np.abs(bands-unfolded)))
assert band_error < 1e-12
lowerM = np.linalg.svd(B, compute_uv=False).sum()
bond_oneparticle = np.block([[A0/2, B], [B.T, A0/2]])
upperM = np.abs(np.linalg.eigvalsh(bond_oneparticle)).sum()/2
assert A <= 2*lowerM
out += ["", f"Ring energy residual <= {err:.3g}; blocked band residual <= {band_error:.3g}.",
        f"Gauge-optimized M lies in [{lowerM:.12f}, {upperM:.12f}].",
        f"Therefore Q <= {A/(math.pi*lowerM*3):.12f}; the example obeys the free theorem.", ""]

def response(T, t=1.0):
    upper = min(2*t/T, 45.0)
    value = quad(lambda x: x*x*np.sqrt(max(0, 4*t*t-T*T*x*x))*expit(-x)*expit(x),
                 0, upper, epsabs=1e-12, epsrel=1e-12)[0]
    return 6/math.pi**2*value

out += ["## XX thermal current", "", "| T/t | R_J(T)/t | (R_J-2t)t/T^2 |", "| :-- | --: | --: |"]
thermal = []
for T in [0.1, 0.05, 0.025, 0.0125]:
    R = response(T)
    thermal.append(R)
    out.append(f"| {T} | {R:.12f} | {(R-2)/T**2:.12f} |")
assert np.all(np.diff(thermal)>0) and thermal[-1] < 2
thermal_coef = -7*math.pi**2/20
assert abs((thermal[-1]-2)/0.0125**2-thermal_coef) < .003
T = 0.3
k_integral = quad(lambda k: ((-2*np.cos(k))*(2*np.sin(k)))**2
                  *expit(2*np.cos(k)/T)*expit(-2*np.cos(k)/T),
                  -math.pi, math.pi, epsabs=1e-12, epsrel=1e-12)[0]/(2*math.pi)
cross_error = abs(response(T)-3*k_integral/(math.pi*T**3))
assert cross_error < 1e-11
out += ["", f"Predicted coefficient = {thermal_coef:.12f}; independent k integral residual = {cross_error:.3g}.",
        f"R_J(10t)/t = {response(10):.12f}; high-temperature leading value = {3/(2*math.pi*10**3):.12f}.", ""]

out += ["## Evaluations of inherited C101 formulas", "",
        "[Numerical] These evaluate established formulas and do not repeat Bethe ansatz or many-body simulations.", "",
        "| Family | Q |", "| :-- | --: |"]
samples = [("XX", 2/math.pi), ("Heisenberg", 1.0)]
for n, k in [(3,1),(4,1),(2,2),(3,2),(8,2)]:
    hk = sum(1/j for j in range(1,k+1))
    d = math.comb(n+k-1,k)
    q = 2*k*(n*n-1)/(n*(n+k)*hk*math.log2(d))
    samples.append((f"SU({n}), symmetric k={k}", q))
for r in [2,3,4,8]:
    samples.append((f"Sp({2*r}), assumed sectors", (2*(r-1)/r+1/(r+1))/math.log2(2*r)))
samples += [("Hubbard repulsive, upper bound", 3/(2*math.pi)),
            ("Hubbard attractive, upper bound", 1/(2*math.pi)),
            ("t-J q=2, upper bound", 3/(2*math.log2(3))),
            ("t-J q=3, upper bound", 1.0)]
for label, q in samples:
    assert q <= 1+1e-12
    out.append(f"| {label} | {q:.12f} |")
out += ["", "## Optical-moment identity on small XXZ rings", "",
        "[Numerical] Derivatives are taken on the finite-ring ground branch before any size limit.", "",
        "| L | Delta | K | D | R | 2s0^2/s1 |", "| :-- | :-- | --: | --: | --: | --: |"]

def xxz_matrices(L, delta, phi=0.0):
    dim = 2**L
    h = np.zeros((dim, dim), complex)
    j = np.zeros_like(h)
    diam = np.zeros_like(h)
    for state in range(dim):
        for a in range(L):
            b = (a+1) % L
            za, zb = (state>>a)&1, (state>>b)&1
            h[state, state] += delta*(za-.5)*(zb-.5)
            if za != zb:
                other = state ^ (1<<a) ^ (1<<b)
                direction = 1 if za == 0 else -1
                h[other, state] += .5*np.exp(1j*direction*phi)
                j[other, state] += .5j*direction
                diam[other, state] -= .5
    return h, j, diam

moment_error = 0.0
curvature_error = 0.0
for L, delta in [(6, .3), (6, 1.0), (8, 1.0)]:
    h, j, diam = xxz_matrices(L, delta)
    evals, vecs = np.linalg.eigh(h)
    gs = vecs[:, 0]
    gaps = evals[1:]-evals[0]
    weights = np.abs(vecs[:, 1:].conj().T @ (j @ gs))**2/L
    K = np.vdot(gs, diam @ gs).real/L
    R = 2*np.sum(weights/gaps)
    D = K-R
    s0, s1 = weights.sum(), np.dot(weights, gaps)
    Rlower = 2*s0*s0/s1
    comm = j @ (h@j-j@h) - (h@j-j@h) @ j
    s1comm = np.vdot(gs, comm@gs).real/(2*L)
    moment_error = max(moment_error, abs(s1-s1comm))
    step = 2e-4
    eplus = np.linalg.eigvalsh(xxz_matrices(L, delta, step)[0])[0]
    eminus = np.linalg.eigvalsh(xxz_matrices(L, delta, -step)[0])[0]
    finiteD = (eplus+eminus-2*evals[0])/(step**2*L)
    curvature_error = max(curvature_error, abs(finiteD-D))
    assert 0 <= Rlower <= R+1e-12
    out.append(f"| {L} | {delta} | {K:.10f} | {D:.10f} | {R:.10f} | {Rlower:.10f} |")
assert moment_error < 1e-11
assert curvature_error < 1e-6
out += ["", f"Double-commutator moment residual <= {moment_error:.3g}; finite-difference stiffness residual <= {curvature_error:.3g}.",
        "", "[Numerical] All assertions passed. No lattice counterexample to the proposed A bound was found.", ""]
Path("CHECKS.md").write_text("\n".join(out))
print("All C103 checks passed; wrote CHECKS.md")
