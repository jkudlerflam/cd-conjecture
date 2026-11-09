"""C180: checks for the bosonization of the two-valley gas of C173 (3.2).

(1) Band facts at phi = pi/2: valley positions, k_a + k_b = pi, mass, cubic term.
(2) Enumeration of vertex operators O = exp(i(n_a th_a + n_b th_b)) exp(-i(p_a ph_a + p_b ph_b))
    (fermion convention, p_i = n_i mod 2) whose lattice momentum vanishes identically on an
    open set of (|w|, n), at the Sigma point (k_a + k_b = pi, n_a = n_b) and off it.
(3) Scaling dimensions from the Gaussian ground-state covariance, with and without the
    Sigma-allowed drift terms d(ph_+)d(th_-) and d(th_+)d(ph_-), checked against the
    decoupled formula.
Run: python3 optable.py  (seconds).
"""
import numpy as np, itertools
from scipy.linalg import sqrtm

# ---------- (1) band ----------
def band(k, wa=1.0, phi=np.pi/2):
    return 2*wa*np.cos(k-phi) - 2*np.cos(2*k)

def valleys(wa=1.0, phi=np.pi/2):
    ks = np.linspace(-np.pi, np.pi, 400001)[:-1]
    e = band(ks, wa, phi)
    idx = [i for i in range(len(ks)) if e[i] < e[i-1] and e[i] < e[(i+1) % len(ks)]]
    return sorted(ks[idx]), [e[i] for i in idx]

for wa in (0.5, 1.0, 2.0):
    kv, ev = valleys(wa)
    k0 = np.arccos(-wa/4)                      # d-gauge valley, c_j = i^j d_j
    print(f"|w|={wa}: valleys k={np.round(kv,6)}, eps={np.round(ev,6)}, "
          f"sum mod 2pi = {np.mod(sum(kv), 2*np.pi):.6f} (pi = {np.pi:.6f}); "
          f"predicted pi/2 -+ k0 = {np.round(np.mod(np.array([np.pi/2-k0, np.pi/2+k0])+np.pi, 2*np.pi)-np.pi,6)}, "
          f"eps_min predicted {-2-wa**2/4:.6f}")
    # d-gauge dispersion eps_d(q) = 2|w| cos q + 2 cos 2q around q = k0
    e2 = -2*wa*np.cos(k0) - 8*np.cos(2*k0)
    e3 = 2*wa*np.sin(k0) + 16*np.sin(2*k0)
    print(f"      eps''(k0) = {e2:.4f} (1/m), eps'''(k0) = {e3:.4f}")

# ---------- (2) operator enumeration ----------
def momentum(na, nb, pa, pb, wa, n, sigma=True, delta=0.0):
    """Lattice momentum (c gauge) of O; valley a at pi/2 - k0, b at pi/2 + k0 (+ offset delta
    off the Sigma point, which also unbalances the densities)."""
    k0 = np.arccos(-wa/4)
    ka, kb = np.pi/2 - k0, np.pi/2 + k0 + (0 if sigma else delta)
    rho_a, rho_b = (n/2, n/2) if sigma else (n/2*(1+0.3), n/2*(1-0.3))
    return na*ka + nb*kb + np.pi*(pa*rho_a + pb*rho_b)

def identically_zero(op, sigma=True):
    rng = np.random.default_rng(1)
    for _ in range(40):
        wa, n = rng.uniform(0.6, 1.6), rng.uniform(0.05, 0.3)
        P = momentum(*op, wa, n, sigma, delta=rng.uniform(0.05, 0.2))
        if abs(np.mod(P + np.pi, 2*np.pi) - np.pi) > 1e-9:
            return False
    return True

ops = [op for op in itertools.product(range(-4, 5), repeat=4)
       if (op[2]-op[0]) % 2 == 0 and (op[3]-op[1]) % 2 == 0 and any(op)]
for sigma in (True, False):
    allowed = [op for op in ops if identically_zero(op, sigma)]
    print(f"\nSigma point = {sigma}: momentum-zero vertex operators with |n|,|p| <= 4 "
          f"(n_a, n_b, p_a, p_b):")
    print("   ", allowed)

# ---------- (3) dimensions ----------
def dims(Kp, Km, up, um, lam1=0.0, lam2=0.0):
    """H = (1/2pi) X^T Hm X, X = (dph_+, dph_-, dth_+, dth_-). Returns covariance matrix G with
    Delta(e^{i(m.ph + n.th)}) = (1/4) v^T G v, v = (m_+, m_-, n_+, n_-)."""
    Hm = np.zeros((4, 4))
    Hm[0, 0], Hm[1, 1] = up/Kp, um/Km
    Hm[2, 2], Hm[3, 3] = up*Kp, um*Km
    Hm[0, 3] = Hm[3, 0] = lam1       # dph_+ dth_-
    Hm[2, 1] = Hm[1, 2] = lam2       # dth_+ dph_-
    J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    Hh = np.real(sqrtm(Hm)); Hi = np.linalg.inv(Hh)
    A = Hh @ (1j*J) @ Hh
    w, U = np.linalg.eigh(A)
    absA = (U * np.abs(w)) @ U.conj().T
    return np.real(Hi @ absA @ Hi)

def Delta(G, na, nb, pa, pb):
    # e^{i(n_a th_a + n_b th_b)} e^{-i(p_a ph_a + p_b ph_b)}; th_a = (th_+ + th_-)/sqrt2 etc.
    v = np.array([-(pa+pb), -(pa-pb), (na+nb), (na-nb)])/np.sqrt(2)
    return 0.25*v @ G @ v

def Delta_closed(Kp, Km, na, nb, pa, pb):
    return ((pa+pb)**2*Kp + (pa-pb)**2*Km + (na+nb)**2/Kp + (na-nb)**2/Km)/8

print("\nDimension checks (decoupled covariance vs closed form):")
Kp, Km = 0.35, 1.5
G0 = dims(Kp, Km, 4.0, 0.9)
for op, name in [((0, 0, 2, -2), "cos 2sqrt2 ph_-  (inter-valley backscattering)"),
                 ((2, 2, 0, 0), "e^{2sqrt2 i th_+} (Delta N = 4, two ab pairs)"),
                 ((2, 2, 2, -2), "e^{2sqrt2 i th_+} cos 2sqrt2 ph_-"),
                 ((0, 0, 4, -4), "cos 4sqrt2 ph_-"),
                 ((-1, 1, 1, 1), "valley transfer, phi_+ dressed (resonant only)"),
                 ((-1, 1, 1, -1), "valley transfer, phi_- dressed (resonant only)"),
                 ((1, 1, 1, -1), "ab pair, Delta N = 2 (resonant only)"),
                 ((0, 0, 2, 2), "total-density 4kF (momentum 2 pi n, not in H)")]:
    print(f"   {name:50s} Delta = {Delta(G0,*op):.4f}  closed form {Delta_closed(Kp,Km,*op):.4f}")
print("\nDrift terms lam1 = lam2 = lam (Sigma-allowed), K+ = 0.35, K- = 1.5, u+ = 4.0, u- = 0.9:")
print(f"   stability: (dth+, dph-) block needs lam^2 < (u+K+)(u-/K-) = {4.0*Kp*0.9/Km:.4f}; "
      f"(dph+, dth-) block needs lam^2 < (u+/K+)(u-K-) = {4.0/Kp*0.9*Km:.4f}")
for lam in (0.0, 0.3, 0.6, 0.83, 0.9, 1.2):
    if lam**2 >= 4.0*Kp*0.9/Km:
        print(f"   lam = {lam:4.2f}: quadratic form not positive (pseudospin Schur complement "
              f"u-/K- - lam^2/(u+K+) = {0.9/Km - lam**2/(4.0*Kp):.4f} < 0), valley-chiral instability")
        continue
    G = dims(Kp, Km, 4.0, 0.9, lam, lam)
    print(f"   lam = {lam:4.2f}: Delta(cos 2sqrt2 ph_-) = {Delta(G,0,0,2,-2):.4f}, "
          f"Delta(e^(2sqrt2 i th_+)) = {Delta(G,2,2,0,0):.4f}, "
          f"min eig of H-matrix = {np.linalg.eigvalsh(np.diag([4/Kp,0.9/Km,4*Kp,0.9*Km]) + lam*(np.eye(4)[[0,3]].T@np.eye(4)[[3,0]] + np.eye(4)[[2,1]].T@np.eye(4)[[1,2]])).min():.4f}")
