"""Operator-norm side of the stiffness budget, no diagonalization of chains.

For each reference bond h on C^d x C^d this computes
  M            bare half-spread (gauge does not lower it for these bonds,
               checked in fable_bound_proof/bond_scale_M.py and
               fable_bound_freefermion/gauge_M.py)
  B_torus      ||h~ - pinch_1(h~)||, the maximal-torus stiffness budget:
               sum over an orthonormal basis of conserved diagonal charges of
               D_alpha <= 2 B_torus (<= 4 M always)
  B_q          ||[q,[q,h~]]|| for a single charge q, and sigma_q^2 M (Bernstein)
  B_G          ||C_G h~||, C_G = sum_alpha ad_{t^alpha}^2 the Casimir
               superoperator on site 1 for the Lie algebra G acting on the site,
               generators orthonormal in tr(t t) = 1 (SU(N)) or the spin
               normalization S^alpha (SU(2), spin S)
and the implied bounds on A for the chains whose CFT content is known.

Also checks the general inequalities
  ||h~ - pinch(h~)|| <= 2 ||h~||        and     ||sum_alpha ad_{q^alpha}^2 h|| <= 4 ||h||
on random Hermitian h, and finds the worst ratio of the first (it exceeds 1).
Only numpy. Run: python casimir_budget.py
"""
import numpy as np
from stiffness_ed import spin_matrices, bond_xxz, bond_swap, bond_bilinear_biquadratic, centre, half_spread, double_commutator

rng = np.random.default_rng(1)


def pinch1(h, d):
    """pinching of the site-1 index: keep only components diagonal on site 1"""
    hh = h.reshape(d, d, d, d).copy()
    for a in range(d):
        for ap in range(d):
            if a != ap:
                hh[ap, :, a, :] = 0
    return hh.reshape(d * d, d * d)


def ad2(Q1, h):
    return Q1 @ Q1 @ h + h @ Q1 @ Q1 - 2 * Q1 @ h @ Q1


def casimir_su2(h, S):
    d = int(round(2 * S + 1))
    out = np.zeros_like(h, dtype=complex)
    for Sa in spin_matrices(S):
        Q1 = np.kron(Sa, np.eye(d))
        out += ad2(Q1, h)
    return out


def su_n_generators(N):
    """orthonormal basis of traceless Hermitian N x N matrices, tr(t t) = 1"""
    gens = []
    for a in range(N):
        for b in range(a + 1, N):
            X = np.zeros((N, N), complex); X[a, b] = X[b, a] = 1 / np.sqrt(2); gens.append(X)
            Y = np.zeros((N, N), complex); Y[a, b] = -1j / np.sqrt(2); Y[b, a] = 1j / np.sqrt(2); gens.append(Y)
    for k in range(1, N):
        v = np.zeros(N); v[:k] = 1; v[k] = -k
        gens.append(np.diag(v / np.linalg.norm(v)).astype(complex))
    return gens


def casimir_sun(h, N):
    out = np.zeros_like(h, dtype=complex)
    for t in su_n_generators(N):
        out += ad2(np.kron(t, np.eye(N)), h)
    return out


def opnorm(X):
    return np.linalg.norm(X, 2)


def report(name, h, d, charge, extra=""):
    ht = centre(h)
    M = half_spread(h)
    Bt = opnorm(ht - pinch1(ht, d))
    line = f"{name}: d = {d}, M = {M:.4f}, ||h~ - pinch h~|| = {Bt:.4f} (ratio to M {Bt/M:.3f})"
    if charge is not None:
        Bq = opnorm(double_commutator(ht, charge))
        sig = charge.max() - charge.min()
        line += f", ||[q,[q,h~]]|| = {Bq:.4f}, sigma_q^2 M = {sig**2*M:.4f}"
    print(line + ("; " + extra if extra else ""))
    return M, Bt


print("Random-bond checks of the two general operator inequalities")
worst_pinch, worst_torus = 0, 0
for d in (2, 3, 4):
    for _ in range(300):
        X = rng.normal(size=(d * d, d * d)) + 1j * rng.normal(size=(d * d, d * d))
        h = X + X.conj().T
        ht = centre(h)
        worst_pinch = max(worst_pinch, opnorm(ht - pinch1(ht, d)) / opnorm(ht))
        # maximal torus, orthonormal traceless diagonal charges
        tot = np.zeros_like(ht)
        for t in su_n_generators(d)[-(d - 1):]:
            tot += ad2(np.kron(t, np.eye(d)), ht)
        worst_torus = max(worst_torus, opnorm(tot) / opnorm(ht))
print(f"  max ||h~ - pinch h~|| / ||h~|| over random bonds = {worst_pinch:.4f} (proved <= 2)")
print(f"  max ||sum_alpha ad_alpha^2 h~|| / ||h~|| over random bonds = {worst_torus:.4f} (proved <= 4; equals 2 x pinch ratio)")
# adversarial: 3 x 3 off-diagonal J - I on site 1 tensor identity gives ratio 4/3
d = 3
O = np.kron(np.ones((3, 3)) - np.eye(3), np.eye(3))
ht = centre(O)
print(f"  explicit bond (J - 1) x 1 at d = 3: ratio = {opnorm(ht - pinch1(ht, d))/opnorm(ht):.4f} (= 4/3)")
print()

print("Reference bonds")
h, q = bond_xxz(1.0)
M, Bt = report("Heisenberg S.S", h, 2, q)
C = casimir_su2(centre(h), 0.5)
print(f"  SU(2) Casimir budget ||C h~|| = {opnorm(C):.4f} (= 2 ||S.S|| = 3 M); sum_{{x,y,z}} D_alpha = 3 v K/pi = 3/4;"
      f" implied A = c v <= (2 pi/(k+2)) ||C h~|| = {2*np.pi/3*opnorm(C):.4f} against A = pi/2 = {np.pi/2:.4f}")

for Delta in (0.5, 0.0, -0.5):
    h, q = bond_xxz(Delta)
    g = np.arccos(Delta); K = np.pi / (2 * (np.pi - g)); v = np.pi * np.sin(g) / (2 * g)
    M, Bt = report(f"XXZ Delta = {Delta}", h, 2, q,
                   f"single-twist bound A <= (pi/K) ||h''|| = {np.pi/K*opnorm(double_commutator(centre(h), q)):.4f} against A = {v:.4f}")

for N in (2, 3, 4, 5, 8):
    h, _ = bond_swap(N)
    gens = su_n_generators(N)
    cart = gens[-(N - 1):]
    M, Bt = report(f"SU({N}) swap P", h, N, np.diag(cart[0]).real)
    tot = sum(ad2(np.kron(t, np.eye(N)), centre(h)) for t in cart)
    C = casimir_sun(centre(h), N)
    A = 2 * np.pi * (N - 1) / N
    print(f"  Cartan budget ||sum ad^2 h~|| = {opnorm(tot):.4f}; sum_Cartan D = A/pi = {A/np.pi:.4f};"
          f" full Casimir ||C_SU(N) h~|| = {opnorm(C):.4f} = 2(N+1); implied A <= pi ||C h~||/(N+1) = {np.pi*opnorm(C)/(N+1):.4f},"
          f" A exact = {A:.4f}, conjecture pi M log2 N = {np.pi*M*np.log2(N):.4f}")

for beta, lab, c, v in ((-1.0, "TB", 1.5, 2 * np.pi), (1.0, "ULS", 2.0, 2 * np.pi / 3)):
    h, q = bond_bilinear_biquadratic(beta)
    M, Bt = report(f"spin-1 {lab} (beta = {beta})", h, 3, q)
    C = casimir_su2(centre(h), 1.0)
    k = 2 if beta < 0 else None
    if k:
        print(f"  SU(2)_2: ||C_SU(2) h~|| = {opnorm(C):.4f}; sum_{{x,y,z}} D = 3 k v/(2 pi) = {3*k*v/(2*np.pi):.4f};"
              f" implied A <= (2 pi/(k+2)) ||C h~|| = {2*np.pi/(k+2)*opnorm(C):.4f}, A exact = c v = {c*v:.4f}, conjecture = {np.pi*M*np.log2(3):.4f}")
    else:
        print(f"  (same as SU(3) swap up to a constant; S^z = sqrt2 x orthonormal Cartan charge)")

# spin-1 bilinear-biquadratic family: how the Casimir budget compares with the
# conjectured bound along the gapless SU(2)_2 -> SU(3)_1 arc theta in [-pi/4, pi/4]
print()
print("spin-1 bilinear-biquadratic cos(theta) S.S + sin(theta) (S.S)^2, Casimir and torus budgets vs M")
for th in np.linspace(-np.pi / 4, np.pi / 4, 7):
    Sx, Sy, Sz = spin_matrices(1.0)
    SS = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)).real
    h = np.cos(th) * SS + np.sin(th) * SS @ SS
    ht = centre(h); M = half_spread(h)
    print(f"  theta = {th:+.4f}: M = {M:.4f}, ||C_SU(2) h~||/M = {opnorm(casimir_su2(ht,1.0))/M:.4f}, "
          f"||[Sz,[Sz,h~]]||/M = {opnorm(double_commutator(ht, np.diag(Sz).real))/M:.4f}, ||h~ - pinch||/M = {opnorm(ht - pinch1(ht,3))/M:.4f}")
