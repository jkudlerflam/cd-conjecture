# A calibrated pair observable and exhaustive stability test

`diagnostics.py` supplies functions for a future converged measurement.
Its default execution uses analytic controls only. It never opens saved
V4 states or correlations and does not fit their inadequate short windows.

## Connected momentum-pi pair

Choose and record a physical spin-pair operator P_j=b_j b_(j+d), with no
unstated Jordan--Wigner string. A separation-two operator can overlap the
inter-valley pair of the proposed V4 field content, but that overlap must
be checked; a different analytic control can require a different d.
Supply the complex raw correlation and both measured means:

\[
C(r)=\langle P_0^\dagger P_r\rangle
 -\langle P_0\rangle^*\langle P_r\rangle.
\]

`connected_pi_pair` subtracts this part, multiplies by exp(-i pi r), then
applies a normalized binomial filter of radius four. The filter suppresses
other momentum carriers and preserves complex values. It does not prove
that the leading pi harmonic is the desired dressed field; several pair
separations, bulk anchors and momentum filters should agree. Nearby
harmonics require enough spatial resolution to separate them.

`pair_exponent` returns a window estimate and logarithmic fit residual.
It rejects windows with fewer than 24 points or less than a factor-two
distance span. These are only basic estimator conditions. An actual
phase interpretation additionally requires an established bulk scaling
window, length and bond-dimension convergence, stable filtering/anchor
choices and known operator overlap. No automated c=2 verdict is returned.

For the leading bosonic inter-valley pair exp(i sqrt(2) theta_c), the
logarithmic covariance convention of the locking audit gives

\[
\eta_{\rm pair}=a_G,\qquad x_{\rm quartet}=2a_G=2\eta_{\rm pair}.
\]

This relation uses the correct chiral Gaussian covariance and survives
drift and unequal velocities. It does not use a density Ward identity for
the raw microscopic N. If eta_pair<1 and the matching coefficient is
nonzero, the quartet is weakly relevant. If eta_pair>1, the other integer
vertices and the Gaussian stability still need testing.

## Exact physical correlation control

Take two independent ferromagnetic XX legs on even and odd sites at half
filling, then b_j -> exp(i pi j/2)b_j. This is a range-two diagnostic
control, **not** a nearest-neighbour candidate for this project.
For P_j=b_j b_(j+1), the exact infinite-chain correlation is

\[
\langle P_0^\dagger P_r\rangle=(-1)^r
 g_{\lfloor r/2\rfloor}g_{\lceil r/2\rceil},
\quad
g_m=\tfrac12\det[G_{i-j+1}]_{i,j=0}^{m-1},\quad
G_n=2\sin(\pi n/2)/(\pi n),\quad G_0=0,
\]

with g_0=1/2. Wick's theorem for the XX chain gives this Toeplitz
determinant; g_1=1/pi and g_2=2/pi^2 are independently asserted in code.
Each g_m has exponent 1/2, so the inter-leg pair exponent is 1 and the
quartet is marginal at x=2. A within-leg separation-two pair is a different
operator in this control and must not be assigned exponent 1 by analogy.

The calibration adds an explicit nonzero disconnected expectation and a
slower-decaying contaminating carrier at momentum 0.31. It then performs
the same subtraction and filter as the future observable. The resulting
eta_pair estimates are:

| Distance window | Estimate | Exact |
|---|---:|---:|
| 32 to 128 | 1.0011874 | 1 |
| 48 to 160 | 1.0006196 | 1 |
| 64 to 192 | 1.0003833 | 1 |

Changing the disconnected background leaves the extracted correlator
unchanged to 1.1e-17. This calibrates the operator handling and estimator
bias in a long analytic window; it does not manufacture a V4 scaling
window or validate a specific V4 pair separation.

## Full integer lattice, including uncertainty bounds

`lattice_stability(Q)` takes the symmetric positive covariance Q on the
commuting fields (theta_c,phi_s). It evaluates

\[
x(r,s)=2(a_Gr^2+b_Gs^2-2h_Grs)
\]

for **all** nonzero integer pairs by reducing to a finite proven bound.
If lambda is the smallest eigenvalue of Q, the quadratic form is at least
lambda(r^2+s^2), and the best pure trial has value min(a_G,b_G). No shortest
vector can lie outside the circle of squared radius min(a_G,b_G)/lambda.
The implementation uses a slightly larger integer box, with both signs of
the mixed coordinate. This is not a cutoff asserted to be exhaustive.

Optional entrywise errors da,db,dh bound the lower eigenvalue by
lambda-||[[da,dh],[dh,db]]||_2. When this is positive, only vectors inside
the circle of squared radius 1/lambda_lower can possibly violate strict
irrelevance. Inside it, exact worst-case entry bounds give
da*r^2+db*s^2+2dh*|rs|. The routine distinguishes a strict relevant witness,
a robust all-irrelevant conditional lattice, a marginal point and an
undetermined uncertainty interval. These numerical margins are not
interval-arithmetic rounding certificates; inputs extremely close to a
threshold should be treated as marginal or unresolved.

Analytic Gaussian controls include:

- Opposite free drift of magnitude 0.6: the correct covariance is I, so
  both pure vertices have x=2. This rejects the old C180 oscillator-metric
  formula, consistently with C191 and the locking audit.
- Q=diag(1.5,1.2): all generic vertices are strictly irrelevant, including
  stated entry errors (0.05,0.05,0.02).
- Q=[[1.1,0.9],[0.9,1.1]]: pure dimensions are 2.2 but (1,1) has x=0.8.
- Q=[[2,6],[6,18.2]]: pure dimensions are 4 and 36.4 but (3,1) has x=0.4.

These are calibration matrices, **not proposed microscopic Hamiltonians**.
They demonstrate why measuring only the pair exponent or checking only
the two pure vertices is insufficient. A physical application must also
establish the Gaussian regime, the allowed operator lattice at the actual
density/momenta, nonzero amplitudes, positive mode energies, and absence
of extra commensurate vertices. The generic lattice alone does not cover
all special resonances or a symmetry-broken state.
