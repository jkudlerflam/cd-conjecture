# A response calculation that identifies the missing medium input

Everything in sections 1--4 is conditional on a two-component compact
bosonic regime with independently identified valley fields. Neither the
existence of that regime nor the matching constants is assumed as a
conclusion about p153c.

## 1. Complete asymmetric convention

Set `Phi=(phi_A,phi_B,theta_A,theta_B)` and

\[
H_G=\frac1{2\pi}\int dx\, (\partial_x\Phi)^T H(\partial_x\Phi),
\quad [\Phi_i(x),\Phi_j(y)]=\frac{i\pi}{2}\eta_{ij}\operatorname{sgn}(y-x),
\quad \eta=\begin{pmatrix}0&I\\I&0\end{pmatrix}.                 \tag{1}
\]

H is an arbitrary real symmetric positive 4 by 4 matrix: ten coefficients,
including density-density, phase-phase and all density-phase mixing.
No equality of valley densities or left/right speeds is imposed. In the
canonical p153c frame the surviving antiunitary inversion sends each field
to minus itself at -x, so every quadratic gradient product is permitted.

Define the four smooth coordinates and their thermodynamic sources by

\[
X=(\delta\rho_A,\delta\rho_B,w_A,w_B)^T
  =D\partial_x\Phi/\pi,\quad D=\operatorname{diag}(-1,-1,1,1),
\quad \delta e=-h^T X.
\]

Here w is a phase gradient divided by pi, **not** the physical transport
current; the latter also depends on H. Around an equilibrium state,

\[
e(X)=e_0+\frac\pi2 X^T DHDX+\cdots,
\quad \chi_X=\frac{\partial X}{\partial h}=\frac1\pi DH^{-1}D,
\quad \boxed{H=\frac1\pi D\chi_X^{-1}D}.                       \tag{2}
\]

This is a genuine full finite-density response identity, not a two-body
scattering approximation. It requires the thermodynamic static limit.
For a finite exactly number-conserving model, a uniform chemical potential
does not change its state between number-sector crossings. Its finite-ring
zero-temperature Kubo response of N is therefore zero. Take a bulk nonzero-q
static response before q->0, or obtain the appropriately normalized
thermodynamic density/winding curvature; do not invert that finite-ring zero.

Writing `H=[[A,C],[C^T,B]]`, the density curvature at fixed winding is pi A.
If winding is allowed to relax it is
`pi(A-C B^(-1) C^T)`. The corresponding density-relaxed winding stiffness
is `pi(B-C^T A^(-1) C)`. These Schur complements specify what is held fixed
when testing total/relative compressibility and stiffness. Positivity of a
raw microscopic susceptibility alone does not certify either two occupied
valleys or a gapless relative boson.

## 2. The measurable kernel and its matching

For four Hermitian microscopic probe densities O_a, define the bulk static
connected response by the imaginary-time integral

\[
\chi^O_{ab}(q,0)=\frac1L\int_{-\infty}^{\infty}d\tau\,
\langle T_\tau O_a(q,\tau)O_b(-q,0)\rangle_c.                  \tag{3}
\]

The real symmetric bulk q->0 limit is used in (2). Both time orderings / both
momentum signs matter with unequal chiral speeds. Equivalently use Hermitian
sine/cosine probes and their symmetrized static spectral sums. A single
positive-frequency residue multiplied by two is not generally enough.
For open chains use central windows with normalized envelope weight and
many wavelengths inside the window, before varying q toward zero.

At the stated matching scale, suppose
`O=Z X + higher-dimension or high-energy operators`. Then

\[
\chi^O=Z\chi_X Z^T+\chi_{\rm UV},\qquad
\chi_X=Z^{-1}(\chi^O-\chi_{\rm UV})Z^{-T}.                     \tag{4}
\]

**The missing quantitative input is this rank-four matched kernel**, together
with its error and compact field/charge assignment. The high-energy analytic
response chi_UV and normalization matrix Z cannot silently be set to zero
and identity. [DRESSING.md](DRESSING.md) already gives one explicit source
of mismatch for raw N. Measuring (3) without an independent normalization
or low-energy spectral separation does not remove the ambiguity in Z.
The equal-time nonanalytic slope provides a second check:

\[
S_X^{\rm sym}(q)=\frac{|q|}{2\pi}DGD+\text{analytic/subleading terms}.
                                                               \tag{5}
\]

Raw spectra can have analytic backgrounds. A susceptibility is an
energy-denominator-weighted integral, not an equal-time correlator.

## 3. A rank obstruction and a concrete fourth probe

**[proved within the leading projected family]** Its independent coefficient
directions can be written `{N,T1_real,T1_imag,K}`, where K contains the fixed
ratio of the second-neighbour hop, D2 and D3 terms. One linear combination
of these four extensive probes is H itself. Multiplying H by a positive
constant leaves its zero-temperature ground state unchanged. Therefore
the static response to these four directions has an exact null vector
and rank at most three (away from a singular change of ground state).
The same statement follows from homogeneity of the ground-state energy.

This is a limitation of that set of **uniform parameter derivatives**;
responses of local energy at finite q can have subleading terms, which
do not provide the missing leading four-coordinate normalization.
In the full microscopic family there are four shape parameters after scale.
The additional epsilon dependence first changes the leading fixed ratios
through higher-order interactions and number-changing terms. It need not
vanish exactly, but relying on it requires those smaller terms and their
operator matching to be resolved. Changing notation does not create a
fourth leading response channel.

A concrete additional observable is the projected, dressed distance-two
current `J2=i(T2-T2^dagger)` (its overall sign is immaterial for rank).
Measure its cross responses with N and the two T1 components. The original
Hamiltonian remains nearest-neighbour: correlation observables may have
range two, and their response can be obtained from (3) without adding a
permanent term to the candidate. Projection and SW dressing must still be
matched to the physical spin observable.

An **exact calibration** of free spinless fermions, deliberately omitting
the target's hard constraint and interactions, uses the dispersion in
DATA.md at mu_eff=-2.1426.
The four Fermi endpoints are

```
-3.1088838363, -2.6618638333, -0.3433539632, -0.1690836744.
```

At those endpoints the old probe form factors
`{1,2cos k,2sin k,-2cos 2k}` have the exact relation
`-mu + cos(1.53)*2cos k + sin(1.53)*2sin k - 2cos 2k=0`.
Their singular values are **4.0231, 3.7905, 0.8366, about 2e-14**.
Replacing the last probe by `2sin 2k` gives
**4.2550, 2.2567, 1.1893, 0.3671**, hence full rank.

For each endpoint the exact density of states is `1/(2pi |v|)`.
Transforming its response to X recovers H with max error 1.4e-15 and G=I
with error 1.3e-15. The four speeds have magnitudes
**1.73312, 1.54041, 0.62682, 0.65653**, so this calibration tests unequal
left/right speeds as well as the normalization factors. Its densities
(0.07115,0.02774) and matching Z are **not p153c medium estimates**: the
actual raw total density is near 0.215, and the interactions are substantial.
Only the response identity and probe-design lesson transfer. Whether J2
really supplies an adequately conditioned fourth in-medium channel is open.

## 4. Covariance, velocities and errors

Let `B_H=H^(1/2) eta H^(1/2)`. The correct spatial chiral result is

\[
\boxed{G=H^{-1/2}|B_H|H^{-1/2}},\qquad
\langle[\Phi(x)-\Phi(0)][\Phi(x)-\Phi(0)]^T\rangle_{\rm sym}
=G\log(|x|/a_0)+O(1).                                        \tag{6}
\]

It obeys `G eta G=eta`; a vertex `exp(i ell.Phi)` has dimension
`x=ell^T G ell/4`. The eigenvalues of B_H (equivalently eta H) are the
signed propagation coefficients. With (1), `partial_t Phi=eta H partial_x Phi`;
if the Fourier phase is `exp[i(qx-omega t)]`, the slopes omega/q are the
**negatives** of those eigenvalues. This fixes the otherwise conventional
right/left sign. Speed units are energy times one physical lattice spacing,
with hbar=1. Multiply effective-unit speeds by |t2|=1/80 to obtain physical
speeds; no target speed has yet been extracted in either unit.

Do not replace eta by the antisymmetric finite-dimensional oscillator
matrix. `calibrate.py` independently transforms (6) to C191's signed chiral
charges and obtains agreement to 1.1e-14, including mixed covariance entries.
G alone does not set absolute velocities: scaling H changes the velocities
and static response while leaving G unchanged. A positive canonical G has
four independent parameters, but a general positive H has ten.

For a supplied spectral-norm bound `||delta chi_X||<=e_chi<lambda_min(chi_X)`,

\[
\|\delta H\|\le
\frac{e_\chi}{\pi\lambda_{\min}(\chi_X)
                         [\lambda_{\min}(\chi_X)-e_\chi]}.    \tag{7}
\]

`gaussian_error_bounds` further propagates an H error e<m=lambda_min(H)
through (6). Set M=lambda_max(H),

```
dS = e/(sqrt(m)+sqrt(m-e))
dR = dS/sqrt(m*(m-e))
dB = dS*(sqrt(M)+sqrt(M+e))
dAbs = (2M+e)*dB/(2m-e)
dG = dR*((M+e)/sqrt(m-e)+M/sqrt(m)) + dAbs/sqrt(m*(m-e)).
```

Then `||delta G||<=dG`, and ordered signed velocities change by at most dB.
These follow from the positive square-root inequality, the bound
`min|eig B_H|>=m`, and Weyl's inequality. They are deliberately conservative
and do not include input matching errors unless those have already been
included in e_chi. Floating-point evaluation is not interval arithmetic.

## 5. Interacting calibration, not an arbitrary fitted covariance

For each of two independent continuum hard-rod gases, the exact energy
density at density rho and winding gradient w is

\[
e(\rho,w)=\frac{\pi^2\rho^3}{6m(1-a\rho)^2}
          +\frac{\pi^2\rho w^2}{2m}.
\]

It follows that `K=(1-a rho)^2`, `u=pi rho/[m(1-a rho)^2]`,
`H=diag(u_A/K_A,u_B/K_B,u_A K_A,u_B K_B)`.
This equation of state follows by squeezing the available length to L-aN
and filling the free hard-core momenta; its curvature is exact for this
control, not an asymptotic-Bethe-ansatz fit to the target.

Use `(rho_A,rho_B)=(.16,.09)`, `(a_A,a_B)=(.6,1.3)`, `(m_A,m_B)=(1,1.4)`.
Then make the **coordinate change** `theta'=theta+F phi`,
`F=[[0,2],[-2,0]]`. It preserves eta and the integer bosonic vertex lattice.
This creates mixed entries in the known G without inventing a microscopic
p153c covariance. Ten finite-difference curvatures of the exact energy,
with steps .001, .0005, .00025, recover H and G with decreasing errors;
at the last step the maximum errors are 3.01e-7 and 1.28e-6.

A separately declared nonidentity probe matrix and additive analytic
contact are then applied and removed using (4); H is recovered to 1.5e-15.
A norm-1e-6 perturbation of the matched response obeys (7) and the covariance
bounds. This calibrates the proposed extraction, including its ill-conditioning
warnings. It supplies **no target occupation, compressibility or K value**.
