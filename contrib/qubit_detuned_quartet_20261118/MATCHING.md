# Connected, amputated matching in the interacting medium

All momenta are radians per **physical site**. Energies in the spectral
formulas are divided by `t=Delta^2/(2V)=1/80`; times are therefore in `1/t`.
Status: exact finite-system identities and double-precision evaluations
in the stated leading effective model. No thermodynamic theorem is claimed.

## 1. Fixed reference and asymmetric sources

Let P forbid adjacent occupied sites. In the bosonic frame used by C195,

\[
\begin{split}
K_0={}&P\sum_j\big[(r+i\kappa)b^\dagger_{j+1}b_j+\mathrm{h.c.}
 -(b^\dagger_{j+2}b_j+\mathrm{h.c.})\\
 &\hspace{24mm}-2n_jn_{j+2}-{2\over3}n_jn_{j+3}\big]P-\bar\mu N,\\
r={}&\cos(1.53),\quad\kappa=\sin(1.53),\quad\bar\mu=-2.1426.
\end{split}                                                    \tag{1}
\]

The microscopic values are Delta=1, V=40, mu=0.0482175 and
`bar_mu=(mu-3 Delta^2/V)/t`. The retained interactions are of order one
in (1); replacing their resolvent by a free denominator is unjustified.

Write the leading source as

\[
K(\lambda)=K_0+\lambda(X^\dagger+X),\qquad
X^\dagger=\gamma_2B_2^\dagger+\gamma_3B_3^\dagger,\quad
B_d^\dagger=P\sum_jb_j^\dagger b_{j+d}^\dagger P,
                                                               \tag{2}
\]
\[
\gamma_2=-2\epsilon r=-0.002039250562079552,\qquad
\gamma_3={4\epsilon\over3}=0.03333333333333333,\quad\epsilon=1/40.
\]

The mixed adjacent-pair/hop process gives the physical separation-two
coefficient `-2 Delta Re(w)/V`; the three-pair process gives
`2 Delta^3/(3V^2)` at separation three. Both are order epsilon in K.
Consequently the new coefficient is homogeneous quadratic in these
sources, order `epsilon^2 t = Delta^4/(2V^3)` physically. The formal
parameter lambda counts source insertions and is set to one only after
extracting this coefficient. It is not an artificial physical U(1)
symmetry assumption: K0 is the number-conserving reference, whereas (2)
and the microscopic model preserve only parity.

We compute this particular source-squared term. The canonical expansion
on the whole blockaded space and the continuum Wilson coefficient are
different objects. Earlier direct-number-four path cancellations are
not substituted for the medium response, and no uniform error from the
omitted microscopic operators is assigned at this finite epsilon.

For each L search **all** constrained N and all translation sectors.
The lowest state is nondegenerate on the rings used here. Write it as
`|Omega>`, its energy as E0, and

\[
R=Q(K_0-E_0)^{-1}Q,\qquad Q=1-|\Omega\rangle\langle\Omega|.     \tag{3}
\]

This is a grand-ground-state resolvent, with positive excitation energies.
It is not an excited-sector or crossing denominator. The field is held
fixed as L grows: N=3,3,5 at L=16,20,24, respectively. The smallest global
gaps are 0.00325740839, 0.00851002152 and 0.02894583805 in units t.
Their nonmonotonicity and the finite-density steps are explicit limits
on thermodynamic use.

## 2. The exact four-leg object and its counterterms

Define projected single Fourier operators

\[
a_\nu={1\over\sqrt L}P\sum_j e^{-ik_\nu j}b_jP,\quad
k_A=-\pi+{2\pi\over L},\quad k_B=-{2\pi\over L}.               \tag{4}
\]

They are **probe momenta**, not measured dressed valley centres. L=20,24
places them within trial patches of half-width 1/8 around both the bare
centres and the available finite-window carriers. That trial width is
only a specified momentum projection, not an established matching scale.
The L=16 check falls outside those patches. Four legs can conserve lattice
momentum even when the dressed-centre quartet is detuned: they then sample
small nonzero residual momenta inside those patches.

Introduce static sources through
`K(lambda,J)=K(lambda)+sum_nu(J_nu^* a_nu+J_nu a_nu^dagger)`.
Let W be minus the ground energy, and take its Legendre transform only
with respect to these **two complex modes**. All other modes have already
been integrated into W. This distinction matters: it is not the full
lattice 1PI vertex before eliminating the opposite momenta.

For a general choice of fields the exact parity-symmetric identity is

\[
\Gamma^{(4)}_\lambda
 =-(G_\lambda^{-1})^{\otimes4}W^{(4)}_{\lambda,c},\qquad
W_{4,c}=\langle bbbb\rangle-\sum_{3\ \mathrm{pairings}}
                  \langle bb\rangle\langle bb\rangle .         \tag{5}
\]

Here the indices include creation/annihilation, mode, and frequency, and
G is the *interacting* connected two-point kernel in that same source
space. Parity kills the three-leg terms. In general taking the lambda^2
coefficient of (5) also differentiates G inverse: these are external-leg
counterterms, and anomalous two-point functions must be retained.

The particular projection (4) makes those subtractions exactly tractable:

- No two of the external annihilation momenta sum to zero modulo 2pi.
  Thus all three anomalous Wick products vanish by exact translation,
  even at nonzero lambda. The two-mode G has no anomalous entries; its
  normal off-diagonal entries also vanish.
- Charge counting makes the all-annihilation four-leg function start
  at lambda^2. Normal external-leg corrections start at lambda^2 too,
  so they do not alter its leading coefficient. Full Nambu amputation
  outside the selected source space is not being asserted.
- Every nonempty proper subset of the multiset `XXa_Aa_Aa_Ba_B` has
  nonzero charge or momentum. A charge-neutral proper subset would
  contain one X dagger and two annihilators; its momentum would be
  2kA, kA+kB, or 2kB, none of which is zero. All folded energy terms
  therefore vanish, including the higher-order normalization terms.

Define c6 by the energy monomial

\[
E(\lambda,J)=\cdots+
 c_6\lambda^2(J_A^*)^2(J_B^*)^2+\mathrm{c.c.}+\cdots .          \tag{6}
\]

There are 90 distinct ordered words w of `BBaabb`, where B=X dagger.
The rightmost action is denoted O1. Ordinary nondegenerate perturbation
theory and the preceding selection rule give the exact coefficient

\[
\boxed{c_6=-\sum_{w}\langle\Omega|O_6R O_5R O_4R O_3R O_2R O_1
 |\Omega\rangle.}                                             \tag{7}
\]

All N+/-1, N+/-2, and other intermediate sectors allowed by a word use
the complete interacting K0. This is a finite-density six-insertion
connected response, not a two-to-six particle crossing amplitude.

The static external susceptibilities are

\[
\chi_\nu=\langle a_\nu R a_\nu^\dagger
                    +a_\nu^\dagger R a_\nu\rangle.             \tag{8}
\]

They are dimensionless here; in physical units they equal t times the
physical susceptibilities. The collective Legendre coefficient is
`c6/(chi_A^2 chi_B^2)`. Each repeated source occurs twice in (6), so the
fourth derivative contributes 2!2!, which also occurs when extracting
the effective-action monomial: there is **no extra factor four**.
For comparison to the normalization of a spatially local vertex define

\[
\boxed{\frac{g_{4,L}}t=\frac{L c_6}{\chi_A^2\chi_B^2}.}        \tag{9}
\]

The factor L follows from four Fourier legs normalized by 1/sqrt(L)
and one spatial sum. Equation (9) is a normalization, **not a proof
that this static finite 1PI response is local**.
Throughout, g4,L denotes the coefficient of lambda^2. It is not the full
response at lambda=1, whose analytic branch and higher terms also need
control.

## 3. Source-projected bound and numerical term

For each word set

\[
z_w=O_3R O_2R O_1|\Omega\rangle,\quad
y_w=O_4^\dagger R O_5^\dagger R O_6^\dagger|\Omega\rangle.
\]

Because R is positive,

\[
|c_6|\le\mathcal B_6:={\sum_w}
 \sqrt{\langle z_w|R|z_w\rangle\langle y_w|R|y_w\rangle},\qquad
\frac{|g_{4,L}|}t\le\frac{L\mathcal B_6}{\chi_A^2\chi_B^2}.     \tag{10}
\]

This uses the actual source overlaps and interacting spectra, not a global
operator norm divided by a guessed gap. It remains an exact algebraic
bound at every finite nondegenerate reference and for the projection in
the next section. Its floating-point evaluations are:

| L | chi_A | chi_B | c6 | B6 | `g4,L/t` | bound in (10) |
|---|---|---|---|---|---|---|
|20|45.68918586|149.28130524|-719.47284665|2720.71330151|-0.000309319086|0.001169701615|
|24|48.43918541|63.33489283|-483.87747464|1140.39139659|-0.001233862897|0.002907939935|

The huge unamputated c6 is not itself a huge local interaction; the
external susceptibilities contain soft poles. Conversely amputation does
not certify convergence. The factor-four change in (9) and different
selected densities forbid treating either row as a continuum C4 value.

The calculation also resolves
`c6=gamma2^2 C22+2 gamma2 gamma3 C23+gamma3^2 C33`:

| L | C22 | C23 | C33 |
|---|---|---|---|
|20|189206.543878|-98494.323762|-660284.979456|
|24|385860.905259|-14712.465274|-438734.028563|

The independently evaluated mixed-source response reconstructs this
quadratic expression to below 4e-13 absolute. The source interference is
included rather than inferred from bare pair form factors. Imaginary
remainders in the real-origin convention are numerical roundoff.

## 4. Precisely defined cutoff, and what it does not accomplish

For a diagnostic hard energy cutoff Lambda in units t, define

\[
\Pi_\Lambda=|\Omega\rangle\langle\Omega|
 +\mathbf1_{K_0-E_0\ge\Lambda,\ K_0-E_0>0},\qquad
R_\Lambda=\frac{\mathbf1_{K_0-E_0\ge\Lambda,\ K_0-E_0>0}}
                   {K_0-E_0}.                                \tag{11}
\]

Compress every source and probe by Pi_Lambda, and replace every R,
including those in chi, by R_Lambda. Equations (5)--(10) then describe
the exact response of this *spectrally projected finite model*. Pi commutes
with number and translation, so the selection rules still hold.

| L | Lambda (physical cutoff/t) | chi_A | chi_B | normalized coefficient (9) |
|---|---|---|---|---|
|20|0.25|0.02377920|0.02337350|+12379.4316|
|20|0.50|0.02377920|0.02337350|+10779.5697|
|24|0.25|0.03907458|0.02670584|-13553.7624|
|24|0.50|0.03907458|0.02670584|-8427.8677|

Here the table label denotes the physical energy cutoff divided by t;
Lambda in (11) is that same dimensionless number. Removing the soft
states also removes nearly all external-leg weight. The resulting large
inverse-response coefficients are an explicit warning: a global energy
projection is not a local Wilsonian elimination retaining the low-energy
external fields. They neither establish an instability nor a matching
plateau. Varying this cutoff is a diagnostic, not physical RG running.

## 5. Why the composite pair shortcut fails

For comparison, `medium.py` computes the fourth-energy coefficient of
two X dagger insertions and two
`P_pi=L^(-1/2) sum_j (-1)^j b_j b_(j+2)` insertions. The six `BBPP`
words have no folded term. Dividing by the pair susceptibility squared
gives a legitimate *composite-field* inverse response. It is not (9).

An exact counterexample is a stable quadratic oscillator

\[
H=\omega a^\dagger a+{\lambda g\over2}(a^{\dagger2}+a^2),
\quad P=a^2,\quad |\lambda g|<\omega.
\]

Its first-order anomalous propagator is
`F1(tau)=-g exp(-omega |tau|)/(2 omega)`. Wick's theorem gives
`<P(tau)P(0)>_c=2 lambda^2 F1(tau)^2+...`, even after subtracting
`<P>^2`. With a probe `eta P^dagger+eta^*P`, exact ground-energy
expansion gives

\[
\chi_P=1/\omega,\qquad
[\lambda^2(\eta^*)^2]E=-g^2/(4\omega^3),\qquad
c/\chi_P^2=-g^2/(4\omega)\ne0.                                \tag{12}
\]

Yet W4,c of four single oscillator fields and its amputated vertex vanish
identically because the model is Gaussian. The overlooked terms are
products of anomalous *single-field* propagators, not disconnected
composite means. This explains why the extra four-leg construction is
necessary, even though the composite finite-field check passes.

## 6. The remaining continuum matching step

If two compact bosons are established, use
`b_j ~ sum_nu Z_nu exp(i k_nu j) :exp(i theta_nu):_(a0)+harmonics`.
The desired normalized operator is

\[
V_4=:e^{2i\theta_A+2i\theta_B}:_{a_0},\qquad
H_4=\int dx\,[\lambda_4 e^{iq_4x}V_4+\mathrm{h.c.}],\quad
q_4=\mathrm{wrap}[2(k_A+k_B)].                                 \tag{13}
\]

This requires, at the same fixed physical point:

1. An interacting low-energy projection and a **local** scale-dependent
   matching functional with the chosen static/dynamic order of limits.
   Its connected four-point kernel and full external inverse two-point
   functions obey (5). Source derivatives and analytic UV contact terms
   must be evaluated in the same regulator and subtracted consistently.
2. Dressed single-field residues, harmonic identification, and the
   normal-order/fusion map from the local four-field kernel to V4.
   Single-leg residues alone do not fix that composite matching map.
3. Convergence in size, external residual momenta and frequencies, a0,
   and regulator, with an error for the microscopic elimination and
   every omitted higher source order. Soft reference denominators make
   an error estimate based only on epsilon inadequate.

In a convention where the local kernel can be matched multiplicatively,
write its total matching factor as M4(a0), including the residues and
fusion normalization; do not silently set it to one or to Z_A^2 Z_B^2.
If additional matching work bounds the additive error by e_L in units t,
the finite bound would imply the **conditional** inequality

\[
|\lambda_4|/t\le |M_4|\,L\mathcal B_6/(\chi_A^2\chi_B^2)+e_L.   \tag{14}
\]

No bound on M4 or e_L is supplied by the current rings or saved equal-time
records. Equation (14) is an explicit measurement/error target, not an
assertion that a multiplicative local match exists at these L. If operator
mixing is needed, replace M4 by the corresponding row of the measured
mixing matrix and bound every contributing kernel. This is where the
finite result must be strengthened before it enters a locking length.
