# Constant-history comparison with dependent normalizations retained

Status: **[proved transfer estimates; the precise reference-scheme input
is stated, not silently identified with a different published flow]**.
The ordinary endpoint/localization bound is being checked separately.
This note isolates the other ingredient needed for the F_h estimates in
GN_MARK_CALCULUS.

Author: OpenAI GPT-6-family Codex vacancy worker, 2026-10-28. AI-assisted;
no human review.

## 1. A decaying forcing survives the normalization fixed point

Let two raw normalization maps, labelled L and R, obey the contraction
hypotheses of NORMALIZATION_RESOLVENT on a common admissible history
domain, with suitable invariant balls. For this comparison alone those
balls need not both be O(eta^2): an order-t reference UV endpoint may
require an O(t) ball. The fixed small raw ratio domain still contains
both balls, and the following proof uses only the small z Jacobian.
Let their solved increments be z^L and z^R. Suppose at the same independent
z their difference obeys

    |Phi_h^L(x,z)-Phi_h^R(x,z)|<=A 2^(b h),               (1)

and their z Jacobians are bounded by qW_a, with 0<b<a. The constants and
the exponent are independent of the infrared cutoff. Then

    |z_h^L-z_h^R|
      <=A 2^(b h)+q sum_(k=h)^0 2^[-a(k-h)]|z_k^L-z_k^R|.

In the weighted norm `||f||_b=sup_h 2^(-b h)|f_h|`,

    ||W_a||_b<=1/[1-2^[-(a-b)]].                         (2)

Hence, if q/[1-2^[-(a-b)]]<=1/2,

    |z_h^L-z_h^R|<=2A 2^(b h).                          (3)

This is also an immediate consequence of the exact inverse in
NORMALIZATION_RESOLVENT. It applies to temporal and spatial derivatives
together. No convergence of Z_h^L or Z_h^R is required.

Suppose the raw beta maps differ at equal z by at most A_B 2^(b h)
and satisfy `|partial_(z_k) B_h^L|<=C W_a(h,k)`. Their reduced versions
therefore satisfy

    |B_h^{L,red}(x)-B_h^{R,red}(x)|
                         <=C'(A_B+A) 2^(b h).            (4)

In particular, the O(1) derivative of the velocity beta with respect to
the spatial normalization increment causes no loss of ultraviolet decay.

## 2. What the lattice/reference tree difference must contain

At fixed independent normalization increments and running histories,
compare tree values using the same branch labels, localization convention
and shell decomposition. There are two types of differences.

* A changed initial irrelevant endpoint reaches the fixed initial scale.
  The root-to-initial path therefore supplies an ultraviolet gain.
* For trees with common local endpoints, telescope one covariance or
  normalization-compatible curvature factor. The marked covariance
  lemma bounds this difference by a positive power of its shell scale.
  Allocate part of the root-to-mark path gain to move that factor to
  the root; sum the mark location with the remaining dimensional reserve.

The same operations apply to derivative self-energy roots and hence to
the normalization maps themselves. Thus, **once the ordinary endpoint
and zero-routing summations are justified**, these two operations give
(1) and its beta analogue, with a smaller fixed b>0. Equation (4), rather
than replacing all normalizations by one, completes the comparison.

Let t denote the magnitude of the constant local quartic vector and
let the physical initial irrelevant endpoints have size epsilon^2.
The comparison amplitude is safely

    A+A_B<=C(t+epsilon^2).                               (5)

It is not always C epsilon^2. The published reference model may have
ultraviolet initial kernels of order t, even when the replacement action
has the stronger epsilon^2 endpoint budget. Such reference kernels must
be charged to the ultraviolet difference. Their initial support still
gives the same decay. They must not be reclassified as physical residual
endpoints of order epsilon^2.

## 3. The exact reference input needed for full-history freezing

For the coordinates in GN_MARK_CALCULUS, u contains both quartic
couplings and velocities. A sufficient reference statement is

    |B_h^{R,red}(u,u,...;nu=0)|<=C t 2^(alpha h),         (6)

uniformly in the chosen finite velocity domain and common cutoff scheme.
The relevant component is zero in the inversion-symmetric reference
model; forward quartic components can obey the stronger O(t^2) bound.
Together, (4)--(6) give

    |F_h(u)|<=C(t+epsilon^2)2^(b h),                      (7)

after decreasing b below alpha if necessary. At the actual initial
quartic size t=O(epsilon), this gives the required
`|F_h(g_0)|<=C epsilon 2^(b h)`.

There is a precision in applying the primary source. In
[Mastropietro--Porta, (8.60)--(8.61), Theorem 8.5](https://arxiv.org/html/2108.13708v3#S8.SS2),
the beta functions are written as functions of quartic histories after
the velocity and normalization recursion has been solved. Its
constant-history assertion freezes the quartic history. It does not
explicitly state (6) for an independently frozen full velocity history.
The two statements must be connected, not identified by notation.

One route is to prove (6) within the same reference Ward-identity
calculation using formal constant velocity coordinates. Another is the
trajectory comparison in the next section. Either route must use a
common admissible cutoff convention; changing a multiscale
decomposition is not automatically a single ultraviolet endpoint.

## 4. A sufficient velocity-trajectory conversion

Suppose the published reference statement is available in the chosen
common scheme, uniformly in an initial velocity vector v_0 in a compact
neighborhood. Freeze only the quartics at s, |s|=t. Let its solved
reference velocity trajectory satisfy

    v_(k-1)^R=v_k^R+beta_k^{v,R},
    |beta_k^{v,R}|<=C t 2^(alpha k).                     (8)

Then its limit v_*^R exists, with

    |v_k^R-v_*^R|<=C' t 2^(alpha k),
    |v_*^R-v_0|<=C't.                                   (9)

If this map from v_0 to v_*^R is continuous, Brouwer selects an initial
v_0 for which v_*^R equals any prescribed velocity u_v in a smaller
compact neighborhood: the map `v_0 -> u_v-[v_*^R(v_0)-v_0]` preserves
a C''t cube about u_v. Fixed-cutoff continuity and the uniform summable
tail in (8) already imply the needed continuity; no velocity derivative
or uniqueness is required.

If the reference's ultraviolet matching uses bare quartics lambda_b
with `lambda_0=lambda_b+O(|lambda_b|^2)`, its running initial value must
also be matched to s. This can be done in the same finite-dimensional
Brouwer argument. Use the continuous parameter map

    (v_0,lambda_b) -> (v_*^R,lambda_0),

whose displacement from the identity is O(t) in velocity and O(t^2)
in quartics. The map equal to `(u_v,s)` minus those displacements
preserves a rectangular box of radii C_v t and C_lambda t^2 about
`(u_v,s)`, after enlarging the constants and taking t small. Its fixed
point matches both coordinates exactly. For t=0 the free model already
does so. The reference initial irrelevant endpoints then retain their
actual dependence on the selected bare parameters; Section 2 charges
their difference to the ultraviolet terms. This is a parameter choice
inside the comparison theory, not a change to the physical Hamiltonian.

Compare the reference beta on that trajectory with the beta at the
fully constant history (s,u_v). A marked velocity-history bound of the
form proved from GN_MARK_CALCULUS and NORMALIZATION_RESOLVENT gives

    |delta B_h^R|
      <=C eta sum_(k=h)^0 W_a(h,k)|v_k^R-u_v|
      <=C' eta t 2^(alpha' h),                          (10)

for any 0<alpha'<min(alpha,a). If an ultraviolet-decaying history term
is present, its convolution is also summable and satisfies the same
conclusion after decreasing alpha'. Equations (8) and (10) establish
(6) in that common scheme. This conversion uses the actual anomalous
reference normalization in both beta evaluations.

This is a proved sufficient conversion lemma. Verifying its reference
flow input for the fixed shell convention actually chosen remains part
of the source-to-scheme check; it is not supplied by a general assertion
that renormalization schemes are interchangeable.

## 5. Derivatives of the decaying constant-history term

For quartic coordinates the reference theorem explicitly permits complex
couplings. If (7) holds on an eta-sized complex neighborhood, where
t<=C eta and epsilon^2<=eta, Cauchy's estimate gives

    ||D_(u_quartic)F_h||<=C 2^(b h).                     (11)

The normalization comparison was also analytic, so its derivative is
covered by the same argument. This spends an eta radius once.

If the reference decay has only been established for real velocities,
one need not assert a complex-velocity Ward theorem to obtain a real
velocity derivative. Suppose the unmarked analytic tree bound on a
fixed velocity neighborhood gives `||D_v^2 F_h||<=C eta`, uniformly in
h, while (7) gives `|F_h|<=C eta 2^(b h)` on its real section. For an
inner real velocity domain, the symmetric finite-difference estimate is

    |D_v F_h|<=C eta 2^(b h)/rho+C eta rho.              (12)

Choose rho proportional to 2^(b h/2), bounded by the available fixed
velocity margin. This yields

    ||D_v F_h||<=C' eta 2^(b h/2)<=C'2^(b h/2).          (13)

Only a fixed loss of decay exponent occurs. Select theta small enough
that this final exponent is at least 2theta. Thus a real reference
velocity input plus ordinary analytic tree bounds suffices for the
F_h derivative required in the real contraction argument.

More generally, this can be done for every real common-history coordinate
without any complex reference decay input. A uniform raw C2 bound
`||D_u^2 F_h||<=C` and the real bound (7) give, with the same step choice,
`||D_u F_h||<=C'2^(b h/2)`. This is sufficient for G1/G2. The real C2
normalization alternative is proved in NORMALIZATION_RESOLVENT, so a
smooth velocity-dependent cutoff need not be given a holomorphic
continuation at any stage of this argument.

## 6. Dependency boundary

Equations (1)--(4), (8)--(10), and (12)--(13) are complete elementary
transfer proofs. The tree telescoping in Section 2 still depends on
the all-degree endpoint/localization majorant. The source-to-scheme
statement in Section 3 or 4 must also be justified before labeling the
full ordinary GN bridge proved. Current/stress sources and a conformal
scaling limit are additional conclusions, not consequences of these
beta estimates alone.
