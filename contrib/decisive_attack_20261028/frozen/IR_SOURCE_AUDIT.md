# Independent audit of the infrared source comparison

Reviewed: [`exact/IR_SOURCE_REVIEW.md`](../exact/IR_SOURCE_REVIEW.md),
particularly sections 3-8, together with the initial branch-source and
finite-cell matching statements. Reviewer: Codex GPT-6-family frozen-sector
agent, 2026-10-28. Read-only audit of those files; no human review.

## Verdict

**[checked conditional statements]** The normalized source recursion,
Hölder remainder summation, two-ray resolvent factorization, and the
common-velocity CFT implication with Gaussian four-point functions are
correct under their stated additional hypotheses. I found no counterexample
to those conditional claims. The text correctly keeps the missing marked
RG estimates and distributional/CFT hypotheses separate from the finite-cell
matching output.

In particular, neither the action bounds alone nor the new frozen n=23
matching statement proves hypothesis (S), the two-source shell estimates
(6), a common-velocity conformal limit, or Gaussian four-point distributions.
The distinction is substantive. The certified n=23 band velocities are not
common, so section 8 is a conditional interpretation theorem for a narrower
setting, not yet an application to the candidate.

## 1. The quantity delta and Gaussian source preprocessing

The distinction in section 3 is necessary and correct. A finite-range
branch filter has fixed nonzero momentum derivatives even as microscopic
interaction couplings vanish. Its raw irrelevant Taylor norm therefore
cannot be put into a parameter that is asserted to vanish with lambda.

At first order in a source, write its Gaussian quadratic vertex using
low and high fields. Integrating the high fields leaves the low-low
bilinear unchanged, kills the two mixed terms, and contracts the high-high
term to a scalar. Thus a change in the localized two-external-fermion
source vertex due to its irrelevant quadratic remainder needs an
interaction endpoint. With a marked-tree norm estimate which retains that
fact, the proposed decomposition

    delta = delta_matching + C_probe epsilon

has the right small-coupling counting. This is not, by itself, the marked-tree
estimate: its constants, localization, memory gains, and field-normalization
terms still need proof. The text says so explicitly.

Fermi-center shifts can alter the bare cardinal values and must either be
included in U_0 or bounded in preprocessing. The text does this through its
reference to the explicit branch-source estimate. No derivative in Fermi
centers is required by the subsequent sequence argument itself.

The convention `chi_b=sqrt(Z_b) psi_b` gives
`psi_b^dagger psi_b=chi_b^dagger chi_b/Z_b`, hence `U_ab=Y_ab/Z_b`.
This is correct. An RG scheme already using normalized fields must not
apply this division a second time. The finite reference normalization
Zbar in the covariance is correctly distinguished from running Z_h.

## 2. Source-flow theorem

The recurrence `U_(h-1)=U_h+d_h` gives

    U_h=U_0+sum_(j=h+1)^0 d_j.

The sign and endpoints in the proof are correct. With
`D=sup_(h<=0) 2^(-theta h)||d_h||`, the history-difference term obeys

    sum_(k=h)^0 2^[-4theta(k-h)] ||U_k-U_h||
      <= D 2^(theta h)
          [sum_(r>=1)2^(-3theta r)]/[1-2^(-4theta)].

This follows by exchanging the positive finite sums before taking any
cutoff limit. Together with the maximum-history estimate it gives exactly
the stated absorption inequality. Small epsilon makes its D coefficient
less than one. Absolute summability gives U_* and the tail estimate;
Weyl's singular-value inequality gives invertibility when the displacement
is below sigma_min(U_0).

The proof needs (S) for canonically normalized vertices, uniformly in the
cutoffs. Bounded running action couplings and finite-scale invertibility
are insufficient. The example U_h=2^(epsilon h)I is a valid counterexample
to that weaker inference. No unproved cancellation is hidden in the
sequence algebra; the cancellation is explicitly an input to (S).

## 3. Two-source comparison and contact terms

The shell assumptions imply

    ||E_h(p)-E_h(0)||
       <= const min(2^(gamma h), |p|2^((gamma-1)h)).

Splitting at 2^h comparable to |p| gives O(|p|^gamma), including the
zero-momentum sum and the Lipschitz ultraviolet term. This proves (7)
under (6). A needlessly small ultraviolet constant is not required.

The unknown constant contact matrix B_0 is retained and then removed by
the ray difference. This matters: simply reading the rank of one raw
susceptibility matrix could count a contact term rather than gapless
responses. The stated hypotheses exclude extra direction-dependent
order-one errors; such errors would not be removable as B_0 and would
invalidate the inference if (6) were weakened accordingly.

## 4. The two-ray identity and full rank

For `F(A)=(I+AL)^(-1)A`, direct multiplication on both sides verifies

    F(A)-F(B)=(I+AL)^(-1)(A-B)(I+LB)^(-1).

This identity does not require A, B or L to commute. Its use with the
normalization Zbar on the far right of (10) is correct.

For real Euclidean momentum and real nonzero v,
`|(-ip0+vp1)/(ip0+vp1)|=1`. Therefore the sufficient inverse bound using
`||Lambda_Z||/(4 pi min|v|)<1` is valid in spectral norm. For the stated
rays, independent substitution gives

    B_b(1,xi)-B_b(1,0)
      =v_b xi/[2 pi |v_b|(i+v_b xi)].

Every diagonal entry is nonzero when xi and v_b are nonzero. Invertibility
of all the factors, including both source matrices, proves full rank.
This remains true with both positive and negative branch velocities. A
quantitative smallest-singular-value bound follows from the product
inequality and does not require normality of the forward-coupling matrix.

The result counts independent directional scaling responses. It does not
by itself identify a current OPE, positivity, or a stress tensor.

## 5. Common-velocity CFT inference

Under the explicitly added local compact unitary CFT hypotheses, the
bounded-susceptibility argument is sound. A nonconstant scalar primary
with 0<Delta<1 has a positive Gram contribution proportional to
`|p|^(2Delta-2)`; a dimension-one scalar contributes a positive Gram
coefficient times `-log|p|`. The full Hermitian probe matrix detects any
nonzero such coefficient in some linear combination. Fixed analytic
contacts cannot cancel it. Integer spin and nonnegative weights leave
only (1,0) and (0,1) at dimension one, and their opposite derivatives are
null in a unitary theory. Higher-dimension terms vanish in the ray
difference under the assumed controlled scaling expansion.

Here the connected dimension-zero sector is understood in the usual
extremal/clustering CFT vacuum sense. If "pure vacuum" were instead meant
to allow a coherent superposition of locally distinguishable vacuum
sectors, a dimension-zero classical contribution at exactly zero momentum
would need separate removal. This is a convention clarification, not a
counterexample to the ordinary pure CFT-vacuum interpretation used in the
note. It is especially relevant when discussing frozen ordered sectors.

Two-point rank gives only `rank(D)<=r_R+r_L`. The SU(2)_1 example correctly
shows that these dimensions are not an Abelian rank. Vanishing three-point
functions within a visible two-dimensional subspace also cannot rule out
a bracket into an invisible third direction.

The four-point repair handles that exact issue. Let epsilon=z1-z2 and
eta=z3-z4, with z2-z4 fixed and nonzero. The double simple-pole coefficient
in the current four-point OPE is

    -k([A,B],[A,B]) / [epsilon eta (z2-z4)^2].

Within-pair Wick contractions give only epsilon^-2 eta^-2. Cross-pair
Wick contractions are analytic in both epsilon and eta near zero.
Therefore a Wick four-point distribution has no epsilon^-1 eta^-1
coefficient. Positivity forces the full bracket, including its component
outside the visible span, to be null. Removing null fields proves that
the visible current algebra in each chirality is Abelian.

The positive-level Abelian Sugawara tensor has central charge equal to
its rank. The unitary orthogonal stress sector has nonnegative central
charge. This proves c_R>=r_R and c_L>=r_L, hence the stated rank/2 bound.
The final c_R=c_L identification is explicitly restricted to the
nonanomalous common-velocity case.

Crucially, the note assumes Gaussian four-point distributions and their
passage to chiral projections. Fixed-ray nonexceptional momentum limits
alone do not imply those collision limits or exclude additional contact
and exceptional-momentum contributions. The text correctly preserves
this additional obligation. Its higher-point scaling exponent in section
7 is dimensionally consistent: an m-point Fourier current cumulant with
the overall momentum delta removed has degree 2-m.

## Consequence for the ongoing proof attempt

The note is safe to use as a conditional bridge. Its concrete remaining
inputs are the normalized one-source estimate (S), the two-source shell
bounds (6), and adequate local distributional/current-algebra control.
For the actual unequal-velocity n=23 model, the conformal interpretation
must also be established separately for its mode/sector decomposition.
No step in the checked algebra converts these outstanding inputs into
proved properties of the microscopic chain.
