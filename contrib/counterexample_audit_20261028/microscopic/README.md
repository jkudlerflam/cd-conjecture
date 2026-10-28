# Independent microscopic and local-normal-form audit

Status: **verified for the microscopic claims in this audit's scope**. No
false lemma or physical-bond mismatch was found in the literal 23/24-state
construction, the frozen 23-state sector selection, the n=23 quadratic
certificate, or the finite-order local Schrieffer–Wolff reduction. This is
not an independent verdict on the vacancy infrared argument or on the
complete c=24 conclusion; those require the separate charge/defect/RG audit.

This audit was performed on 2026-10-28 by an independent GPT-6-family Codex
worker. Prior AI endorsements were not used as evidence. The physical
term list, implementation, certificate arithmetic and cited primary source
were inspected. The new script constructs a full small-cell calculation
from the literal words, independently of the carrier/source formula being
tested. It also regenerates the large certificate. No human review is
claimed.

## 1. Literal bond and exact code

[checked] The h24 matrix embeds all of h23 on occupied pairs. The distributed
onsite terms are consequently contact weighted after squeezing, exactly as
the campaign states. Replacing them by bare occupied-site terms would make
a different model. The explicit lift is Hermitian and conserves occupied
number, internal charge and marker number. Our moderate-parameter test
also includes the physical chemical-potential term; the recovered occupied
block agrees to 1.5e-18 absolute error after floating-point subtraction.

[checked] Every nonzero two-site transition in the literal kinetic and port
operators starts and ends on a matched pair, preserving the two external
arrows. The one-site conversions preserve both arrows. This is the needed
strong statement that **each** forbidden-bond projector commutes with the
full frozen Hamiltonian, rather than merely the sum of its penalties.
The marker indicator and total charge are conserved as well.

[proved from the literal terms] With the physical data field treated as an
order-lambda-squared perturbation, the lambda=0 vacuum block of a selected
cell is the length-24 path tensored with the identity on 23 data bits. Its
ground-space rank is exactly 2^23. The first vacuum excitation has gap

    omega_24 = 2 cos(pi/25) - 2 cos(2pi/25)
             = 0.04706308... > 0.047.

The excited scalar graph is a unit-weight path plus an edge of magnitude
t_d. Its adjacency norm is at most 2+t_d, including after deleting vertices
or edges. The conversion matrix consists of disjoint two-state blocks of
norm g. Therefore, at every cell length, its energy before the marker
potential is bounded below by

    kappa = Delta - 2 - t_d - g > 0.

At the selected length the excited energy relative to the vacuum ground is
at least 2-t_d-g. Thus the vacuum-walker gap is a valid uniform cell gap
Delta_0. This verifies the tensor-product unperturbed space needed by local
perturbation theory. The small field does not belong in the exactly
degenerate H0; putting it in V2, as C33 does, resolves the apparent wording
ambiguity in statements about the code of the full finite-lambda model.

## 2. Frozen marker-sector selection

[proved/checked] The global-sector comparison in
`counterexample_campaign/route_a/README.md` survives examination of the
literal enlarged alphabet. Its use of exact interacting open-chain ground
energies is essential and valid.

For the homogeneous selected-cell chain, concatenation by a port of norm
B=sqrt(3)|lambda| gives |A_m - m u| <= B, with u the infinite-chain energy per
cell. Cutting the at most 2k port bonds touching k wrong cells yields at
most k selected runs. Consequently the excess over the selected periodic
sector is at least k(gamma_b-3B)-2B, and hence k(gamma_b-5B).

For r illegal bonds, the arrow grammar gives r legal intervals, at most 2r
end fragments, and at most k+r selected runs. Every end fragment has at
most one token and energy at least -2; its nonnegative marker, data and
express potentials do not spoil this bound. Removing at most 2k+4r ports
gives the stated excess bound

    k(gamma_b-5B) + r(K-4-7B).

The cut does not discard half of an onsite term: onsite operators remain
inside the intervals, while the forbidden pair contributes its scalar K
and has no active pair transition.

The all-length minimization follows from strict discrete convexity of
E_l=mu_sel-2cos(pi/(l+1)). For real l>=1 its second derivative is positive,
so the minimum positive margin is attained at a neighbor of b. At the
midpoint selector the two neighboring margins are equal to
(delta_23-delta_24)/48. A new rational Machin-series/cosine-series
calculation gives

    gamma_24 > 7.5667157091258e-5,
    kappa > 0.0157705973710,
    gamma_24 - 5*(2*2^-20) > 6.6130413927195e-5,
    K - 4 - 7*(2*2^-20) > 5.99998664855.

Here 2 bounds sqrt(3) from above. These statements concern the frozen
material chain on compatible rings. Transferring the defect margin to h24
with its dynamic contact projector is a further argument, not a consequence
of this sector proof alone.

## 3. Independent literal-to-C33 normalization test

[numerical] `independent_audit.py` constructs legal cell words using the
literal bond implementation and forms the complete two-cell source matrix
S=V(G tensor G). At zero conversion it constructs A=H_exc tensor I +
I tensor H_exc - 2E0 and the two-cell conversion C. The complete leading
coefficient is calculated as

    X = A^-1 S,
    H2[g^2] = -(CX)^dagger A^-1 (CX).

All code states are retained. Independent explicit Jordan–Wigner matrices
and normalized traces extract the quadratic coefficients and residual;
this does not reuse C33's polynomial-size source contraction. The n=1,2
tests use the physical phases phi_x,0=pi/(2n), phi_d,0=pi and
phi_d,+ = phi_d,- = pi/2.

The literal result matches the separate integer-loop formula as follows:

| Test | n=1 | n=2 |
| --- | ---: | ---: |
| Relative hopping error | 4.2e-16 | 5.2e-16 |
| Absolute onsite-matrix error | 5.5e-20 | 1.4e-20 |
| Nonquadratic norm at order g^2 | 0 | 1.4e-20 |

The n=1 onsite matrix is exactly zero in the reference, so a relative error
is not meaningful there. At n=2 the onsite imaginary entry is 1.0437e-14;
its small residual is dominated by floating-point cancellation in the
literal scalar background. The test checks the phase orientation and the
factors alpha^4, beta^2 and 2 in

    E_unit = 2 alpha_23^4 beta^2 g^2.

The common physical field is lambda^2 E_unit 2^-214. There is no missing
lambda^2 or factor of two in the certificate-to-bond identification.

## 4. The n=23 integer certificate

[checked] In a new temporary workspace, the full-loop matrices were
regenerated at 512 bits and then all rational endpoint and pair-product
tests were rerun. Every saved F, M, EF, EM entry and metadata value agrees
exactly with the original. Runtime was about 26 seconds including the
small literal tests. The results are:

    roots = 46,
    delta_2 > 0.0140193864521,
    delta_4 > 3.15587645409e-7,
    smallest endpoint congruence margin > 1.93052640174e-8,
    field/E_unit = 2^-214,
    certified finite conversion g = 2^-365.

[checked] The following potential failure points were examined rather than
inferred from a successful program exit:

* The Jacobi solve uses the correct common imaginary loop orientation and
  a contraction below 1/2+t_d/4. The complex-modulus error bound 10/S is
  sufficient for the stated product radius, including its 100/S^2 term.
* The onsite matrix is antisymmetrized with the matching summed radius.
  Its stored denominator and the hopping denominator both give the
  physical normalization above.
* The proposed numerical eigenvectors are only used to choose an exact
  rational congruence. Integer positivity of its Gram matrix and a
  diagonal-dominance sign test supply the acceptance decision.
* The Laurent determinant degree is at most 46. The 46 disjoint intervals
  with a one-unit inertia change therefore establish exactly 46 simple
  roots, not merely at least 46 crossings.
* Every unordered pair with replacement is included. The interval-width
  subtraction follows from the Lipschitz bound for 2 atan(t). Adjacent
  product phases give the global minimum circular separation, so pair
  transfer between distinct bands is not omitted.
* With A>=3 and ||C||<=2, the even-conversion remainder is bounded by
  alpha^4 g^4/4 for g<=1/2. Quadratic projection and two bond directions
  give the 4n alpha^4 g^4 Bloch bound. The conditioning trace budget is
  exact and passes at g=2^-365; smaller g improves it when the field is
  rescaled with g^2.

The certificate concerns the quadratic projection of H2. Its acceptance
does not itself discard the nonquadratic g^4 remainder or establish c=23.

## 5. Local Schrieffer–Wolff applicability

[read/checked] The primary source was opened directly:
[Bravyi, DiVincenzo and Loss, Sections 4.1–4.4](https://arxiv.org/pdf/1105.0675).
It permits degenerate onsite ground spaces. Equations (4.23)–(4.25)
construct the exact conjugation by a truncated local generator. Lemma 4.2
bounds coefficient interaction strengths uniformly in volume; connected
support growth makes the discarded series exponentially local for a
sufficiently small fixed lambda. Lemma 4.1 and Corollary 4.1 cover small
locally block-diagonal perturbations. The paper does not supply the final
multicomponent CFT stability statement.

[proved application] Blocking the selected cells gives finite onsite
dimension, rank-2^23 code projectors P_j, gap Delta_0>0, and neighboring
perturbations of norm sqrt(3). The order-lambda-squared field can be
included as a polynomial perturbation. Therefore a fixed-order identity
U H U^dagger=D_m+R_m with local ||R_m||=O(lambda^(m+1)) is appropriate;
an extensive many-body band separation is not a required hypothesis.

The auxiliary-threshold refinement is also sound. Write D_m=H0+W and
replace H0 by H0-(Delta_0/2)sum Q_j. It has the same product ground space
and gap at least Delta_0/2. At sufficiently small lambda the unchanged W
satisfies the block-diagonal stability condition. Its ground energy is
the ground energy E_eff of the P restriction. Consequently

    D_m >= E_eff I + (Delta_0/2)sum Q_j.

This uses an operator comparison and does not require D_m to commute
with each Q_j. Since it commutes with global P, the complement is at least
Delta_0/2 above E_eff. High excited energies within P can still overlap
the complement; that does not invalidate the threshold above the ground.

At fixed g, the second-order logical energy scale is positive. Choosing
lambda small makes the order-lambda-cubed local remainder negligible
relative to lambda^2 times even the smallest retained velocity. The
standard infrared stability inference is a separate assumption, but the
claimed local normal form supplies its advertised microscopic input.

## Reproduction and scope

Run from the repository root with NumPy and SciPy:

```sh
PYTHONPATH=/private/tmp/cd-campaign-python /Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 contrib/counterexample_audit_20261028/microscopic/independent_audit.py
```

The script writes `small_checks.json` and `audit_results.json` only in this
folder; C33 generation occurs in a temporary directory. The numerical
tests verify finite statements and normalization. The analytic arguments
above supply the all-length and local-normal-form checks. Overall status
within this scope is **verified**, with no identified refutation of the
model. A complete verdict on c=24 must additionally assess the dynamical
vacancy matching, defect dressing and allowed continuum operators.
