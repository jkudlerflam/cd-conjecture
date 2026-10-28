# Uniform vacancy matching: exact cancellations and an all-order route

Status: **[proved algebraic reductions; proposed controlled RG completion]**.
This note does not claim that a microscopic uniform defect inequality has
already been proved. It identifies an all-order cancellation specific to the
cell construction, gives an exact Gaussian energy inequality, and separates
a gap-free ultraviolet normal form from the remaining infrared estimate.
It is more specific than an assertion that irrelevant operators are harmless.

Author: OpenAI GPT-6-family Codex worker, 2026-10-28. No human review.
The cell-gradient identity was independently found by the parallel
`matching_attack` worker. The Gaussian completion-of-squares bound was
independently found by the coordinating worker. Only this directory was
written by this worker. No numerical experiment is asserted here.

## 1. What must be estimated

Work first at zero port coupling. For a collection of cells with lengths
at most a fixed cutoff, the data space is an exact spectator and the
vacuum walker has a finite-dimensional ground state. Write

    H = H_ch + J sum_(c,r) A_(c,r) h_(c,r),
    sum_r h_(c,r) = H_c,
    H_c |0_c> = e_c |0_c>.

The material terms include the original distributed on-site halves. A bath
bond at a cell boundary can multiply contributions from both neighboring
cells; this does not spoil the cell-local decomposition. The allowed finite
cell list has a positive minimum gap omega_* in unscaled material units.

For finite-cell defects or interfaces, the desired error is a *spatially
weighted connected energy norm*, not the ordinary energy density norm.
If a linked energy interaction Phi(X) is available, a sufficient bound is

    sup_x sum_(X contains x) (1+diam X) ||Phi(X)|| <= J epsilon(s),
    s=J/t,                 epsilon(s) -> 0.                 (1)

After identical translated bulk clusters are canceled, a cluster of diameter
R can meet at most O(R) possible positions relative to one interface. This
is the reason for the extra spatial moment. A rank-one canonical density
term must be separated before (1) is applied.

## 2. An all-order vacuum-excursion identity

**[proved]** Let m be the number of charge contacts assigned to one cell,
including boundary contacts carrying distributed on-site halves. Set

    f_r = h_r-H_c/m,
    B_r = sum_(j<=r) f_j,
    B_0=B_m=0.

Discrete summation by parts gives the exact operator identity

    sum_(r=1)^m A_r h_r
      = Abar H_c + sum_(r=1)^(m-1) (A_r-A_(r+1)) B_r,
    Abar = (sum_r A_r)/m.                                  (2)

In particular, for every excited eigenstate |a> of H_c,

    sum_r <a|h_r|0> = 0,
    <a|sum_r exp(i q r) h_r|0> = O(q).                     (3)

There is a stronger time-dependent statement. For any scalar function a(tau),
the Hamiltonian J a(tau) H_c has mutually commuting values at different times.
It cannot create a vacuum-to-excited transition. In a time-ordered vacuum
matrix element of (2), every excursion out of |0> must therefore begin and
end with the second term of (2). Each excursion has at least two external
charge-difference factors. Insertions of Abar H_c can occur between its ends,
but cannot replace either end. Taking the logarithm to select connected
excursions preserves this property: products of excursions already contain
at least two difference factors each.
The scalar term J sum_r A_r <0|h_r|0> is the linear vacuum profile and is
subtracted before a contribution is called a nontrivial excursion vertex.

This statement is an identity of the unintegrated cell influence functional.
It does **not** mean that every renormalized charge vertex at a later RG scale
still carries two explicit derivatives. Contracting high-energy legs can
generate scalar sources and quadratic charge terms, which must be retained
or resummed in the effective action.

**[Gaussian infrared consequence]** For the smooth density component of the
contact operator, A-alpha is proportional to rho_ch, of scaling dimension
one. A difference in (2) adds one derivative. In the derivative expansion
at external frequencies below J alpha omega_*, the irreducible, uncontracted
n-leg walker vertex therefore has long-wavelength dimension at least n+2,
hence at least four. The corresponding smooth-charge transition spectral
density starts as omega^3 instead of the generic omega of an uncanceled
density coupling. This consequence uses the Gaussian density expansion;
the underlying identities (2)-(3) do not. No assertion about the 23 gapless
payload fields is needed for this zero-port calculation.

If common-mode insertions are expanded around Abar=alpha, the simple
resolvent power count for an ordered n-insertion excursion is

    C^n J^n/(J alpha omega_*)^(n-1)
      = J C^n/(alpha omega_*)^(n-1).                        (4)

Equation (4) counts denominators, not the combinatorial sum of all ordered
excursions or functional derivatives; no convergence theorem is implied.
The powers of 1/J do not by themselves make this large: at a charge momentum
cutoff Q, each smooth density leg has typical amplitude O(Q), and two legs
have an additional O(Q) derivative. This suggests an infrared expansion
parameter controlled by Q/omega_*, at fixed cell list. This is useful power
counting; replacing a typical amplitude by an operator-norm estimate would
be incorrect.

## 3. Exact all-order Gaussian volume bound

**[proved]** Let a finite collection of Gaussian bath modes have canonical
coordinates q_k,p_k and positive frequencies omega_k. Let internal Hermitian
operators O_k commute with the bath but need not commute with one another.
For any internal Hamiltonian H_s,

    H = H_s + sum_k (p_k^2+omega_k^2 q_k^2)/2
              + J sum_k q_k O_k

satisfies

    H >= E_bath + H_s - (J^2/2) sum_k O_k^2/omega_k^2.     (5)

To prove this, put Q_k=q_k+J O_k/omega_k^2 and complete each square. Since
[p_k,Q_k]=-i, each shifted oscillator is bounded below by omega_k/2.
Mutual commutators of Q_k and Q_l are irrelevant to the sum of these
individual operator inequalities. A product of the unshifted bath vacuum
and an internal ground state gives the complementary variational upper bound

    E(H) <= E_bath+E(H_s).                                 (6)

Suppose O_k=sum_j C_kj S_j and the static covariance matrix obeys

    C^T diag(omega_k^-2) C <= (C_0/t) I.

The corresponding operator quadratic form is bounded by
(C_0/t) sum_j S_j^2 even for noncommuting S_j: factor the positive matrix
and express its quadratic form as a sum of Hermitian squares. Thus bounded
local sources give an energy lowering at most C J^2 n/t on n source sites,
with no internal gap assumption. This is suited to the linear-cost long-cell
part of the sector comparison.

It is not a per-interface bound. Applying (5) separately in two sectors
leaves an extensive uncertainty. The exactness of (5) must not be used to
claim that the common extensive contributions have canceled. It is also an
exact statement about a Gaussian bath, not an exact lattice identity for the
full contact projector A.

## 4. A gap-free ultraviolet normal form

There is a way to avoid assuming an all-order microscopic charge cumulant
bound merely to remove high transfer frequencies. The phase-model hopping
has bounded finite-range terms because the phase lowering operator has norm
one, despite its infinite on-site alphabet.

Use dimensionless H_0=H_ch/t, V=sum A_j h_j, and s=J/t. For delta=Lambda/t,
let

    V_low = integral_R g_delta(u) exp(iu H_0) V exp(-iu H_0) du,
    Fourier[g_delta](nu) = exp[-(nu/delta)^2].              (7)

Define an anti-Hermitian filtered inverse-Liouvillian S_delta by its energy
matrix elements

    (S_delta)_ab = [(1-exp[-((E_a-E_b)/delta)^2])
                         /(E_a-E_b)] V_ab.                (8)

The value at zero is zero. Hermiticity of V makes S_delta anti-Hermitian,
and [S_delta,H_0]=-(V-V_low). The exact conjugation therefore has the form

    exp(s S_delta)(H_0+sV)exp(-s S_delta)
           = H_0+s V_low+R_delta.                         (9)

The real-time kernel of (8) is an imaginary multiple of sign(u) times an
erfc function. It has Gaussian tails on time scale delta^-1. A Lieb-Robinson
bound consequently gives a quasi-local S_delta with range O(delta^-1).
The proposed next estimate is an interaction-norm bound on the exact
integral remainder in (9), of order s^2 times a finite inverse power of
delta, with one spatial moment included. It must be derived in a specified
local interaction norm; the finite-volume BCH identity alone does not give
it. The filtered-locality construction needs no many-body gap. The
spectral-flow literature provides its tools, but the gapped-phase theorem
is not being applied to this gapless problem.

A conservative implementation may choose delta=s^a with a sufficiently
small positive a. If the chosen interaction norm gives

    ||R_delta||_weighted <= C s^2 delta^-p,                (10)

then a<1/p gives t||R_delta||_weighted/J -> 0. The best exponent is not
needed for existence. A dimensional estimate giving J^2/Lambda for the
ordinary interaction strength and range t/Lambda would yield the sharper
J^2 t/Lambda^2 weighted estimate, but that estimate should not be substituted
for a specified interaction-norm proof: support multiplicities can add
powers of t/Lambda. In particular, Lambda proportional to J is a poor choice
for this deliberately crude first-moment estimate.

Equations (7)-(9) are exact finite-volume algebra. This note has not fixed
a Banach interaction norm and proved a numerical exponent p in (10).
That is a local-dynamics technical task, and no global spectral gap is
required to perform it.

## 5. The low-frequency step and the bridge cancellation

There is a crucial distinction after (9): a low *energy-transfer-frequency*
operator is not a projection onto low absolute-energy charge states.
In particular,

    ||A_low-alpha|| = O(Lambda/t)

does not follow from the small vacuum spectral weight. A high-energy charge
state can have small-frequency matrix elements. Treating (7) as a bounded
small-amplitude Gaussian field in the full Hilbert space would be a concrete
error. A ground-state RG or a relative-form estimate is still needed.

In the low-energy Gaussian representation there is, however, a useful
all-order organization of the infrared graphs:

1. Resum every static one-charge-leg vacuum blob as a static forward source.
   Its Gaussian energy is an exact quadratic functional of that source.
2. Retain the local quadratic charge terms in the Gaussian action, provided
   their relative stiffness remains small. Their coefficients may vary in
   space according to the cell pattern.
3. Expand the remaining walker vertices and charge descendants/harmonics.
   The walker vertices obey the endpoint identity in Section 2 before their
   high-frequency contractions. The first charge harmonic has dimension
   four at the reference density.

Why step 1 is necessary: a single smooth-density covariance has an absolute
time integral proportional to 1/|x|, so its absolute first spatial moment
does not converge. If a charge line is a bridge separating a connected
vacuum graph into two parts, time-translation invariance of each part makes
that line carry zero frequency. The exact smooth static covariance is local
in space, apart from the separately handled canonical zero mode. Such bridges
must be collapsed using that identity before absolute-value estimates.

The nonlocal remainder has no one-line bridges. A connected graph with V
collapsed spacetime vertices and E smooth covariance lines then has E>=V.
At frequencies below J alpha omega_*, where finite-cell vertices admit a
local derivative expansion, common long-distance power counting in 1+1
dimensions gives at least two powers of infrared decay for its vacuum
energy per unit time and length. The extra power needed to anchor it at an
interface leaves a convergent dyadic sum. The endpoint derivatives and
dimension-four harmonics improve this counting. This argument does not
control the intermediate scales at which the cell vertices remain retarded.
For instance, a two-endpoint covariance obeys schematically

    |K(x,tau)| <= C/(1+x^2+t^2 tau^2)^2,
    sum_x integral d tau (1+|x|)|K(x,tau)| < infinity.       (11)

This bridge argument is a diagrammatic explanation of common-bulk
subtraction at all orders. It is not, by itself, a bound on the combinatorial
sum of all graphs. A constructive RG completion must establish a uniform
small running interaction norm and control overlapping subgraphs; simply
asserting E>=V would not do that.

## 6. A precise sufficient RG output

The following output would complete the finite-cell part. It states a
stronger, concrete intermediate estimate than a claim of phase stability.
Let ell be the charge length scale, measured in units of the matching
cutoff. After forward and quadratic resummation, let W_ell(X) be the
connected vacuum energy assigned to scale ell. It is sufficient to prove

    sup_x sum_(X contains x) ||W_ell(X)||
             <= C J eta(s) ell^(-2-eta_0),
    diam(X) = O(ell) in a summable localization norm,       (12)

with eta(s)->0 and eta_0>=0. Even eta_0=0 is sufficient: multiplication by
the interface anchor ell leaves sum_(ell=2^n) ell^-1 finite. The bound must
hold for every placement of the finitely many cell types, with constants
depending on their fixed norm/gap bounds but not their number or separation.
The Gaussian canonical zero mode is an explicit separate term. At most k
interfaces or altered bounded cells then give

    |Delta E_remainder| <= C J k eta(s).                  (13)

The ultraviolet term (10), port cutting/restoration norms, and the Gaussian
volume bound for long cells can then be combined with (13). One obtains
the desired sector inequality with a small function epsilon(s), plus the
explicit port budget and finite-ring rounding terms.

The all-order endpoint identity supplies a strong reason for (12), and
the spectral filter isolates ultraviolet locality from the charge CFT input.
I have not derived (12) for the exact phase-model contact influence
functional. In particular, the single-impurity theorems below do not imply
the uniform spatially weighted array estimate. Thus the present result is
an explicit route and partial proof, not a completed verification of M.

## 7. Primary sources actually opened

- [Abdesselam, *The Ground State Energy of the Massless Spin-Boson Model*,
  arXiv:1005.4366](https://arxiv.org/pdf/1005.4366), Sections 1, 4, 5.
  It constructs a convergent energy expansion by reorganizing linked graphs;
  the displayed convergence radius depends on temporal kernel norms. This
  is a single two-level impurity and supplies no spatial first-moment array
  estimate. It cannot simply be cited as the desired theorem.
- [Hasler and Herbst, *Ground States in the Spin Boson Model*,
  arXiv:1003.5923](https://arxiv.org/pdf/1003.5923), Theorem 1 and Remark 2.
  Their analytic massless-bath result relies on a particular off-diagonal
  coupling structure. The finite-dimensional extension explicitly requires
  P A P=0 and (1-P)A(1-P)=0. Generic cell bond operators do not satisfy that
  latter condition, so direct application would be invalid.
- [Bachmann, Michalakis, Nachtergaele and Sims,
  *Automorphic Equivalence within Gapped Phases of Quantum Lattice Systems*,
  arXiv:1102.0842](https://arxiv.org/pdf/1102.0842), Sections 2-4.
  The filtered locality machinery accommodates infinite-dimensional on-site
  spaces with suitable bounded interactions. The main ground-state
  equivalence theorem assumes a gap; (7)-(9) only use the filter construction
  and local dynamics, and do not invoke that gapped conclusion.
- [Dimock and Hurd, *Sine-Gordon Revisited*,
  arXiv:math-ph/9907017](https://arxiv.org/pdf/math-ph/9907017),
  Sections 2-4, especially Theorem 18. This is an actual convergent infrared
  polymer-RG construction in the irrelevant-cosine regime. Its polymer norm
  contracts at each blocking step, and it controls the extracted vacuum
  energy and resummed quadratic stiffness. Theorem 18 treats translation
  invariant even scalar sine-Gordon interactions. Extending its estimates to
  the inhomogeneous retarded cell vertices here requires work; the theorem
  cannot be quoted as an already established arbitrary-pattern result.

The source statements above are [read]. The algebra and scaling arguments
in Sections 1-6 are independently derived, with their scopes stated where
used. No literature claim of a pre-existing theorem for the full d=24
vacancy model is made.
