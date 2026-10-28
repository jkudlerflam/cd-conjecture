# Independent local perturbation audit and a dynamical response bound

## Claim and status

**[physics inference]** I find no concrete instability in the proposed
one-site invariant, nearest-neighbour, Hermitian d=24 chain at fixed,
sufficiently separated nonzero couplings. The finite-order local change
of basis used for its frozen 23-state component is legitimate. The
additional charge fluid has a bounded smooth-density response even at
nonzero frequency, which strengthens the stated reason that its weak
coupling does not overwhelm the internal scale.

This supports the claimed total Virasoro c=24 at the project's stated
theoretical-physics standard. It is not a new finite-size determination
of c, a rigorous microscopic thermodynamic theorem, or an independent
reconstruction of the 23-band certificate. The result depends on the
existing finite bond, band and material-sector certificates, and on the
usual compact-boson matching and weak-coupling RG assumptions.

## Scope of this independent check

The following were read: the complete d=24 synthesis, the global material
sector proof, the local SW application and source review, the physical
harmonic audit, the phase-model argument, and the tagged-defect response
analysis. I independently opened the primary SW paper rather than relying
on the earlier workers' descriptions. I did not repeat their matrix
inertia certificates or the original tagged-defect numerics.

**[read]** Bravyi, DiVincenzo and Loss allow a tensor product of degenerate
on-site ground spaces with a positive on-site gap. Section 4.4 constructs
a finite polynomial local anti-Hermitian generator. Equations (4.23)-(4.25)
split the exactly conjugated Hamiltonian into a block-diagonal truncation
and a remainder. Lemma 4.2 bounds coefficient interaction strengths, before
conversion to an extensive operator norm. Lemma 4.1 and Corollary 4.1
place a ground state of a sufficiently weak, locally block-diagonal
perturbation in the product code. These statements are visible in the
[primary paper](https://arxiv.org/html/1105.0675#S4.SS4), including the
[block-diagonal section](https://arxiv.org/html/1105.0675#S4.SS2).

A finite polynomial generator is essential here. No convergence as the
truncation order tends to infinity is being used.

## 1. The frozen-cell application passes

**[deduced from the cited local construction]** At fixed n=23, subtract
one vacuum-cell energy from each blocked cell. The cell space is finite,
the code projector P_j has fixed rank, and H_0,j >= Delta_0 Q_j with
Q_j=1-P_j and Delta_0>0. All coefficients, gaps, and dimensions are held
fixed before increasing the number of cells. The order-lambda-squared
data field can be included in the analytic local perturbation.

For any fixed truncation order m, write the exact transformed model as
D_m+R_m. Its retained restriction agrees with the computed order-two
Hamiltonian through that order. On a one-dimensional nearest-cell graph,
a coefficient of order r has a connected declared support of diameter
at most r. Consequently a factor exp(mu diameter) can be absorbed by
shrinking the permitted lambda. The local remainder has an exponentially
weighted strength bounded by C_(m,mu) |lambda|^(m+1), with constants
independent of the number of cells. This is a full-space interaction
statement, and does not assume isolation of the entire extensive code
band.

**[proved deduction]** Subtract a sum_j Q_j from D_m, with
0<a<Delta_0. The new on-site gap is at least Delta_0-a; sufficiently
small lambda still satisfies block-diagonal ground-sector stability.
The restriction to the global code is unchanged, with ground energy
E_eff. By the variational characterization of the lowest eigenvalue,

    D_m - a sum_j Q_j >= E_eff I,
    D_m >= E_eff I + a sum_j Q_j.

Commutation of D_m with each individual Q_j is unnecessary. It commutes
with the global code projector. On its complement sum_j Q_j>=1, so the
reference auxiliary threshold is at least a. This justifies using a
massive complement above the ground while allowing extensive code-band
overlap at high energies.

**[physics inference]** At fixed g the smallest retained speed is
lambda^2 w_*, with w_*>0. The normalized remainder at m=2 is bounded by
(C/w_*) |lambda|. It can be made smaller than any fixed Gaussian stability
budget. This uses stability of the certified phase under small local
perturbations, not a CFT theorem from the SW paper.

## 2. New check: the smooth charge response is bounded dynamically

The static cancellation in the existing defect analysis is correct,
but one can state a stronger continuum fact. It applies to quantum
internal sources as an influence-action kernel, as well as to static
scalar defects.

**[derived in the Gaussian charge theory]** Use material coordinate x,
charge sound velocity v, and the gap-number field phi with action

    S_gap = (1/(2 pi K)) integral dx d tau
            [ (partial_tau phi)^2/v + v(partial_x phi)^2 ],
    delta r = -partial_x phi/pi.

With the spatial zero mode omitted in the canonical sector, the Euclidean
density response is

    G_rr(k,omega) = (K v/pi) k^2/(omega^2+v^2 k^2),
    0 <= G_rr(k,omega) <= K/(pi v) = 1/e_gap''.

The upper bound is uniform in both momentum and frequency. In particular,
there is no denominator proportional to the inverse slow payload speed.
The noncommuting limits at k=omega=0 do not violate it; the canonical
zero mode is simply absent.

For A=alpha+a delta r+..., a=d alpha/dr, and coupling J A O to a
bounded local internal operator O, integrating this Gaussian part
formally gives

    delta S = -(J^2 a^2/2) integral_(k,omega)
                  G_rr(k,omega) O(-k,-omega) O(k,omega).

This expression is an influence action, not a claim that the charge
oscillator has disappeared. The charge mode must be retained when
counting the final central charge. The bounded kernel shows directly
that its smooth response has coefficient O(J^2/t), including dynamical
sources. After time is expressed in internal units of order J, its
relative coefficient is O(J/t), with the fixed internal norms and
slowest velocity included in the allowed tolerance.

For the exact half-filled phase-model reference,

    K=4,  v=t,  e_gap''=pi t/4,
    a=-(1/4+1/pi^2),
    a^2/(2 e_gap'') = 0.07857579749.../t.

**[derived canonical static specialization]** Setting omega=0 at nonzero
k gives a constant, so the finite-ring susceptibility is

    chi_forward(i,j) = C(delta_ij-1/N),
    C = a^2/(2 e_gap'').

The energy is -J^2 C sum_i s_i^2 + J^2 C (sum_i s_i)^2/N for static
sources s_i. A change of k bounded sources changes the rank-one term by
O(J^2 k), since k<=N. Its absolute first spatial moment need not stay
bounded. It therefore must remain separate from the decaying kernel,
as the current synthesis already does.

## 3. What this bound does and does not settle

**[standard continuum power counting]** The remaining local contact
operators have charge dimension at least two, apart from the Gaussian
current. Their connected two-point functions are integrable at long
spacetime separation. The dimension-two contribution has static tail
no slower than |x|^-3; the dimension-four first phase-model harmonic has
|x|^-7. These have a finite first spatial moment. The smooth current's
apparently dangerous absolute-integral tail is handled by its exact
bounded Fourier kernel above.

This is a useful local matching argument for the operator-valued charge
coupling. It does not alone prove a uniform microscopic linked-defect
estimate, since a continuum influence action is not an operator-norm
inequality for the lattice Hamiltonian. The established frozen-sector
comparison and the usual controlled continuum matching are still
required. A finite-size scalar-defect calculation also cannot replace
that assumption.

The corrections are small relative to J times a fixed positive material
defect cost when J/t is chosen last. The nonforward momentum thresholds,
fractional physical marker harmonics, and full unequal-velocity Gaussian
stability budget must remain in the argument. They cannot be replaced
by U(1) conservation alone. The synthesis retains all of these inputs.

**[audit verdict]** No false source attribution, illicit extensive-band
projection, missing marginal locking term, or contradictory scale order
was found in the reviewed argument. The absence of an exact all-order
block diagonalization is not a physical objection at the agreed
standard. The end-to-end assertion should continue to be described as a
controlled theoretical-physics construction with explicit continuum
inputs, not as a measured microscopic c=24 or a rigorous CFT-limit result.

## Provenance and reproduction

Independent GPT-6-family Codex worker, 2026-10-28. No human review is
claimed. The equations above can be checked by Gaussian integration and
ordinary differentiation; no new numerical evidence is represented as a
proof. Run `python3 contrib/vacancy_continuation/local_audit/checks.py`
from the repository root for the response-normalization and canonical
source-energy consistency checks. Only Python's standard library is
needed. The checks are diagnostic and do not establish the continuum
or microscopic hypotheses.
