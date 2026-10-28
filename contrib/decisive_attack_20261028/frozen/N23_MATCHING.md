# Matching the actual frozen n=23 bond to a weak 23-channel action

## Result and scope

**[proved reduction, conditional only on the repository's finite-cell
massive-matching theorem]** Fix the scalar rails, weak loop, selector and
phases of C33's `n=23` certificate. For every prescribed positive tolerance
`eta`, and every fixed finite choice of output field weights, exponential
space/time moments and source-derivative order, there are **fixed positive**
conversion and port couplings `g, lambda` such that the selected frozen
sector of the literal 23-letter bond has an exact zero-temperature logical
action of the form

    S = S_scalar + S_star + V_static + V_retarded,
    ||V_static|| + ||V_retarded|| <= eta,                 (1)

after time is measured in its actual positive payload energy unit. The
norm in (1) is on nonconstant kernels. `S_scalar` retains the reference
energy density exactly and the scalar matching correction with its
separate local-marker bound; its reference part need not be small in
payload units. `S_star`
is the fixed certified 23-band quadratic reference. It has 46 simple Fermi
points and C33's positive quadratic/quartic momentum margins. All 46
branch-resolving sources have physical quasi-local representatives and
controlled matching corrections. The constants may be enormous but are
independent of the number of cells and can be chosen uniformly before
reducing `g` in the interval specified below.

This verifies the finite-cell and parameter hypotheses for the actual
`n=23` bond. It is not a proof that an arbitrary weak retarded action of
this type flows to the proposed 23-boson CFT. That logical RG and
source/CFT step remains the substantive infrared assumption.

The conditional dependency is specifically
[`FINITE_CELL_MATCHING.md`](../../direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md),
including its coefficient norms, auxiliary-vacuum limit and exact action
identity. This note does not replace its proof. The finite-order local
change of basis below is applied directly to the fixed frozen cells; no
mobile-wrapper hypothesis is imported.

## 1. Fixed parameters and the two couplings that remain free

Use the literal bond in
[`explicit_candidate/README.md`](../../explicit_candidate/README.md), with

    n=23, b=24, t=t_x=t_switch=1,
    D_x=D_M=1000, Delta_phys=4-2 cos(pi/25), delta=0,
    phi_x,0=pi/46, phi_x,+ = phi_x,- = 0,
    phi_port=pi/2,
    phi_d,0=pi, phi_d,+ = phi_d,- = pi/2.

Stub and switch phases vanish. Let

    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1),
    beta_rail=1/d_24,
    t_d=beta_rail 2^-136,
    alpha=sqrt(2/25) sin(pi/25),
    c_23=2 alpha^4 beta_rail^2 > 0.

Choose the already certified midpoint marker selector `mu_sel`, and `K=10`.
Only `g_+=g_-=g` and the port coupling `lambda` will be reduced. Set

    0<g<=g_star=2^-365,
    epsilon(g)=c_23 g^2 2^-214,
    physical data field=lambda^2 epsilon(g).              (2)

Thus the physical common field changes with these chosen microscopic
couplings, exactly as in C33. Once selected, every entry of the bond is
fixed before taking volume, temperature or infrared limits.

The defect-sector proof applies for all these choices and
`0<lambda<=2^-20`. Its positive excited-cell, selector and illegal-sector
bounds improve or remain valid when `g` and `lambda` decrease.

## 2. Exact finite cell, code and gap

Work first in one conserved legal marker sector with every cell of length
24. Its complete cell Hilbert space has dimension

    dim K=(7n+4)2^n=1,384,120,320.                          (3)

This is the full legal space, not only port-accessible excitations. There
are `n+1` vacuum-token positions and `3(2n+1)` excited-token positions, each
carrying an arbitrary n-bit data word.

Put the physical field to zero in the unperturbed cell Hamiltonian and
subtract its vacuum energy

    E_b=mu_sel-2 cos(pi/25).

Call the resulting positive operator `A_g`. Its exact kernel is

    G |a_1...a_n> = sum_(j=0)^n v_j
                     |vacuum token at j; a_1...a_n>,
    v_j=sqrt(2/25) sin(pi(j+1)/25).                        (4)

Neither the code isometry nor the code projection depends on `g` or
`lambda`. The vacuum walker is a length-24 open path. Its next level lies

    omega_v=2[cos(pi/25)-cos(2pi/25)]>0

above the code. In the excited-token sector the free scalar matrix relative
to `E_b` has diagonal 4 on the stub and 1004 on the express/marker rail. Its
unit-hop path adjacency has norm at most two, the extra loop edge has norm
`t_d`, and the conversion operator has norm `g`: at a fixed token position
its two nonzero conversion blocks are disjoint two-state blocks. Therefore

    A_g restricted to excited token >= 2-t_d-g >= 1.

The last inequality holds throughout (2). Consequently, uniformly for
`0<=g<=g_star`,

    ker A_g=ran G, dim ker A_g=2^23,
    A_g >= Delta_0 (1-GG^dagger),
    Delta_0=min(omega_v,1)>0.                              (5)

The marker potential cancels in this subtraction. The field in (2) is not
put into `A_g`, so no false assertion of an exactly degenerate finite-field
code is needed.

The cell charge is data occupation plus the token-species charge
`q_v=q_0=0`, `q_+=1`, `q_-=-1`. Every term of `A_g` preserves it. Its parity
is `(-1)^Q`, and on (4) charge and parity become `N_log` and `(-1)^N_log`.
The physical vacuum-code vector with empty data has charge zero. The
intercell port conserves total charge and is even.

### The global Jordan-Wigner dictionary, including the left-cell parity

This step needs more than an on-cell Fock label. Let `f_(j,i)` denote the
on-cell annihilator, with data sites ordered left to right, and put
`P_j=(-1)^(N_log,j)`. On the entire ordered chain define

    C_(j,i)=(product_(k<j) P_k) f_(j,i).

Then, in the ordinary tensor-product cell representation,

    C_(j,i)^dagger C_(j+1,l)
       = f_(j,i)^dagger P_j tensor f_(j+1,l).             (5a)

The parity on the left cell is essential. Without it the claimed
multichannel hopping would not be the displayed free-fermion Hamiltonian.

It is present in this physical construction. On each neutral-token rail
use the data gauge `U_i=product_(r<i) Z_r`, while leaving charged rails
unchanged. Use the same U_i on the stub and express vertex at data position
i, and `U_n=P_j` at the marker. A neutral hop carries the physical factor
Z_i and `U_(i+1)=U_i Z_i`; these cancel. Junction switches have the same
gauge on both endpoints. The weak loop from the last data site to the
marker has the same cancellation. Thus this gauge is valid on the actual
folded graph including its loop, not only on C22's original path.
Neutral-to-charged conversions become `g f_i` or `g f_i^dagger`.

The left neutral source is at i=0 and has gauge I; the right neutral
source at the marker has gauge P_j. The additional scalar rail-phase gauge
at C33's point makes that right source `-i P_j`; charged sources retain
only scalar phases. Therefore the terms with one conversion on each cell
contain precisely the factor in (5a), or its adjoint. Same-cell
second-order terms contain two conversions and are quadratic and even;
the zero-conversion terms have two identical source parities and are
scalar. This is the explicit parity mechanism behind (12), not an assumed
replacement of an ordinary tensor product by a graded one.

This agrees with C22 RESULT equations (10)-(12),(16), C33 PROOF sections
1,2,9, and the literal comparison code: `ycell.second_order` builds the
ordinary two-cell matrix, while `ycell.logical` fits against
`jw_ops(2*n)` in the *global* cell-then-data ordering. The finite tests
reported by C33 compare the source formula with that literal matrix.
The standalone check in this folder separately verifies the string
gauge and (5a) on occupation bases. Boundary closure retains the usual
fermion-parity-dependent flux on the ring. No mobile coordinate map is
used anywhere in this argument.

## 3. The actual perturbation and a local change of basis

The exact selected-sector Hamiltonian, after subtracting `M E_b`, is

    H(g,lambda)=sum_j A_g,j + lambda sum_j V_(j,j+1)
                         + lambda^2 sum_j F_g,j,
    ||V_(j,j+1)|| <= sqrt(3),
    F_g=epsilon(g) N_data,  ||F_g||<=n epsilon(g).         (6)

The junction switch acts identically on the preceding marker state and is
already an internal cell term. Thus no hidden intercell term has been
omitted from (6). Each port changes a vacuum token into an excited one,
so `PVP=0`. The data number operator preserves the vacuum-walker ground
state and never changes token species, hence `F_g` is reducing for the
local code and `P F_g P=epsilon(g) N_log`.

For an interval support X put `P_X=product_j P_j`, `Q_X=1-P_X` and
`R_X=Q_X (sum_j A_g,j)^(-1) Q_X`. Equations (5) give

    L_X(O)=R_X O P_X-P_X O R_X,
    [L_X(O),sum_j A_g,j]=-O^off,
    ||L_X(O)||<=2||O||/Delta_0.                            (7)

Here off denotes the part between `P_X` and `Q_X`. Set, support by support,

    S1=L(V),
    G2=F_g+[S1,V^diag]+[S1,V^off]/2,
    S2=L(G2), S=lambda S1+lambda^2 S2.

The fixed-polynomial local-SW calculation in C45 Section 3 uses only (7),
bounded one- and two-cell interactions, and `PVP=0`. These hypotheses have
now been verified directly in (3)-(6). Its proof therefore yields the
exact identity

    exp(S) H exp(-S)=A+B+R,
    B=lambda V^diag+lambda^2 G2^diag,
    P B P=lambda^2 sum_j [H2(g)_(j,j+1)+epsilon(g) N_log,j]. (8)

Each declared term of B is locally code-block diagonal. Cross-bond
second-order contractions vanish: an outer cell touched once by a port
remains excited and cannot return to its vacuum code under the unperturbed
resolvent. No global inverse with a presumed isolated extensive code band
is used.

For the interaction norm
`||O||_a=sup_j sum_(X containing j) exp(a|X|)||O_X||`, choose bounds

    v0=max(1,2 sqrt(3) exp(2a)),
    f0=n epsilon(g_star) exp(a),
    G0=f0+48 v0^2/Delta_0,
    C_SW=2^16(v0 G0/Delta_0+v0^3/Delta_0^2).               (9)

These bounds are uniform for `0<=g<=g_star`. The exact conjugation satisfies

    ||B||_a <= 2 lambda v0+2 lambda^2 G0,
    ||R||_a <= C_SW lambda^3,                             (10)

provided

    0<lambda<=min(v0/G0,Delta_0/(192v0)).                  (11)

The BCH series converges for this fixed S with volume-independent bounds;
this is not an assertion that the full SW perturbation series converges.
The reference also has a ground-relative machinery threshold by the
locally block-diagonal comparison in C45 Section 4. That threshold does
not separate its whole extensive excited spectrum.

## 4. Why reducing g preserves the certificate and controls interactions

C33's Section 9 uses the fixed scalar rail and loop above. Its even
conversion expansion gives

    H2(g)=a00 I+g^2[a02 I+Q2]+D4(g),
    ||D4(g)||<=alpha^4 g^4/4,       0<=g<=1/2.             (12)

`Q2` is a number-conserving quadratic two-cell operator. The constants and
Q2 do not depend on g. At zero conversions the source compression is a
scalar. At two conversions at most two logical ladder operators occur,
including the finite Jordan-Wigner strings required by the graded local
bond, so only a scalar and a quadratic term occur. C33 computes their loop
values directly; (12) bounds the *entire* remaining operator, including its
quadratic and scalar pieces.

Define the g-independent quadratic reference

    H_star=sum_j Q2_(j,j+1)/c_23 + 2^-214 N_log,
    W=c_23 g^2 lambda^2 >0.                              (13)

Every field-independent term, including scalar pieces of D4 and the
matched residual, belongs to S_scalar and is retained separately
throughout. The nonconstant coefficient norm is bounded by the full
coefficient sum used below. Scalar energy density is not confused with a
mode or discarded in a proposed pressure comparison.
The normalized Bloch symbol of H_star is exactly the leading symbol in
C33's integer certificate.

That certificate also bounds the exact finite-g quadratic projection.
Its normalized, conditioned error is

    <= 2n g^2 beta_rail^-2 tr(W0 W0^dagger),               (14)

and the exact integer test at `g_star=2^-365` bounds this by `2^-33`.
Everything multiplying g^2 is fixed and positive. Thus the same strict
root brackets, inertia changes, hopping margin and pair-sum margins hold
for **every** `0<g<=g_star`, not only the displayed endpoint. The certified
bounds are

    46 simple roots,
    delta_2>0.0140193,
    delta_4>3.1558e-7.

The leading symbol itself is included with zero error. Its root simplicity
also gives strictly nonzero physical derivatives. It must be used for
velocities before any singular congruence; a conditioned root metric is
not a velocity normalization. At this one fixed symbol all 46 speeds have
a positive minimum, however small it is. No uniformity in n is required.

For output field weight H>=2 and spatial weight a_out>=0, converting the
two-cell operator D4 to logical normal monomials costs at most
`(8H^2)^(2n)`. Including the two bonds incident on an anchor, define

    C_log=2 exp(2a_out) (8H^2)^46.

The normalized static correction therefore obeys

    ||V_static|| <= C_log alpha^4 g^4/(4 c_23 g^2)
                  = C_log g^2/(8 beta_rail^2).            (15)

This directly establishes smallness in a strong coefficient norm, not
merely in operator norm. In particular, for any eta>0 choosing

    0<g<=min(g_star, 2 beta_rail sqrt(eta/C_log))          (16)

makes this term at most eta/2. The reference H_star and its momentum
margins have not changed. This g choice is made before choosing lambda.
Shrinking lambda alone would not improve (15).

## 5. A finite auxiliary embedding with uniform constants

The full cell has equal even/odd dimensions. Subtracting its code gives

    N_even_exc=N_odd_exc=(7n+3)2^(n-1)=687,865,856.

Choose a joint energy/charge/parity eigenbasis of its excited space.
Represent an odd excited eigenvector by one dedicated auxiliary fermion,
and an even one by two dedicated auxiliary fermions, with their energies
summing to the original excitation energy. This requires

    m_aux=2 N_even_exc+N_odd_exc=2,063,597,568,
    F=n+m_aux=2,063,597,591.                               (17)

These astronomical numbers only bound an analytic representation. They
are not an enlargement of the microscopic physical dimension d=23, and
we do not propose constructing their Fock matrices numerically.

Every auxiliary energy is at least `g_mass=Delta_0/2`. Odd original
integer charge C is assigned to its one auxiliary. For even C, assign
charges 1 and C-1 to its two auxiliaries. Each is odd modulo two and their
sum is correct. The even embedding therefore intertwines the physical
U(1) and parity. The n logical modes carry charge 1. Jordan-Wigner strings
of even interval operators remain inside the interval; the usual two
parity-dependent ring boundary conditions are retained.

Off the auxiliary vacuum, the Gaussian extension uses the *global* logical
CAR generators, whose preceding-cell parity includes both logical and
auxiliary parity. In an ordinary cell tensor representation an adjacent
hopping is therefore `f_j^dagger Gamma_log,j Gamma_aux,j tensor f_(j+1)`,
not the code hopping tensored naively with auxiliary identity. Graded
regrouping on each support interval makes this extension unitarily
equivalent to the logical operator tensored with identity, preserving its
norm. It agrees with the code operator on the auxiliary vacuum. This
choice makes H_log auxiliary-free in logical-CAR symbols; the difference
from a naive extension is vacuum reducing and is absorbed in K_D. The
thermal addendum states this dictionary in detail, including ring parity.

Neither the dimensions in (17) nor any coefficient conversion bound below
depends on g. Eigenvectors may vary or become degenerate at g=0, but the
bounds depend only on their orthonormality, operator norms and the uniform
gap (5). No differentiable eigenbasis in g is needed: choose g first and
then any joint eigenbasis, using the same numerical bounds.

Choose a time weight `kappa_t=Delta_0/4`, set

    delta_det=2 sqrt(m_aux),
    H'=H+2 delta_det,
    alpha_det=2/(g_mass-kappa_t)=8/Delta_0,
    omega_det=2 alpha_det/delta_det^2=4/(m_aux Delta_0).

Choose the physical SW weight a in (9) strictly greater than

    a_out+F log(8(H')^2).                                (18)

The operator-to-coefficient conversion of the finite-cell theorem then
turns (10) into valid strong coefficient bounds

    k<=4 lambda v0+4 lambda^2 G0,
    r<=C_SW lambda^3,                                    (19)

for the vacuum-annihilating block interaction and the nonconstant leakage.
Here the factor two in k bounds the embedded block and its logical
extension. An extra fixed factor may be inserted for a different support
anchoring convention; doing so changes only the finite constants and
lambda thresholds. Scalar leakage is tracked by the theorem's local
marker bound.

For example, in addition to (11), impose

    lambda<=1/(128 omega_det v0),
    lambda<=(1/(16 omega_det C_SW))^(1/3).                (20)

Because (11) implies `lambda^2 G0<=lambda v0`, (19)-(20) give
`k+r<=1/(8 omega_det)`, leaving half the theorem's analyticity allowance
for sources. Also impose `C_SW lambda^3<=Delta_0/8`. Equations (10)-(11)
then make the full interaction strength smaller than `g_mass`, so the
embedding's invalid states have a positive threshold above the valid
sector. They are not an extra low-energy branch.

The finite-cell theorem now applies at this actual cell dimension. It
gives an exact zero-temperature retarded action correction bounded by
`4 C_SW lambda^3`, with the chosen exponential time/space norms, including
fixed additional moment losses when those are requested. Since a pinned
kernel norm of every field degree acquires the common factor `1/W` under
`tau=W t`,

    ||V_retarded|| <= 4 C_SW lambda/(c_23 g^2).            (21)

The temporal exponential rate becomes `kappa_t/W`; for `W<=1` this
preserves any requested weaker fixed time weight. Thus it suffices also
to choose

    lambda<=eta c_23 g^2/(8 C_SW),   W<=1,               (22)

to obtain the other eta/2 in (1). Every bound (11), (20), (22), the
invalid-sector bound, and `lambda<=2^-20` is a strict positive upper
bound once g is fixed. A fixed positive dyadic lambda can satisfy them
simultaneously. No lower bound on lambda competes with them.

## 6. Physical sources and their constants

Use the 46 leading-reference roots and the branch polynomials in
[`BRANCH_SOURCES.md`](../sources/BRANCH_SOURCES.md). They have R=45 and

    D_a=product_(b != a) sin^2[(k_a-k_b)/2]>0,
    ||O_a||<=n/D_a.

Each full density occupies at most R_src=91 cells, or 2184 physical sites.
Lift it with the product code isometry, extended by zero outside the code.
In the transformed frame this is a finite-range even charge-conserving
operator, reducing for the local all-code projection. In the original
physical frame use `exp(-S) O_phys exp(S)`. It is quasi-local by the same
finite-time local evolution used in (8).

For a completely explicit, deliberately loose coefficient bound on the
source's vacuum-annihilating embedded part one may take

    s_a=2 R_src exp(a_out R_src)
             (8(H')^2)^(F R_src) n/D_a.                  (23)

The factor two covers the embedded operator and its logical extension;
R_src bounds the number of translated supports through an anchored cell.
These constants are independent of g and lambda. Fix source radii

    rho_a=1/(8 omega_det 46 s_a).

Their total norm budget is exactly `1/(8 omega_det)`, so together with
(20) the source lemma has the required reserve. For every fixed
multi-index nu, its leakage-dependent action coefficient satisfies

    ||partial_z^nu F_match(0)||
       <=4 C_SW lambda^3 nu! product_a rho_a^(-nu_a).      (24)

Pure logical source vertices are retained exactly and are not included in
this small correction. Formula (24) is a bound on effective-action
coefficients, not yet on gapless correlation functions.

If sources are normalized to unit coupling in tau time, couple the
physical Hamiltonian to `W j_a(Wt) O_phys,a`. Equation (24), the chain
rule and the common `1/W` action rescaling give, with `p=|nu|`,

    ||partial_j^nu F_match'(0)||
       <=4 C_SW lambda^3 W^(p-1)
                         nu! product_a rho_a^(-nu_a).     (25)

For any prescribed finite source-derivative order these bounds can be
made arbitrarily small by taking the same fixed lambda still smaller.
In particular the first derivative is O(lambda^3), with finite constants.
This does not silently identify an undressed original observable with its
code projection.

## 7. Endpoint degree control is also finite

The instantaneous static correction acts on two cells with 23 logical
modes each, so its maximum normal-monomial degree is `P_static=92`.
Equation (15) can therefore be made small enough for any finite-degree
endpoint hierarchy by further reducing g.

For fixed g the exact residual action is analytic in lambda in the
finite-cell theorem's ball, starts at order lambda^3, and obeys a finite
order-degree bound. A coefficient at physical order q in the local-SW
expansion has support on at most q+1 cells. A product of vertices of total
order q has summed support at most 2q. After the embedding, its degree is
at most `4Fq`; auxiliary contractions cannot increase external degree.
Thus the residual satisfies the hypotheses of the elementary degree
reparameterization lemma with `D_degree=4F`. Analyticity supplies a bound
`||K^[q]||<=C A^q` by Cauchy's estimate on a fixed smaller lambda disk.
This recovers the higher-degree endpoint powers after taking lambda
sufficiently small. It does not establish the subsequent marginal beta
function or Ward-identity comparison.

## 8. Exact remaining obligation

The finite microscopic chain, complete finite cell, uniform code gap,
charge grading, positive bandwidth, strong coefficient smallness,
nonresonant reference and physical sources have now been specified in
one consistent parameter order:

1. Fix the certified rails, loop, phases and selector.
2. Fix the desired norm weights, source order and logical tolerance.
3. Fix a positive g satisfying (16); the certificate remains valid.
4. Compute the finite SW/embedding/source constants with these fixed data.
5. Fix a positive lambda satisfying all stated upper bounds.
6. Only then take volume, zero-temperature and infrared limits.

What remains is a logical theorem for this weak, translation-invariant,
U(1)-preserving, high-degree retarded 23-channel action: control the
marginal flow, retain the full rank of the source vertex matrix through
infrared renormalization, and identify 23 independent commuting chiral U(1) currents with a positive-
definite level matrix in a unitary conformal bulk limit (and the corresponding
left-moving sector). Such an Abelian current algebra would imply c>=23 by the existing Sugawara
argument and exceed d-1=22 for the frozen model's
bulk interpretation. Neither generic power counting nor a scalar pressure coefficient alone
proves it. Merely counting independent positive-level currents is also
insufficient if their algebra is non-Abelian: SU(2)_1 has three currents
but c=1. The Abelian/Gaussian operator-product structure is an explicit
remaining obligation, as emphasized in the parallel source review.

The separate [defect-pressure theorem](DEFECT_PRESSURE.md) removes all
other frozen sectors from any polynomial low-temperature coefficient,
conditional on that coefficient existing in the selected sector. It does
not erase the [incompatible-ring spectral obstruction](ALL_RING_OBSTRUCTION.md).
A finite-temperature version of the full source/scalar matching would
still be needed if pressure, rather than zero-temperature current algebra,
is used to establish the central charge.

Provenance: Codex GPT-6-family frozen-sector worker, 2026-10-28. The C33
proof and saved exact certificate were read directly. The coordinating agent independently checked the dimensions, gap, C33
error budget, SW constants, source radii and parameter schedule. The global
Jordan-Wigner dictionary is derived explicitly above and checked by the
standalone script. No human review or solved conjecture is claimed.

Reproduction:

    python3 -B contrib/decisive_attack_20261028/frozen/check_n23_matching.py

The standard-library script verifies 177,796 exact adjacent-cell
Jordan-Wigner identities and 3,586 neutral-wire gauge cases, including a
counterexample when the left-cell parity is omitted. It checks dimension
counts and saved exact C33 certificate fields. It does not rerun C33
inertia certification or perform an infrared calculation. Output is
[n23_matching_checks.json](n23_matching_checks.json).
