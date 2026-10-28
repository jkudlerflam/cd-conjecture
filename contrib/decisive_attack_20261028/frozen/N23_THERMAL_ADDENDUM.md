# Thermal operator-valued matching for the frozen n=23 candidate

## Statement and dependency

**[proved specialization conditional on the cited abstract operator-polymer
theorems]** The actual fixed-cell n=23 model satisfies the gap, heat-trace,
coefficient-norm, grading, source and locality hypotheses of
[`INFINITE_CELL_PROGRESS.md`](../../direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md)
and its
[`TREE_NORM_EXTENSION.md`](../../direct_attack_20261028/infrared/TREE_NORM_EXTENSION.md).
Their conclusions therefore supply a beta-independent, weak retarded
logical endpoint, its thermal periodization with exponentially small error,
and bounded-source derivatives of the massive influence action, including
its scalar density. The output norm includes a terminal spanning-tree
weight, uniformly in external degree.

This is an application of the repository's newly proved abstract theorem,
not an assertion that a published literature theorem already states that
result. Its independent AI audit is
[`INFINITE_CELL_AUDIT.md`](../../direct_attack_20261028/alternatives/INFINITE_CELL_AUDIT.md).
No human review or central-charge conclusion is claimed.

The earlier [dedicated-fermion reduction](N23_MATCHING.md) remains an
independent route conditional on the separate determinant-based matching
theorem. This addendum uses a different auxiliary representation. It avoids
the two-billion-mode coefficient conversion and proves thermal statements
unavailable from that earlier vacuum-only output.

The full thermal influence action here is that of the enlarged comparison
representation. The original selected-cell pressure differs from its
pressure only exponentially at low temperature, by a separate exact
comparison below. This pressure bound alone is not promoted into a theorem
about derivatives of the full gapless pressure or physical correlators.

## 1. Exact microscopic inputs already verified

Use the actual parameter family, code, gap and physical Jordan-Wigner
identity of N23_MATCHING, with

    n=23, b=24,
    dim K=1,384,120,320,
    dim P=2^23=8,388,608,
    d_exc=dim QK=1,375,731,712,
    Delta_0=min(2[cos(pi/25)-cos(2pi/25)],1)>0.

At zero physical common field, `A_g` has exact kernel `G F_23` and
`A_g>=Delta_0 Q`, uniformly for `0<=g<=2^-365`. Its field is retained as
`lambda^2 epsilon(g) N_data` in the perturbation. The source and weak-loop
parity mechanism gives

    f_(j,i)^dagger P_j tensor f_(j+1,l)
       =C_(j,i)^dagger C_(j+1,l)

in the global cell-then-data fermion ordering. Thus the code reference is
actually quadratic in the global fermions. This point is not inferred
merely from having a Fock-space label on each cell.

Only the fixed legal length-24 cell space is used here. There is no
unbounded length source, movable cell boundary, physical-length constraint,
or marker-root quotient to verify inside this sector. On compatible rings,
there are M=L/24 ordinary tensor factors and the usual parity-dependent
fermionic closing bond. Competing microscopic sectors are handled by the
separate [defect-pressure theorem](DEFECT_PRESSURE.md).

## 2. Abstract auxiliary space and explicit heat trace

Take

    E=C Omega direct_sum QK,
    a_aux Omega=0,  a_aux|QK=A_g|QK,
    J G|u>=|u>_log Omega,
    J|s>=|0>_log |s>_aux,  s in QK.                      (1)

No excited eigenbasis or auxiliary Fock expansion is used. This is an even
isometry, and `(I_log tensor a_aux)J=JA_g`. Give the auxiliary states their
original integer charges and the logical modes charge one; the vacuum has
charge zero. Then J intertwines U(1) and parity. Even microscopic interval
operators remain local in the already verified global fermion dictionary.

There is exactly one auxiliary vacuum and a uniform excited gap Delta_0.
Since QK is finite-dimensional,

    Tr_Q exp(-t A_g)<=d_exc exp(-t Delta_0),
    Tr_E exp(-t a_aux)<=1+d_exc exp(-t Delta_0)             (2)

for every t>0, uniformly in the whole g interval. This is the required
heat-trace hypothesis, proved directly for the actual cell. No legal-length
counting estimate or unproved infinite-cell uniformity is imported.

The valid subspace is the code with auxiliary vacuum plus logical vacuum
with an excited auxiliary. Invalid states are the nonempty logical states
with excited auxiliaries; they have energy at least Delta_0. Every embedded
physical interaction vanishes if any cell in its declared support is
invalid. Validity projectors commute with the enlarged Hamiltonian.

### The Gaussian extension off the code includes auxiliary parity

Let Gamma_log,j and Gamma_aux,j be the two local gradings. With factors
ordered by physical cells and logical modes before auxiliaries within a
cell, a global logical annihilator is

    C_(j,i)=(product_(k<j) Gamma_log,k Gamma_aux,k) f_(j,i).

Consequently the Gaussian reference hopping off the code is

    C_(j,i)^dagger C_(j+1,l)
       =f_(j,i)^dagger Gamma_log,j Gamma_aux,j
                                              tensor f_(j+1,l).      (2a)

The extra Gamma_aux,j is essential. We choose H_log^CAR, the extension
formed from these global logical CAR generators. We do not choose the
naive ordinary-spin-tensor extension of the code hopping with identity
on every auxiliary factor.

On any support interval, the graded permutation which groups its logical
factors before its auxiliary factors maps O_log tensor I_aux to the
corresponding CAR extension. This is unitary and preserves the operator
norm. On the auxiliary vacuum all Gamma_aux factors equal one, so the
code compression is exactly the already verified physical code operator.
The difference from the naive extension vanishes on the auxiliary vacuum
from both sides and belongs to K_D if that naive extension was used in
intermediate notation. Thus the norm bound by twice the original block
norm and the vacuum-reducing property are unchanged.

In logical-CAR symbols H_log^CAR has no auxiliary coefficient other than
identity, and therefore factors out of the auxiliary influence functional.
Odd auxiliary coefficients use the graded tensor product; outside-support
Jordan-Wigner strings in an ordinary matrix realization are not discarded.
Their combined even vertices satisfy the locality/factorization convention
of the operator-polymer theorem. The same choice is used for logical
source extensions. Ring closure retains total fermion parity and its
associated boundary flux. No independent choice of a logical-only ring
parity is substituted off the code.

## 3. Fixing the tree norm before selecting couplings

Choose a desired logical field weight H>=2 and spatial metric exponent
`a_out>=0`. Choose, as prescribed in the tree theorem,

    L_weight>=8(a_out+log(256)+4),
    tau>=max(2 L_weight/Delta_0,
             log[d_exc/(e-1)]/Delta_0),
    kappa_tree=L_weight/(32 tau),
    v_*=1/(8 tau).                                      (3)

Equations (2)-(3) give `Tr_E exp(-tau a_aux)<=e` and all the thermal
hypotheses for beta>=2tau. A larger tau is harmless if additional reserves
are desired. These constants are independent of g, lambda, volume and beta.

For external spacetime points X, the output weight is

    exp[T_(a_out,kappa_tree)(X)],

where T is the minimum terminal spanning-tree length in the anisotropic
space/time metric. At finite beta use circular time distance. The theorem
controls this weight directly from occupied spacetime networks; no
field-degree multiplier on a diameter bound is inserted.

Expanding an operator on r cells only in its 23r logical modes gives

    sum_I H^(degree I)||O_I||
       <=(8H^2)^(23r)||O||,                              (4)

with bounded auxiliary-operator coefficients. Thus choose the local-SW
support weight

    a_SW>L_weight+23 log(8H^2).                          (5)

This replaces the dedicated-fermion condition involving F=2,063,597,591.
The auxiliary spectral multiplicity enters through the finite heat time
in (3), not through a monomial expansion of its vacuum projector.

## 4. Local reduction and an explicit positive-coupling schedule

Recompute the fixed local-SW constants at the weight (5):

    v0=max(1,2 sqrt(3) exp(2a_SW)),
    f0=23 epsilon(g_star) exp(a_SW),
    G0=f0+48v0^2/Delta_0,
    C_SW=2^16(v0G0/Delta_0+v0^3/Delta_0^2).

They are uniform for `0<=g<=g_star`. The exact conjugation and its logical
compression are the same formulas (8)-(12) in N23_MATCHING, although its
finite local unitary may use a different admissible lambda because the
norm weight has changed. The strong operator-valued coefficient norms obey

    k<=4lambda v0+4lambda^2G0,
    r<=C_SW lambda^3.                                   (6)

Here K_D is the embedded locally block-diagonal reference minus its
logical CAR extension specified in (2a); every one of its auxiliary coefficients annihilates Omega from
both sides. R includes its scalar coefficient in the strong density norm,
as required by the operator theorem.

Keep the microscopic selector and SW restrictions

    lambda<=2^-20,
    lambda<=min(v0/G0,Delta_0/(192v0)).

Also impose

    lambda<=v_*/(128v0),
    lambda<=(v_*/(16C_SW))^(1/3),
    C_SW lambda^3<=Delta_0/8.                            (7)

Since `lambda^2G0<=lambda v0`, equations (6)-(7) give `k+r<=v_*/8`, leaving
an equal source reserve. The ordinary interaction strength of B+R is then
less than Delta_0/2, which also verifies the invalid-sector comparison's
hypothesis. All right-hand sides are strictly positive and fixed.

Let B_tree be a valid finite output constant of the cited tree theorem,
including its harmless periodization enlargement, and put

    C_match=4 B_tree/v_*.

Its dependencies are the already fixed gap, heat time, norm weights and
support conventions; it does not depend on g, lambda, M or beta. The
thermal influence map then obeys

    ||W_beta[K_D+R]-W_beta[K_D]||<=C_match C_SW lambda^3,
    ||W_beta[K_D]||<=B_tree exp[-beta/(2tau)].             (8)

The first norm includes the separate scalar-density convention. The
second includes the free auxiliary scalar if that scalar has not been
normalized out. Its zero-temperature limit is zero. The logical action
factors out exactly, without introducing an inverse logical bandwidth
in the massive integration.

Choose g first using the static error bound in N23_MATCHING. With
`c_23=2 alpha^4 beta_rail^2`, the physical bandwidth is

    W=c_23 g^2 lambda^2>0.

To make the normalized massive endpoint at most eta/2, it suffices to add

    lambda<=eta c_23 g^2/(2 C_match C_SW).                (9)

Reducing g preserves C33's certificate, and reducing lambda after it
satisfies every upper bound. No microscopic coupling depends on temperature
or ring length.

## 5. A single endpoint at every low temperature

Let E_infinity be the vacuum influence action for K_D+R, including its
scalar density. The baseline cancellation and (8) give

    ||E_infinity||<=C_match C_SW lambda^3.

The thermal periodization theorem gives

    ||W_beta[K_D+R]-Per_beta E_infinity||
       <=C_per B_tree exp[-beta/(2tau)],                 (10)

where C_per is another fixed finite constant, or can be absorbed into
B_tree. The periodization has the antiperiodic signs of odd logical
fields and circular tree metric. This is a common beta-independent initial
endpoint, not a family of uncontrolled beta-dependent interactions.

With normalized time `u=W t` and inverse temperature `beta_log=W beta`,
every pinned kernel norm, irrespective of degree, gains a factor 1/W.
Consequently

    ||E_infinity'||<=C_match C_SW lambda/(c_23 g^2),
    ||E_beta'-Per_(beta_log) E_infinity'||
       <=(C_per B_tree/W)
                      exp[-beta_log/(2tau W)].           (11)

Scalar coefficients are divided by W as energy densities. The original
reference scalar, including the order-lambda^2 term which is large in
bandwidth units, is retained separately and exactly.

For fixed positive g and lambda, (11) is exponentially small as temperature
tends to zero. At the mere physical threshold beta=2tau it need not be
small after division by an extremely narrow W. For any desired normalized
error zeta>0, the sufficient additional temperature condition is

    beta>=max(2tau,2tau log[C_per B_tree/(zeta W)]).

This is a restriction on the range in which the weak logical thermal
endpoint is used, not a temperature-dependent choice of Hamiltonian.

## 6. Bounded branch sources

The 46 branch probes already have finite physical lifts. For the full
91-cell density, with `D_a` as in BRANCH_SOURCES, a conservative strong
coefficient bound in the present representation is

    s_a=2*91 exp(91 L_weight) (8H^2)^(23*91) 23/D_a.       (12)

It includes the embedded lift and its logical extension. Set
`rho_a=v_*/(8*46*s_a)`. The total source allowance is v_*/8; together with
(7) it remains in the theorem's v_*/4 ball.

Each source auxiliary remainder is vacuum reducing. The same thermal
Cauchy argument therefore bounds, for every fixed multi-index nu,

    ||partial_z^nu [W_beta(K_D+D(z)+R)-W_beta(K_D+D(z))]_(z=0)||
      <=C_match C_SW lambda^3 nu! product_a rho_a^(-nu_a).

The baseline and the thermal-to-periodized-vacuum differences have the
same derivative factor multiplying their exponential bounds. These are
bounds on massive effective-action kernels, including scalar density,
with all logical fields still external. In normalized sources coupled
physically as `W j_a(Wt) O_a`, the chain rule adds `W^(|nu|-1)` after
normalizing the action by W.

The original-frame observables are dressed by the local-SW unitary and
are quasi-local with controlled tails. No full interacting logical
correlation function, current level or pressure derivative follows
merely from these coefficient bounds.

## 7. Exact original/enlarged pressure comparison

For the finite cell (1),

    z_invalid(beta)=(2^23-1) Tr_Q exp(-beta A_g),
    z_valid(beta)=2^23+Tr_Q exp(-beta A_g).

Let `w=sup_j sum_(X containing j)||(B+R)_X||<Delta_0/2`. In each fixed
invalid-site set all touching embedded interactions vanish. Cutting and
restoring them costs at most w per invalid site, so the exact trace
comparison gives, per selected cell,

    0<=f_valid(T)-f_enlarged(T)
       <=T log[1+((2^23-1)/2^23) d_exc
                         exp[-(Delta_0-w)/T]].           (13)

The local-SW unitary preserves the original finite-volume trace exactly.
Restore the subtracted E_b and divide by b=24 for physical energy density.
The separate frozen defect-pressure theorem then compares this selected
pressure to the complete 23-letter microscopic pressure, with another
exponentially small correction. Thus all these pressures have the same
polynomial low-temperature coefficients, **if those coefficients are
established for the final interacting logical theory**.

Equation (13) does not imply uniform bounds on derivatives of the full
logical pressure or on arbitrary physical thermal correlators. A gapless
partition function can have complex-source zeros approaching the real
axis as beta grows. The Cauchy bounds above apply to the massive influence
action, whose logical fields are external, and not to that partition
function. Transferring full physical source correlators across invalid
sectors would require its own marked comparison. The zero-temperature
valid ground sector and code-reducing observables do not have this
thermal-sector issue.

## What is now isolated

Subject to the cited operator-polymer theorem, the frozen n=23 example
now has explicit microscopic hypotheses for full thermal massive
matching, a degree-uniform tree norm, one beta-independent retarded
endpoint, physical branch probes, and an exponential transfer of pressure
through all auxiliary and frozen defect sectors. The remaining analytic
work is the logical retarded multichannel RG, its normalized source
comparison, and a sufficient thermal/current-algebra identification of
central charge. Those conclusions are not supplied by the matching
estimate or by smallness alone.

Provenance: Codex GPT-6-family frozen-sector worker, 2026-10-28. This
specialization was derived from the two linked abstract theorems and the
already checked literal n23 cell. It retains the dedicated-fermion
argument independently and is awaiting its own cross-review.
