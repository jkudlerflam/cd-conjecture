# A local dressing for the full atomic cell space

Status: proved one-step normal form and fixed-slow-time comparison.
These estimates are uniform in volume and in the upper on-cell energy.
They do not establish an all-order infrared reduction or a CFT.

## 1. Hypotheses and norms

Consider a finite open chain or ring of M cells. Each cell may have a
separable, infinite-dimensional Hilbert space. Let A_j be positive
self-adjoint on that space, P_j its finite-dimensional kernel, and

    A_j >= Delta (1-P_j),    Delta>0.

The operators A_j act on separate cells. Write A=sum_j A_j and
P=product_j P_j. The perturbation is V=sum_X v_X, where X is a site or
neighboring pair, v_X is bounded and Hermitian, and

    J=sup_j sum_(X contains j) ||v_X|| < infinity.

A scalar energy per cell can be omitted. All statements concern the
specified interaction decomposition. For a general represented
interaction Phi define

    ||Phi||_kappa=sup_j sum_(X contains j) exp(kappa |X|)||Phi_X||.

The unweighted norm is ||Phi||_0. Supports may be connected intervals
on a ring. Repeated supports are allowed and their norms are summed.
Counting representation terms rather than combining them can only make
the upper bounds larger.

## 2. Explicit generator and exact first cancellation

[proved] For every X let

    P_X=product_(j in X) P_j, Q_X=1-P_X, A_X=sum_(j in X) A_j,
    R_X=Q_X (A_X|Q_X)^(-1) Q_X,
    O_X=Q_X v_X P_X+P_X v_X Q_X,
    D_X=P_X v_X P_X+Q_X v_X Q_X,
    S_X=R_X v_X P_X-P_X v_X R_X.

Then S_X is bounded and anti-Hermitian, and

    ||R_X||<=1/Delta, ||S_X||<=||v_X||/Delta,
    ||O_X||<=||v_X||, ||D_X||<=||v_X||,
    [S_X,A]=-O_X.

Proof. In the decomposition P_X plus Q_X, S_X has off-diagonal
blocks B and -B dagger, where B=R_X v_X P_X. Its norm is ||B||.
The same block argument proves the O_X estimate. D_X is the
block-diagonal conditional expectation of v_X and is contractive.
Finally A_X R_X=Q_X and A_X P_X=0 give
S_X A_X=-P_X v_X Q_X and A_X S_X=Q_X v_X P_X.
All A_j outside X commute with S_X.

These identities are valid on Dom(A). The resolvent maps into
Dom(A_X); S_X preserves Dom(A), and its commutator with A has a
bounded extension. On each finite ring, the exponential of the finite
sum S=sum_X S_X preserves Dom(A) and the commutator integral formula
holds. This can also be proved first with bounded spectral cutoffs of
A_j, then using the uniform commutator identities to pass to the
strong resolvent limit. No norm of A_j enters.

Set U=exp(epsilon S), D=sum_X D_X and O=sum_X O_X. Then

    U(A+epsilon V)U dagger = A+epsilon D+epsilon^2 R(epsilon),    (1)

where the exact remainder is

    R(epsilon)=int_0^1 exp(t epsilon ad_S)[S,V] dt
               -int_0^1 (1-t)exp(t epsilon ad_S)[S,O] dt.       (2)

Here ad_S B=[S,B]. Expanding only the first commutator of A proves
(2). The remaining expansion contains bounded interactions only.

[proved] Every D_X preserves the global all-code space P, so [D,P]=0
and PDP=PVP. Thus (1) has the exact first-order target in its P block.
It is not necessary to divide by logical energy differences.

## 3. A volume-independent remainder bound

[proved] For kappa>=0 and 2 exp(kappa)|epsilon|J/Delta<1, (2) has an
interaction representation satisfying

    ||R(epsilon)||_kappa
       <= 12 exp(3 kappa) J^2/Delta
          /[1-2 exp(kappa)|epsilon|J/Delta]^5.                 (3)

Proof with the combinatorics included. Let a represented interaction
Phi have support size at most r and norm J_Phi=||Phi||_0. An intersecting
pair S_X,Phi_Y has support X union Y, size at most r+1. Split the
anchored sum according to whether the anchor is in X or in Y. Using
||[B,C]||<=2||B||||C|| gives

    ||[S,Phi]||_0 <= 2(r+2)(J/Delta) J_Phi.                    (4)

The two contributions use sum_(Y meets X)||Phi_Y||<=2J_Phi and
sum_(X meets Y)||S_X||<=rJ/Delta. Double counting is harmless.
Iterating (4), starting with support size r, yields

    ||ad_S^n Phi||_0
      <= (2J/Delta)^n (r+n+1)!/(r+1)! J_Phi.

Its support size is at most r+n. Consequently the exponential series
has weighted norm at most

    exp(kappa r) J_Phi
      sum_(n>=0) binomial(n+r+1,r+1)
                   [2 exp(kappa)|epsilon|J/Delta]^n
    =exp(kappa r)J_Phi/[1-2 exp(kappa)|epsilon|J/Delta]^(r+2).

Both [S,V] and [S,O] have support size at most three and unweighted
norm at most 8J^2/Delta by (4). The integrals in (2) have total weights
one and one half. Replacing t by one in the bound gives (3).
This also proves absolute convergence of this single conjugation in
the represented interaction norm. It does not prove convergence of
successive block diagonalizations.

For example, |epsilon|J/Delta<=1/8 and kappa=log 2 give the conservative
constant ||R||_(log 2)<=3072 J^2/Delta. The unweighted bound is much
smaller. Neither bound depends on M or the dimension of a cell.

## 4. Charges, physical length, and the C8 application

[proved] Suppose A_X and v_X commute with the total charge on X.
Their kernel projections, reduced inverses and S_X do also. Therefore
U preserves that charge. This applies to total microscopic U(1), total
cell length, and any additional exact charge. Translation covariance
of the local terms makes the sum S covariant as well.

In unprojected_cell_reduction.md, take the complete direct-sum cell
space K, the stated atomic gap Delta_at, and the full bounded endpoint,
field and marker-transfer terms. For the flat wrapper one may use

    J <= |kappa_marker|+B_field+2 B_cut+4u <= 2 B_flat.

Here kappa_marker is the marker chemical potential in the perturbation,
not the norm weight in (3). The first inequality counts each cut term
at both of its incident cells. All quantities are the full-space norms,
not their projected visibilities.

The generator therefore exists on the infinite direct sum without
truncating long cells. It preserves sum ellhat=L. The cyclic sum of its
terms is compatible with the rooted-ring reference action, so (1)
descends to the physical quotient. Equivalently, work at each physical
momentum with its induced seam twist; unitary twists leave the estimates
unchanged. This uses the exact cell representation from run 9.

Illegal-selector sectors are invariant and are not included in this
legal-sector theorem. Comparing their ground energy with the dressed
legal sector remains a separate task.

## 5. The explicit second-order Hamiltonian

[proved, finite order] At epsilon=0 the remainder is

    R(0)=[S,V]-(1/2)[S,O]=[S,D]+(1/2)[S,O].

Its global P compression is the usual second-order expression

    P R(0) P = -P V Q (A|Q)^(-1) Q V P, Q=1-P.                (5)

To see this without a band-gap assumption, note that SP=Q A_Q^(-1)VP,
PS=-PV A_Q^(-1)Q and OP=QVP. Also [D,P]=0, PSP=0, hence
P[S,D]P=0. The two terms of P[S,O]P agree and have the negative sign
in (5). The inverse here is that of the unperturbed A at finite volume;
it is not the dressed bare-Q inverse ruled out in run9_audit.md.
Nonzero terms in (5) involve overlapping insertion supports, so their
cell support has size at most three. This improves the coordinate
bookkeeping, while leaving the finite-order scope unchanged.

## 6. A uniform fixed-time statement in slow time

[proved] Let O_Z be a bounded operator on a fixed set Z of cells. In
the dressed variables, compare the full dynamics generated by the
left-hand side of (1) with A+epsilon D. For t real and
J_R=||R(epsilon)||_0, one has the safe bound

    ||tau_full(t)(O_Z)-tau_bd(t)(O_Z)||
       <=6 |Z| ||O_Z|| epsilon^2 J_R |t|
                                      exp(8|epsilon|J|t|).     (6)

Proof. Pass to the interaction picture for the sum of on-cell A_j.
This changes no support or interaction norm, even for unbounded A_j.
For the nearest-neighbor interaction epsilon D, define the propagation
matrix K_xy=sum_(X contains x,y)||D_X||. With mu=log 2 its weighted
row sum is at most 2 exp(mu)J=4J. Iterating the commutator integral
inequality gives the Lieb-Robinson estimate

    ||[tau_bd(t)(O_Z),B_Y]||
       <=2 ||O_Z||||B_Y|| exp(8|epsilon|J|t|)
                            sum_(z in Z,y in Y) exp(-mu dist(z,y)).

The bound is deliberately also nonzero for initially disjoint supports
at t=0. It remains valid for time-dependent on-cell conjugates.
Sum this bound over the terms R_Y. At each y their norm sum is at
most J_R. On a chain or ring, sum_y exp(-mu dist(z,y))<=3.
The Duhamel comparison of the two evolutions has the single factor
epsilon^2 and a time integral. Bounding the exponential by its endpoint
proves (6). Absolute interaction convergence in (3) justifies the sum.

Compressing the comparison by P gives reference dynamics with
Hamiltonian epsilon PVP and observable PO_ZP, because [A+epsilon D,P]=0.
This applies to every initial state in P, including entangled states.
The physical preparation and observable are dressed by U dagger.
For a bare local observable the additional dressing error obeys

    ||U O_Z U dagger-O_Z||<=2|epsilon||Z|J||O_Z||/Delta,          (7)

by integrating U_s[S,O_Z]U_s dagger and using the local commutator
bound. The dressed observable is quasi-local by the same estimates.

For epsilon>0 and T=epsilon t fixed, (6) is O(epsilon), uniformly
in volume. It is a rigorous slow-dynamics limit. The exponential in
T prevents using it at arbitrarily late times at one fixed epsilon.
It does not compare ground states or establish conformal towers.

## 7. Where this line stops

The whole direct sum over cell lengths no longer obstructs a one-step
local dressing or finite-slow-time dynamics. Only Delta and bounded
interaction norms are used. Iterating the construction is possible to
any prescribed finite order, but (3) does not give a uniform convergence
bound for that iteration: supports and the inverse-problem bookkeeping
grow. No such convergence is asserted here.

The infrared question takes L and times to infinity at one fixed
nonzero epsilon. Neither an O(epsilon^2) per-cell remainder nor (6)
controls that limit. The soluble massive-auxiliary test in
critical_crossover.md shows why this qualification matters, even when
the first-order model is a familiar critical chain.

## 8. Finite checks and numerical conditioning

[numerical] Run local_dressing_checks.py with single-thread BLAS.
The tests use two to four four-state cells, dimensions 16,64,256,
with a two-state local kernel and charge-conserving random bonds.
The upper local excitation energy is varied through 1,17,10^6 while
the gap remains one. The generator identities, conserved charge,
complete second-order compression, quadratic remainder scaling and
linear slow-time error all pass. Ground-energy and atomic-occupation
bounds from defect_energy_bounds.md are tested after dressing as well.

The large-energy fixture exposed an eigensolver precision issue.
Direct lowest-eigenvalue extraction lost a positive energy difference
of order 10^-11 against an upper energy of 10^6, returning a spurious
negative difference of order 10^-10. The test now obtains the ground
vector from a shifted positive inverse and evaluates its Rayleigh
quotient. At epsilon=.01 the difference is 1.12177e-11 and the full
residual is 4.98e-10; the residual is in the large-energy components.
The ordinary-energy fixtures have residuals near machine precision.
This repair concerns the numerical test, not the analytic bounds.
No original microscopic C8 CFT was simulated.
