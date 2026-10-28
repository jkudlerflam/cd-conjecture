# Two charge gradings, separate Jordan-Wigner maps, and auxiliary matching

Status: **[proved algebraic locality and norm statements]**. This removes
the common-Jordan-Wigner qualification in the proposed massive matching
application. It proves no convergence or infrared estimate.

Provenance: OpenAI GPT-6-family Codex agent `exact_escape`, 2026-10-28.
The coordinator identified the required Gaussian extension and requested
its norm audit. No human review. The argument below is finite-volume
algebra, with bounds independent of volume and auxiliary dimension.

## 1. The issue and the required two symmetries

The retained cell is F_n,payload tensor C^2,length. Its reference is
n payload fermion bands plus an independent hard-core length chain.
One interleaved Jordan-Wigner ordering of all physical cell factors
can turn one of these two hopping terms into a correlated parity-string
term. It is incorrect to call both terms quadratic after that map.

There are two independent exact integer charges: payload charge and
total physical length minus b times the marker count. On the code they
are respectively N_payload and N_length, with N_length=0,1. Work with
interactions neutral under both charges separately. Neutrality implies
evenness under each charge parity, not merely their product.

Individual interaction terms may be chosen neutral without increasing
their norm or support: average each term over the two compact U(1)
actions. If the full interaction is invariant, the averaged sum is the
same Hamiltonian. The actual local rules and their local SW reductions
already have this termwise property.

## 2. Include auxiliary charges before applying Jordan-Wigner

Let a physical cell K have finite code F_n tensor C^2 and arbitrary
separable complement QK. Set E=C Omega direct_sum QK. The embedding is

    J|u,eta>_code=|u,eta>_log tensor Omega,
    J|s>_Q=|0,0>_log tensor |s>_aux.                       (1)

The auxiliary vacuum has both charges zero. On QK assign the auxiliary
payload charge to be the original cell charge, and auxiliary length
charge to be ell-b. Negative auxiliary length charges are allowed and
cause no parity problem. Then

    Q_payload,total=N_payload,log+Q_payload,aux,
    L_cell=b+N_length,log+Q_length,aux                      (2)

intertwine the original physical charges on ran J. If the original code
has a constant payload offset q0 per cell, subtract q0 throughout (2)
and restore its marker chemical term separately.

The auxiliary space has commuting parity involutions

    pi=(-1)^Q_payload,aux,    lambda=(-1)^Q_length,aux.

The on-cell auxiliary Hamiltonian commutes with these two charges and
has vacuum Omega. Its dimension can be infinite. A charge-exchange
vertex can be odd in its *logical* payload or length fields, compensated
by odd auxiliary charge. Applying Jordan-Wigner only to logical factors
would therefore leave strings outside its support. The construction
below uses the two total cell parities instead.

## 3. An explicit ordinary CAR and graded auxiliary representation

Work on a finite open interval of cells j=1,...,M. Initially all distinct
cell factors are ordinary tensor factors. Let c_(j,a), a=1,...,n, be
local payload CAR matrices, and b_j the hard-core lowering matrix on
the length qubit. Different cell operators and the two logical families
initially commute. Define the diagonal involutions

    p_j=(-1)^N_payload,log,j,   q_j=(-1)^N_length,log,j,
    P_j=p_j pi_j,              Q_j=q_j lambda_j,
    S_j=product_(k<j) P_k,     T_j=product_(k<j) Q_k,
    K=product_(k=1)^M P_k.                                  (3)

Use separate strings and one global Klein factor:

    psi_(j,a)=S_j c_(j,a),
    chi_j=K T_j b_j.                                       (4)

These are n+1 ordinary CAR modes per cell. Within each family this is
the standard string calculation: a lowering operator anticommutes
with its own total parity and commutes with the other cells' parities.
Between families, all factors commute except that psi anticommutes
with K. Hence psi and chi anticommute, including their adjoints.

Here is the auxiliary representation, needed to make (4) useful for
operator-valued Grassmann kernels. Decompose a bounded auxiliary
operator into its four homogeneous components O_(r,s), r,s in {0,1},
where

    pi O_(r,s) pi=(-1)^r O_(r,s),
    lambda O_(r,s) lambda=(-1)^s O_(r,s).

Define

    Xi_j(O_(r,s))
      =(K pi_j)^s S_j^r T_j^s p_j^r q_j^s O_(r,s),         (5)

and extend linearly. Since K pi_j=p_j product_(k!=j)P_k,
every factor before O in (5) commutes with the entire original
auxiliary algebra at j. The factors multiply as the characters of
Z2 x Z2. Thus Xi_j is a unital *-homomorphism on B(E_j).

It is isometric. Indeed put U_j=S_j p_j and
V_j=(K pi_j)T_j q_j. They are commuting involutions in the commutant
of B(E_j). In a simultaneous U_j,V_j eigenspace, (5) is

    O -> sum_(r,s) u^r v^s O_(r,s),  u,v in {+1,-1}.

This is conjugation of O by a product of pi_j and lambda_j. Every
nonempty block has exactly norm ||O||, proving ||Xi_j(O)||=||O||.
This reasoning works for infinite E and does not expand it into
auxiliary creation operators.

Xi_j(O_(r,s)) has ordinary total parity r+s modulo two. Directly
commuting the factors in (3)-(5) gives

    psi_i Xi_j(O)=(-1)^(r+s) Xi_j(O) psi_i,
    chi_i Xi_j(O)=(-1)^(r+s) Xi_j(O) chi_i,                (6)

including i=j. At different auxiliary sites,

    Xi_i(O_(r,s)) Xi_j(V_(u,v))
      =(-1)^((r+s)(u+v)) Xi_j(V_(u,v)) Xi_i(O_(r,s)).      (7)

For example, the exponent in (7) is the sum
ru+sv+rv+su: the first two terms come from the separate strings and
the last two from the Klein factors. This is the usual single-parity
graded tensor-product rule. Within one auxiliary site Xi preserves
the full original multiplication, which need not be a free-fermion
algebra. Thus the representation is ordinary CAR for the logical
fields and an arbitrary graded operator algebra for the auxiliaries.

The auxiliary Hamiltonian and both auxiliary charge/parity operators
are even under both gradings, so Xi leaves them unchanged. In particular
the product auxiliary vacuum and its propagator are the same ones.
There is no new auxiliary spectral gap or denominator.

The physical-length fugacity is compatible without a new estimate:
logical number operators are unchanged by (4), and auxiliary charges are
unchanged by (5). Thus a tilt s L adds s N_length,log to the logical
Hamiltonian and s Q_length,aux to the auxiliary Hamiltonian, plus sb per
cell. The already derived auxiliary gap and heat-trace window applies.

## 4. Neutral interval operators remain local and retain their norm

Let X=[a,b] be an interval and O_X a bounded operator supported there,
neutral under both total charges. Expand it into homogeneous products
of local logical matrices and auxiliary operators, and use the inverse
of (4)-(5). Strings strictly left of a occur with exponents equal to
the total payload parity and total length parity of the product. Both
are zero. The global factor K occurs with the total length parity,
also zero. All remaining parities have sites in X. Reordering produces
signs but no operator outside X.

Consequently O_X is an operator in the CAR/graded-auxiliary algebra on
the same interval X. This argument includes vertices with logically odd
fields and compensating odd auxiliary operators; only the *total*
charges of the full vertex are required to be neutral. Odd individual
fields need not be local, and are not asserted to be local.

Equations (3)-(5) are exact operator identities in one faithful Hilbert
representation. Therefore the rewritten O_X has exactly its original
operator norm. No coefficient expansion is used for this isometry.
The result passes to the infinite-chain neutral local algebra because
the formula for an interval is independent of the number of exterior
sites after cancellation. Quasi-local interactions follow term by term
with the same support and interaction-norm bounds.

At a finite periodic seam, the separate strings create the usual
charge-parity boundary twists. The physical pressure argument already
cuts that seam at a bounded cost before taking the thermodynamic limit.
This note does not replace the finite-ring gluing rule by an incorrect
unconstrained periodic fermion trace.

## 5. The reference extension must be chosen in CAR variables

The correct quadratic extension is

    H_G=sum_j [psi_j^dagger t psi_(j+1)+h.c.
                  -t_len(chi_j^dagger chi_(j+1)+h.c.)]
          + charge-preserving on-cell quadratic terms.    (8)

It is independent of the auxiliary algebra in this CAR representation.
In the original tensor variables, the same terms are

    psi_(j,a)^dagger psi_(j+1,b)
       =c_(j,a)^dagger p_j pi_j c_(j+1,b),
    chi_j^dagger chi_(j+1)
       =b_j^dagger q_j lambda_j b_(j+1)
       =b_j^dagger lambda_j b_(j+1).                      (9)

The last equality uses b_j^dagger q_j=b_j^dagger. On the auxiliary
vacuum pi_j=lambda_j=1. Thus (9) reduces exactly to the intended payload
Jordan-Wigner hopping and identity-transport hard-core length hopping.
Off that vacuum the bounded auxiliary parity factors are essential.
A naive extension by the identity on the original spin auxiliary factor
would instead produce correlated hoppings in CAR variables.

This extension can be made isometric for every separately neutral code
operator, including the complete retained interacting Hamiltonian. Let
O_X be such a code operator. First write it using the separate code-only
CAR maps, then insert the same expression using (4). Call the result
E_JW(O_X). It is diagonal in all auxiliary parity sectors. Fix the
auxiliary eigenvalues pi_j=alpha_j and lambda_j=beta_j. Relative to
the code-only CAR operators, the payload fields acquire signs

    s_j=product_(a<=k<j)alpha_k,

and the length fields signs

    t_j=product_(a<=k<j)beta_k.

There is also a spatially constant Klein sign in the length fields;
it cancels from every separately neutral O_X. Exterior prefix signs
cancel for the same reason. In this auxiliary parity sector,

    E_JW(O_X)=U_(alpha,beta) O_X U_(alpha,beta)^dagger,    (10)

where U is the diagonal logical gauge unitary multiplying a logical
occupation basis vector by

    product_(j in X) s_j^(N_payload,j) t_j^(N_length,j).

Every block in (10) is unitarily equivalent to O_X, so

    ||E_JW(O_X)||=||O_X||.                               (11)

The construction is support preserving on intervals and agrees with the
original code operator at auxiliary vacuum. It is a genuine extension,
not a coefficientwise norm bound with an exponential loss.

For the locally block-diagonal reference B_X, set O_X=P_X B_X P_X and

    K_D,X=J_X B_X J_X^dagger-E_JW(O_X).

Local block diagonality makes K_D,X annihilate the complete auxiliary
vacuum from both sides. Equations (11) and ||O_X||<=||B_X|| give

    ||K_D,X||<=2||B_X||.                                 (12)

Thus the existing small-input bound for massive matching is preserved.
Only its reference-extension convention changes. The retained Hamiltonian
in CAR variables is exactly the previously computed code Hamiltonian;
its quadratic part is (8), with n payload bands and one length band.

## 6. Logical coefficient bounds remain dimension independent

On an interval of r cells, let m=(n+1)r be the number of logical modes.
After the isometric local rewriting above, expand O_X in logical matrix
units with auxiliary-operator coefficients. There are 4^m matrix blocks,
each of operator norm at most ||O_X||. Each logical matrix unit expands
through the logical vacuum projector into at most 2^m normally ordered
monomials, each of degree at most 2m. Consequently, for H>=1,

    sum_I H^(degree I)||O_(X,I)||
       <=(8 H^2)^((n+1)r)||O_X||.                        (13)

To express a graded auxiliary coefficient in an ordinary tensor matrix
block, one may multiply its logical matrix unit by the logical parity.
On a matrix unit this only changes a sign, so it does not alter (13).
The two auxiliary parity decompositions introduce no dimension factor:
neutrality fixes each block's compensating auxiliary charge. The bound
depends on the finite logical code and support length, not dim E or the
maximum physical cell length.

For K_D, its two-sided vacuum-annihilation implies the same property for
every coefficient, by linear independence of the logical matrix-unit
basis. Completed auxiliary-vacuum polymers have zero charge in both
families, since the vacuum has zero charge and the vertices conserve both.
The effective logical action therefore retains both U(1)'s. The same
neutrality holds for an auxiliary thermal trace by charge invariance.

These facts provide the algebraic and norm inputs to the proposed
operator-valued polymer theorem. They do not prove that theorem's
connected expansion or its scalar thermal remainder estimates.

## 7. Checks and application boundary

`jw_checks.py` verifies the CAR, graded auxiliary relations, locality of
charge-exchange vertices, and the sectorwise isometric reference
extension on a two-cell example with one payload mode, one length mode,
and all four auxiliary charge parities. It uses NumPy only. These are
algebra checks, not a central-charge measurement.

The application requires both exact charges and the parity-twisted
extension (9)-(10). It would not apply unchanged if only total combined
parity were conserved, or if the coefficient operators carrying auxiliary
odd charge were treated as ungraded commuting scalars. Under the stated
hypotheses, no common-Jordan-Wigner locality or Gaussian-reference gap
remains in the matching setup.

The locality statement applies to neutral Hamiltonian terms and neutral
source insertions such as densities and currents. A physical charged
single-field source generally acquires Jordan-Wigner strings and needs
its usual string-observable treatment. Arbitrary charged physical sources
must not be silently included in a local neutral-source theorem.
