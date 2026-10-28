# Nonnegative joint Clifford codes and classical readout histories

Status: proved finite-dimensional code and history-Hamiltonian statements.
The circuit interface uses reversible permutations only. A uniform
nearest-neighbor Hamiltonian realizing the entire initialized computation
has not been supplied. These results do not settle the c-d conjecture.

## 1. A joint code with nonnegative amplitudes

Let n>=1, D=2^n, and choose 2n distinct unordered pairs of letters from
{0,...,q-1}. This is possible precisely when

    2n <= q(q-1)/2.

Label the pairs (a_{f,s},b_{f,s}), a<b, by f=1,...,n and s=+1,-1.
An arbitrary fixed enumeration is sufficient. Define two Hermitian
q-by-q matrices for each flavor,

    X_f = sum_s s (|a_{f,s}><b_{f,s}| + |b_{f,s}><a_{f,s}|),
    Y_f = sum_s s i (|b_{f,s}><a_{f,s}| - |a_{f,s}><b_{f,s}|).  (1)

The left operators are A_{2f-1}=X_f, A_{2f}=Y_f. The right operators
are B_{2f-1}=Y_f, B_{2f}=-X_f. All have zero diagonal. They are
linearly independent, since distinct flavors use distinct pairs and
the real and imaginary entries on each pair are independent.

For a logical bit string x=(x_1,...,x_n), write z_f(x) for x with bit
f deleted, and define

    s_L(f,x)=(-1)^(sum_{j<f} x_j),
    s_R(f,x)=(-1)^(sum_{j>f} x_j).

Let ell_f(x) be a_{f,s_L} if x_f=0 and b_{f,s_L} if x_f=1.
Define r_f(x) with s_R in place of s_L. The environment has orthonormal
labels |L,f,z> and |R,f,z>, hence dimension nD. Fix any letter o at an
inactive port, for example o=0, and set

    V|x> = (2n)^(-1/2) sum_f [
        |ell_f(x)>_L |o>_R |L,f,z_f(x)>
      + |o>_L |r_f(x)>_R |R,f,z_f(x)> ].                  (2)

This is an isometry C^D -> C^q tensor C^q tensor C^(nD).
Each column has 2n distinct computational basis configurations. Two
columns can share an environment label only if their bit strings differ
at its specified flavor, and then their active port letters differ.
Thus their full supports are disjoint. Every nonzero matrix entry of V
is the same positive number (2n)^(-1/2).

Use the conventional Jordan-Wigner matrices on the logical bit register,

    Gamma_{2f-1}=Z_1 ... Z_{f-1} X_f(logical),
    Gamma_{2f}  =Z_1 ... Z_{f-1} Y_f(logical),
    P=Z_1 ... Z_n,   R_a=i Gamma_a P.                     (3)

The claim is

    V^dag A_a(left) V = Gamma_a/(2n),
    V^dag B_a(right) V = R_a/(2n).                        (4)

Proof of the left map. A boundary operator cannot change the environment
label. Its diagonal contribution is zero by (1). Within the record
|L,f,z>, only the two logical states x with x_f=0 and y=x xor e_f
occur. A_{2f-1} has matrix element s_L from x to y, and A_{2f} has
matrix element i s_L. These are the corresponding Gamma matrix
elements. Each is multiplied by the product of the two amplitudes,
1/(2n). The reverse entries are their complex conjugates. Other
flavors have zero matrix element on this pair, since the pairs are
distinct. In right records the left port is o, and its zero diagonal
again gives zero. This proves the first identity.

For the right map, take x_f=0. The two desired entries are

    <y|i Gamma_{2f-1} P|x> = i s_L (-1)^(sum_j x_j) = i s_R,
    <y|i Gamma_{2f} P|x>   = -s_L (-1)^(sum_j x_j) = -s_R.

Those are precisely the entries of B_{2f-1}=Y_f and B_{2f}=-X_f
on the pair labeled by s_R. The same environment and normalization
argument proves the second identity. This also fixes the otherwise
important parity-dressing signs. QED.

It follows that the physical port bond

    v = sum_{a=1}^{2n} B_a(right) tensor A_a(left)          (5)

compresses on two codes to (2n)^(-2) sum_a R_a tensor Gamma_a.
This is the free Majorana hopping target already derived and calibrated
in RESULT.md Section 2, with c_target=n. No new CFT input is used here.

In particular, n=floor(q(q-1)/4) grows quadratically in the port
dimension. At q=5, n=5 exceeds q-1. The physical block in (2), however,
has dimension q^2 n 2^n=4000 at this example. Declaring that entire block
to be one site would not violate c-d. The number q is only its port
dimension, and there is no local q-state parent in this construction.

The flat parent Delta(I-VV^dag), Delta>0, is stoquastic in this basis.
Its off-diagonal entries are nonpositive because VV^dag is entrywise
nonnegative. It has the complete ground code (2) and gap Delta. It is
an unrestricted block Hamiltonian. This proves that neither
nonnegative code amplitudes nor stoquasticity of an unrestricted parent
forces the linear-in-q channel budget. It says nothing about
stoquasticity of the inter-block bond (5).

## 2. Exact grading

Split the q port letters into q_+ even and q_- odd letters. Require
every selected pair to have a even and b odd. There are q_+q_- such
pairs, so the same construction now permits

    n <= floor(q_+ q_- / 2).                              (6)

Every A_a and B_a is odd under the diagonal port grading omega.
Give the environment record |L,f,z> the parity

    omega_E(L,f,z) = (-1)^(sum_j x_j) omega(ell_f(x)) omega(o),

where x is either extension of z by the missing bit. The value is
independent of that choice: flipping x_f reverses both factors that
depend on it. Define the right record parity in the same way using
r_f(x). Then

    (omega_L tensor omega_R tensor omega_E) V = V P.       (7)

The flat parent is even. Thus the positive-amplitude construction
also accommodates an exact total fermion parity, at a factor-two
cost in the number of complex flavors relative to the unrestricted
pair count. For balanced grades n is about q^2/8, which still leaves
a quadratic rather than linear channel budget.

This does not contradict the earlier boundary-Schmidt-rank theorem.
Whenever n>q-1, (4) itself and that theorem imply that every nonzero
code vector is entangled across either active port. Positivity of a
vector's coefficients does not make it a product vector.

## 3. A clean readout history with an exact gap

There is an elementary history realization of the same boundary maps.
It is useful because its propagation gates are permutations, and it
does not require a quantum circuit that prepares the Choi dilation.
The following explicit version uses a global clock and global gates.
It is NOT the missing uniform nearest-neighbor parent.

Add one inactive letter * to each port, and extend every A_a and B_a
by zero on that letter. The data workspace has basis |x>, x in {0,1}^n.
Define baseline configurations

    b_x=|*>_L |*>_R |x>.

For each flavor f define left and right readout configurations

    a^L_{f,x}=|ell_f(x)>_L |*>_R |z_f(x)>,
    a^R_{f,x}=|*>_L |r_f(x)>_R |z_f(x)>.

The shorter bit string z is stored in the first D/2 workspace states.
For each readout label (u,f), let U_{u,f} swap b_x with a^u_{f,x}
for every x, and fix every other basis state. These transpositions
are disjoint. U_{u,f} is a real permutation, Hermitian and unitary.

Execute each U_{u,f} and immediately its inverse, in a fixed order
over the 2n readouts. There are T=4n steps and C=T+1 clock values.
The history at even clock values is the baseline, and each odd clock
value is one designated readout. With W_t=U_t ... U_1 and W_0=I, use

    H_prop = (1/2) sum_{t=0}^{C-2} [
       |t><t| tensor I + |t+1><t+1| tensor I
       - |t+1><t| tensor U_{t+1} - |t><t+1| tensor U_{t+1}^dag ],
    H_init = |0><0| tensor (I-P_baseline).                 (8)

Every term is positive semidefinite. The complete zero space is

    V_hist|x> = C^(-1/2) sum_{t=0}^{C-1} |t> W_t b_x.      (9)

To see completeness, conjugate by the controlled unitary
W=sum_t |t><t| tensor W_t. The propagation becomes the path
Laplacian L_C=(1/2)sum_t(|t>-|t+1>)(<t|-<t+1|), while H_init is
unchanged. A path zero mode is constant, and its data vector must
lie in the baseline space. This proves that there are exactly D
zero states, not just D constructed ones.

The same conjugation gives the entire spectrum. If K=(q+1)^2 D is
the data-space dimension, the two blocks are

    L_C,                         multiplicity D,
    L_C+|0><0|,                  multiplicity K-D.

Their eigenvalues, for j=0,...,C-1, are respectively

    1-cos(pi j/C),
    1-cos(pi(j+1/2)/C).                                   (10)

For the first use the usual cosine eigenvectors
cos(k(t+1/2)) with k=pi j/C. For the second use
sin(k(t+1/2)): the left endpoint equation is the odd-reflection
condition u_{-1}=-u_0, while the right equation is u_C=u_{C-1}.
The latter requires cos(kC)=0. Both give the displayed eigenvalues.
The gap is consequently

    Delta_hist=1-cos(pi/(2C)) > 0.                        (11)

Every intermediate gate and initial vector has nonnegative entries,
so (9) has nonnegative amplitudes. The off-diagonal entries of (8)
are nonpositive. Thus this is an exact frustration-free stoquastic
history parent, with its full ground multiplicity and gap known.

Since the clock is orthogonal, compression of a port operator is
the average of its compressions at each time. It annihilates the
baseline configurations. At its unique readout for flavor f its
compression is the specified Gamma or R. Hence

    V_hist^dag A_a(left) V_hist = Gamma_a/C,
    V_hist^dag B_a(right) V_hist = R_a/C.                  (12)

Unlike a computation that assumes the port operator is unchanged at
all times, this history intentionally uses different port actions at
different times. There is no limiting clock-coupling error in (12).
Longer computation dilutes the amplitude but does not distort the
normalized Clifford maps if only clean readouts activate the ports.

## 4. A reversible circuit interface using a finite gate menu

The global swaps in (8) are not nearest-neighbor gates. Nevertheless,
their action on the initialized inputs admits a clean local circuit
interface. The grading version is especially simple. This statement
concerns a scheduled circuit, not an autonomous uniform Hamiltonian.

Take a workspace bit register, parity-even auxiliary bits, a q-state
partner next to each q+1-state port, and an inactive port letter * of
even parity. Choose reference partner letters e_0 and o_0 of even
and odd parity. For a left readout of bit f do the following.

1. Reversibly compute sum_{j<f} x_j mod 2 into an auxiliary bit whose
   two states both have even grading. The input bits are controls
   and are unchanged. For a right readout use sum_{j>f} x_j instead.
2. Move x_f to the partner neighborhood by adjacent swaps. Transfer
   it into the partner by exchanging the two states (1,e_0) and
   (0,o_0), leaving (0,e_0) fixed. The logical bit position is now
   zero, and the partner letter is e_0 or o_0. The exchange preserves
   total grading.
3. Conditioned on the computed sign, rotate the even letters to send
   e_0 to a_{f,s}, and the odd letters to send o_0 to b_{f,s}. Powers
   of one cyclic permutation on each grade suffice. This uses at
   most q_++q_- cycle operations for each sign branch. The sign
   depends only on z_f(x), so no other register remembers x_f.
4. Exchange the partner and inactive port, extending the partner by
   the same inactive letter if needed. All other registers now
   depend only on z_f(x). Idle for the same number of clock slices
   at every readout, exchange back, and reverse the computation.

While the computation is running the physical port is *, so every
chosen boundary operator is zero there. During readout its only
logical-bit dependence is the designated letter a or b. Consequently
each readout slice contributes precisely the matrix elements used
in Section 1. If there are w active slices per readout in a complete
history of C slices, the compression scale is w/C at both ends.

All elementary operations are permutations and preserve total grading.
Parity can be accumulated by a moving even auxiliary bit, and data
can be routed by adjacent swaps, so the scheduled circuit needs only
local gates and a number of auxiliary registers polynomial in n and q.
A fixed menu suffices: swaps, controlled toggles of an even auxiliary
bit, the transfer above, controlled cycles within each grade, and
the boundary exchange. Additional classical controls can be stored
in a constant number of bits per gate location. No quantum phase
gates or general unitary synthesis are needed.

The fixed sequence selecting f, its sign branch, its cyclic-shift
counts, and the return sweep is still external in this description.
Replacing it by a repeated nearest-neighbor Hamiltonian requires a
reversible controller with counters, initialization, and a proof that
all malformed histories lie above the code. The local boundary
exchange itself also has to be gated by that controller. No theorem
in this note asserts that those obligations have been met.

For the global-clock Hamiltonian, inserting such a finite permutation
circuit and initialization projectors at t=0 preserves the proof of
the complete history code. The simple spectrum (10) applies with a
single projector onto the allowed input space; a sum of commuting
initialization projectors only increases the invalid-sector penalty
and retains the lower bound (11). A physical spatial clock is a
separate construction.

## 5. Resource accounting and the remaining bottlenecks

The static isometry uses environment dimension n2^n and total dimension
q^2 n2^n. Ignoring the missing parent, it can fit into k q-state sites
whenever q^(k-2)>=n2^n. The nonnegative isometry is therefore no more
expensive in leading block length than the earlier abstract dilation.

With an exact port grading it supports n<=floor(q_+q_-/2). If a future
controller and selector have a fixed dimension cost K per site, and
the physical alphabet is K(q+1), the formal inequality

    floor(q^2/8) > K(q+1)-1

holds for sufficiently large even q. This is only resource accounting.
Neither the constant K nor the required uniform parent has been
constructed here. The inequality cannot be reported as a microscopic
counterexample.

The circuit route narrows one open step. A classical reversible
controller with clean boundary exchanges could implement the needed
first-order ports; it need not realize arbitrary quantum gates or
the entire prescribed excited-state spectrum of the second-order
gadget. But three substantial issues remain even after such a
controller is written: selection of the correct initialized cell
sectors, conformal infrared behavior at fixed nonzero inter-cell
coupling, and compatible payload transport on all ring lengths.
The exact code identities do not settle any of these three.

## 6. Numerical checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 classical_readout_checks.py

The script writes classical_readout_checks.json. It constructs the
isometry as an actual sparse matrix, applies the physical boundary
operators, and compares both maps and the complete projected bond.
Cases are (q,n)=(3,1),(4,3),(5,5),(6,6), and the exactly graded
(4,2),(6,4). These include two static codes with n>q-1. Maximum
operator errors are below 3e-16, and the maximum bond Frobenius error
is below 1e-15. The graded intertwining errors vanish at this precision.

The global-clock Hamiltonians at (q,n)=(3,1),(4,2) have dimensions
160 and 900. Their complete spectra agree with (10) within 1.2e-14,
their ground multiplicities are 2 and 4, and their gaps agree with
(11). The code and both physical port maps agree within 4e-16.
These are algebraic checks, not central-charge estimators. The
large-q resource examples and the proposed spatial controller have
not been simulated.

## 7. Independent flavor densities from diagonal two-site readouts

A stable-target proposal needs adjustable Fermi momenta in addition to
the symmetric Majorana hopping. The graded code can supply independent
logical Z_f fields through diagonal operators on a port and its neighbor.
This is another code statement, not a complete local parent.

Take q=2s with even letters a=0,...,s-1 and odd letters b=s,...,2s-1.
Color the complete bipartite graph by

    color(a,b)=(a+b-s) mod s.                             (13)

At a fixed vertex no two incident edges have the same color. Thus the
pair (vertex,color) identifies an edge uniquely. This remains true for
any selected subset of the edges, including the 2n signed pairs.

Augment each readout record in (2) by a register next to its active
port carrying the color of its pair. The inactive neighbor carries a
separate inactive color. The color is a function of f and s_L or s_R,
hence of z_f(x), and is unchanged when x_f flips. All proofs of the
Clifford maps therefore survive: they only used equality of the
remaining record for those two logical inputs.

Define a diagonal operator D_f on the left port and color register.
For each of the two signed pairs belonging to flavor f, give its
(a,color) slot eigenvalue +1 and its (b,color) slot eigenvalue -1.
Set all other slots, including the inactive color, to zero. Proper
edge coloring ensures these assignments are consistent and that the
supports of D_f and D_g are disjoint for f!=g. Thus

    V_col^dag D_f V_col = Z_f/(2n),
    D_f D_g=0 for f!=g,   ||sum_f mu_f D_f||=max_f |mu_f|. (14)

Indeed the left flavor-f readout sees +1 for x_f=0 and -1 for x_f=1,
with probability 1/(2n). Every other left flavor and every right
readout sees zero. This proves the compression and the norm statement.
The operators are diagonal and preserve the exact grading.

The color register is a separate neighboring site, not a second
q-dimensional factor on the port site. In the scheduled circuit of
Section 4, the partner can display the color after the logical bit
has been transferred to the port. For balanced grades there are s
colors and s even partner letters, so this display can preserve
parity. It must be enabled only at completed readouts, for example
by a fixed-size local controller flag; during partial computation
(14) is not automatically the compression of an always-on field.
This is an additional explicit obligation for a spatial controller.

To see which target is accessible, let alpha=1/(2n) for the static
colored code. With global Majoranas gamma_{j,a}, the already-proved
identity R_{j,a} Gamma_{j+1,b}=i gamma_{j,a} gamma_{j+1,b} converts
an arbitrary real coefficient matrix T in the boundary bond into
an arbitrary real inter-site Majorana hopping matrix. For one
complex flavor c_j=(gamma_{j,x}+i gamma_{j,y})/2,

    -t(exp(-i theta)c_j^dag c_{j+1}+h.c.)
      = (it/2)[sin(theta)(gamma_{j,x}gamma_{j+1,x}
                            +gamma_{j,y}gamma_{j+1,y})
                -cos(theta)gamma_{j,x}gamma_{j+1,y}
                +cos(theta)gamma_{j,y}gamma_{j+1,x}].      (15)

Also -mu c_j^dag c_j=(mu/2)Z_j-mu/2. Consequently a physical boundary
bond lambda sum_ab T_ab B_a tensor A_b, together with the intra-cell
field lambda alpha sum_f (mu_f/2)D_f, compresses to lambda alpha^2
times the desired hopping and independent chemical potentials, up
to a scalar. Formula (15) gives T explicitly, flavor by flavor.
For history codes the field and port compression amplitudes may
differ; their known scalar factors can be compensated in the field
coefficient. This removes a separate controllability issue at first
order. It does not prove stability at finite lambda.

The script additionally checks the colored codes at (q,n)=(4,2),(6,4).
It applies the actual boundary matrices and diagonal two-site fields
to the augmented sparse isometry. Both port identities and all
independent density identities hold below 3e-16, the density supports
are disjoint, and the field norm agrees with (14). These tests are
included in classical_readout_checks.json.
