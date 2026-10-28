# Exact parity for the second-order route-A gadget

## Claim

Let the logical cell be the Fock space of n complex fermions, with
N = 2^n and parity P. Any Hermitian two-cell target K satisfying
[K, P tensor P] = 0 can be synthesized at second order through two
qubit ports by identical finite-dimensional cells, with all of the
following properties:

1. The cell has an exactly N-dimensional ground code and a positive gap.
2. Every traceless single-port operator compresses to zero on that code.
3. A product parity U_cell = Z_left tensor Z_right tensor U_aux commutes
   with the complete cell Hamiltonian and restricts to P on the code.
4. The inter-cell bond is the ordinary Pauli exchange
   V = Z_right Z_left + X_right X_left + Y_right Y_left. It preserves
   the global product parity exactly, at every perturbative order.
5. Its complete second-order chain Hamiltonian is a scalar plus a
   positive coefficient times the desired K on every neighbouring pair.
   Cross terms between different bonds vanish at this order.

This is an abstract-cell statement. It does not supply a microscopic
one-site translation-invariant nearest-neighbour implementation within
each cell, nor establish a CFT. The logical free-fermion target would
have Virasoro central charge c = n when its bands have two simple Fermi
points; the gadget theorem itself makes no claim about the infrared
limit at finite coupling.

## Status

[proved] The algebraic theorem below extends C9's spectral-weight
construction to exact parity. [numerical] An independently assembled
physical resolvent verifies it on two cells, including mixed even and
odd target operators, and on three cells for an odd target. The
fixed-alphabet local implementation and finite-coupling CFT remain open.

## How it was established and checked

- Author: OpenAI Codex, GPT-6 family (exact runtime model identifier was
  not exposed to this worker), 2026-10-28.
- Inputs read: C9 RESULT, C8 RESULT and PROGRESS, ir_control README,
  minimal_d0 README, C13 RESULT, C18/C19/C20 specifications, and STATUS.
- Checks: direct matrix construction and resolvent evaluation in
  `checks.py`; no human review or independent-agent review is claimed.
- No external theorem is needed for the finite-dimensional proof.

## 1. Why the diagonal C9 recipe needs modification

Give each port the grading Z. Its Pauli channels

    F_0 = Z, F_1 = X, F_2 = Y

have parity signs tau = (+1, -1, -1). If the cell Hamiltonian preserves
parity and its code intertwines the physical and logical gradings, its
spectral weights necessarily satisfy

    P M_e[a,b] P = tau_a tau_b M_e[a,b].                 (1)

This follows by conjugating G^dagger F_a Pi_e F_b G. In particular,
every diagonal weight M_e[a,a] is logically even. C9 Section 2 puts
all target operators into diagonal weights, so that prescription
cannot put parity-odd hopping factors into an exactly graded cell.
This is a limitation of that prescription, not an impossibility
theorem for route A. Off-diagonal even-to-odd syndrome weights repair it.

## 2. A parity-compatible positive spectral-weight construction

Choose Hermitian parity-homogeneous lists A_i and B_i, i = 1,...,m,
with parities s_i in {0,1}, and a real matrix K_ij such that

    K = sum_ij K_ij A_i tensor B_j,
    K_ij = 0 whenever s_i != s_j.                         (2)

Such lists always exist, by expanding a Hermitian globally even
operator in even and odd Hermitian bases. Missing basis elements can
be padded with zero coefficients. In syndrome space C^3 set

    T_0 = |0><0|,
    T_1 = (|0><1| + |1><0|)/sqrt(2).                    (3)

These are real symmetric, obey

    sum_ab (T_s)_ab (T_t)_ab = delta_st,
    diag(tau) T_s diag(tau) = (-1)^s T_s.               (4)

Choose k = m+2 distinct positive energies e_e and probabilities p_e>0.
Let C_ef = 1/(e_e+e_f). C is positive definite because

    x^T C x = integral_0^infinity (sum_e x_e exp(-e_e t))^2 dt.

The subspace orthogonal to both 1 and Cp has dimension at least m.
Consequently there is a real k-by-m matrix Zc satisfying

    1^T Zc = 0, p^T C Zc = 0, Zc^T C Zc = I_m.          (5)

On each 3N-dimensional port syndrome space, define

    M_e^R = p_e I + rho sum_i (Zc)_ei T_{s_i} tensor A_i,
    M_e^L = p_e I - rho sum_j (Zc K)_ej T_{s_j} tensor B_j. (6)

For each fixed finite target choose rho>0 small enough that both
perturbations in (6) have norm at most p_e/2, for every e. Then
M_e^L and M_e^R are positive, their sums are identity, and each
commutes with diag(tau) tensor P. Their direct sum M_e on the two
port syndrome spaces is a parity-covariant POVM.

The contraction in second-order perturbation theory is entrywise,
not a matrix product in syndrome labels:

    H2 = -sum_ef C_ef sum_ab M_e^R[a,b] tensor M_f^L[a,b]. (7)

The baseline contributes -3 p^T C p times identity. Terms linear
in rho vanish by (5). Equations (2), (4) and (5) make the remaining
term +rho^2 K, with no missing sign or transpose. Therefore

    H2 = -3 p^T C p I + rho^2 K.                         (8)

For actual cell energies Delta e_e and inter-cell coupling lambda V,
multiply (8) by lambda^2/Delta. Complex Hermitian logical matrices,
such as Y, cause no problem: the *syndrome* matrices T_s are real
symmetric. This is why their contraction in (4) has the required sign.

## 3. Symmetry-equivariant realization by a physical cell

Use the same Bell-pair code as C9:

    G|i> = (1/2) sum_{a,b=0,1} |a>_L |b>_R |a,b,i>_aux. (9)

The auxiliary space contains these 4N basis states, plus unused states.
Set U_aux|a,b,i> = (-1)^(a+b) P_i |a,b,i>, where P_i is the
logical parity eigenvalue. Assign either sign to any unused auxiliary
state. Then U = Z_L Z_R U_aux obeys UG = GP. For

    S = [F_0^L G,F_1^L G,F_2^L G,F_0^R G,F_1^R G,F_2^R G],

the Bell identities give

    G^dagger S = 0, S^dagger S = I_{6N},
    U S = S T,
    T = (diag(tau) tensor P) direct_sum (diag(tau) tensor P).

Let W x = direct_sum_e sqrt(M_e) x. It is an isometry and obeys
(direct_sum_e T) W = W T. Both parity sectors of its source have
dimension 3N, and both sectors of its dilation have dimension 3kN.

The physical U has D/2 states of each parity, independently of the
unused auxiliary signs, because Z_L Z_R itself has equal multiplicities.
The code uses N/2 states of each parity. Thus, if

    D >= 16N, D-N >= 6kN, and 4 divides D,               (10)

the complement of the code contains at least 3kN states of each parity.
Extend the isometry W -> S separately in the plus and minus sectors.
This gives an isometry J from the dilation into the physical code
complement, with

    J W = S, UJ = J (direct_sum_e T), J^dagger G = 0.

This is an explicit equivariant Naimark extension, not an assumption
that an arbitrary ungraded extension happens to preserve symmetry.
The script constructs it by orthogonal completion in each parity sector.

With Pi_e the coordinate projector of the dilation, set

    H_cell = Delta sum_e e_e J Pi_e J^dagger
             + Delta e_* (I-GG^dagger-JJ^dagger), e_*>0. (11)

It has precisely the code (9) at energy zero, a positive gap, commutes
with U, and realizes S^dagger J Pi_e J^dagger S = M_e. The physical
inter-cell Pauli exchange preserves U tensor U because both copies
of each F_a have the same parity. Equation (1) is exact for (11).

For different bonds, their second-order cross term vanishes: each
outermost cell is acted on just once and G^dagger F_a G = 0. This
includes adjacent bonds. Hence (8) is the full second-order chain
Hamiltonian rather than merely its isolated two-cell contribution.

## 4. Linear-in-n target rank and logarithmic workspace overhead

An n-flavour quadratic nearest-neighbour target has Majoranas
gamma_a, a=1,...,2n, and parity-dressed right factors

    R_a = i gamma_a P.

They are Hermitian and odd. Any parity-even quadratic hopping bond,
including independent hopping phases and flavour mixing, can be
written as sum_ab t_ab R_a tensor gamma_b with real t_ab. Add an
arbitrary even quadratic on-cell term D by distributing it as
(D tensor I + I tensor D)/2. Equation (2) then needs at most

    m = 2n+2, k = 2n+4,
    D_cell = 4 ceil(((12n+25) 2^n)/4).                   (12)

Thus log_2 D_cell = n + O(log n), rather than the 3n+O(1)
qubits of a complete operator basis in C9's general bound. This
does not make the controller uniform by itself: its matrix entries,
spectral-weight normalization, and equivariant completion still have
to be generated and implemented by local rules.

For example, c_j^dagger c_{j+1} = c^dagger P tensor c under
Jordan-Wigner, and

    t exp(i theta) c^dagger P tensor c + h.c.
      = (t/2) [sin(theta)(R_x tensor gamma_x + R_y tensor gamma_y)
                + cos(theta)(R_x tensor gamma_y - R_y tensor gamma_x)].

Independent theta values and chemical potentials therefore fit (12).
Band crossings and the stronger joint nonresonance conditions of
ir_control must still be imposed separately. The examples in the
script test synthesis, not those infrared conditions.

An exact parity on an abstract auxiliary register is also not yet the
required repeated on-site parity of a fixed-alphabet physical chain.
A later compiler must realize the auxiliary grading as a product of
local gradings and use parity-preserving gates. Padding and encoding
can accommodate a finite parity representation; no local compiler is
claimed here.

### 4.1 Closed-form square roots for the hopping-only weights

For hopping alone take m=2n and let all A_i and B_i be Clifford
generators, respectively R_i and gamma_i. Each nontrivial weight then
has the form

    M = p I + rho T_1 tensor Gamma(v),
    Gamma(v) = sum_i v_i Gamma_i, v_i real.

The Clifford relations give Gamma(v)^2 = ||v||^2 I. Set

    E = |0><0|+|1><1|, eta = ||v||/sqrt(2),
    a = [sqrt(p+rho eta)+sqrt(p-rho eta)]/2,
    b = [sqrt(p+rho eta)-sqrt(p-rho eta)]/(2 eta).

For p>rho eta and v != 0, its positive square root is exactly

    sqrt(M) = [sqrt(p)(I-E)+a E] tensor I
                + b T_1 tensor Gamma(v).                 (13)

For v=0 it is sqrt(p)I. Equation (13) follows by splitting into
the T_1=0 sector and the two eigenvalues +/-eta of T_1 tensor
Gamma(v). The minus sign in a left-port weight is included by
replacing v with -v. This eliminates a generic matrix square root
on an exponentially large logical register.

There is also a constructive parity-preserving alignment using
2n-1 adjacent Majorana rotations. Starting from j=2n and working
down to 2, let theta_j = atan2(v_j,v_{j-1}) for the current vector
and conjugate by

    U_j = exp[(theta_j/2) Gamma_{j-1} Gamma_j].

Its action on those coefficients is
(v_{j-1},v_j) -> (sqrt(v_{j-1}^2+v_j^2),0). The product U therefore
satisfies U Gamma(v) U^dagger = ||v|| Gamma_1. Every generator is
even. For Gamma_i=R_i it equals gamma_{i-1} gamma_i, just as for
Gamma_i=gamma_i. Adjacent Majorana rotations are on-site or
nearest-neighbour two-mode Gaussian gates in the conventional
Jordan-Wigner ordering. Thus each alignment requires O(n) such
gates. Controlled uses and the complete Naimark extension are
additional operations; this lemma is not a finished circuit compiler.

The k-by-m coefficient matrix Zc is still n-dependent, and a local
initialized controller must still generate, apply and uncompute the
required rotations and energy labels. Equations (12)-(13) provide
an algorithmically uniform finite-cell prescription. They do not
provide fixed-alphabet spatial transition rules.

### 4.2 Supplying chemical potentials separately

Chemical potentials need not be included in these hopping POVMs.
Given an even logical on-cell operator D, add the abstract cell term

    (lambda^2/Delta) G D G^dagger.

It preserves exact parity and contributes D to the logical chain
at the same order as (8). It leaves the order-lambda^2 exchange
coefficient unchanged: it has no ground-to-excited matrix element,
and its insertion into virtual denominators affects higher orders.
The isolated cell code is intentionally split by the desired D;
the unperturbed H_cell of (11) still has its exact degeneracy.

This observation permits the simpler m=2n hopping weights and then
independent chemical-potential detuning. It does not implement
G D G^dagger as a microscopic local term. That remains a compiler
obligation, just as the hopping gadget does.

## 5. Checks and limitations

`checks_results.json` records these Frobenius-norm errors:

| Target | N | D_cell | Verification | Error |
|---|---:|---:|---|---:|
| X tensor Y | 2 | 40 | two-cell full resolvent | 4.7e-15 |
| X tensor Y | 2 | 40 | three-cell full resolvent | 8.7e-15 |
| X tensor Y + 0.3 Z tensor Z | 2 | 52 | two-cell full resolvent | 3.3e-15 |

The explicitly constructed cell parity commutators and first-order
port compressions are zero in these checks. The equivariant POVM
realization errors are below 9e-16; the gap is 1 within 8e-14.
Free-target spectral weights are checked at n=1,2,3. At n=3 the
exchange coefficient is only 1.7e-10, so subtracting the scalar
offset limits relative numerical precision to about 2.2e-6. This is
reported rather than mistaken for high relative precision.
The Jordan-Wigner hopping identity in Section 4 is also checked
directly, including the hopping phase and Hermitian conjugate.
Equation (13) and the adjacent-Majorana alignment are independently
compared with matrix multiplication and diagonalization at n=1,2,3,4.

No finite-coupling central charge is extracted. Exact microscopic
parity permits interacting higher-order corrections; it is a
necessary compatibility condition for the proposed stability argument,
not a proof of stability. Neither a uniform alphabet nor a microscopic
nearest-neighbour realization follows from a whole-cell dimension.

## Reproducing

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py` in this
directory with NumPy installed. The script writes `checks_results.json`.
Only NumPy and the Python standard library are needed. The recorded
run took under one second; memory is dominated by the 64,000-dimensional
three-cell vectors. It constructs the full unperturbed resolvent by
diagonalizing each physical cell, without a truncation of its excited
space, and checks the projected second-order answer.

## Relation to the rest of the repository

The follow-up [binary-clock interface audit](binary_clock_interface.md)
extends its initialization theorem to commuting format checks and proves
the precise edge-use budget. It explains why an exponential clock path
cannot automatically be treated as repeated complete computation sweeps.
The [full-sweep interface](sweeping_clock_interface.md) formalizes the
existing fifteen-letter remedy, arbitrary ordered two-site sweep gates,
computational acceptance pins, and the remaining endpoint-readout contract.

This repairs the exact-grading compatibility gap between C9 Section 2
and ir_control hypothesis H1, and supplies a smaller target-specific
dimension bound for C19. It strengthens the algebraic case for route A.
It does not close the initialized fixed-alphabet controller problem,
the boundary-local simulation problem, or the finite-coupling infrared
and mobile-cell composition problems.
