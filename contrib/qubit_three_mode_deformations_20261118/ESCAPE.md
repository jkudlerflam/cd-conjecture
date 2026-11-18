# A number-changing escape: explicit six-term family and bounded test

Status: **[exact microscopic specification]**, **[numerical finite leading-
model diagnostics]**. This escapes the hypotheses of the saturation
certificate by taking a nonzero parent displacement together with a
normal hopping deformation and a finite chemical potential. It is not
another examination of the alpha-only zero-chemical parent kernel.

Choose the following dimensionless controls of the complete operator
(7) in [NORMAL_SPACE.md](NORMAL_SPACE.md):

```
a=2, zeta=1/100, eta=0, mu=1/2.                             (1)
```

The D3 coefficient is now **+2**, a repulsion, and the indispensable
charge-changing facilitation has amplitude **2**. Both are kept. This
is a useful sharp test of the tempting repulsion argument: it retains
the entire physical companion while placing the chemical potential
above all three bare one-particle wells. Real zeta keeps the literal
microscopic Pauli count smaller than a chiral hopping would.

## 1. Literal nearest-neighbour qubit realization

Set alpha=16, beta=0 and all higher pair jets to zero. For any fixed
positive epsilon, define

```
H_micro = D1+epsilon F+u N+[t T1+p A+h.c.],
p = 3 epsilon²+16 epsilon³,
u = 2 epsilon²-34 epsilon⁴-256 epsilon⁵-(1012/3)epsilon⁶,
t = -epsilon²+8 epsilon⁴+32 epsilon⁵-(5968/75)epsilon⁶.      (2)
```

This is a literal polynomial bond, not an effective range-three
microscopic interaction. The effective scale is `(128/3)epsilon⁶`.
Epsilon and every coefficient are fixed before L grows. The odd fifth
jets in (2) are essential; omitting them would invalidate its leading
operator. Equation (2) is Hermitian and one-site translation invariant
at physical d=2, with no markers, blocking or encoded dimension.

After constants and telescoping terms are removed, its six literal
Pauli monomials are

```
(1/4) ZZ -(1+u)/2 Z +(epsilon/2) XZ -(epsilon/2) ZX
                         +(t+p)/2 XX +(t-p)/2 YY.           (3)
```

All entries denote uniform sums. The real family obtained before fixing
(1) has **four independent dimensionless controls** epsilon, real a,
real zeta and mu (beta=eta=0 and the overall scale removed), not six
independent coefficients. The slice (1) has only epsilon remaining.
The two oriented XZ and ZX terms count as two literal monomials.

| Fixed epsilon | Max/min nonzero absolute Pauli coefficient | Sixth kinetic scale |
|---|---:|---:|
| 1/20 | 142.3341640 | 6.6666667e-7 |
| 1/5 | 4.6104533 | 0.0027306667 |

The exact fractions and a literal-bond/Pauli identity check are in
`escape_results.json`. Common positive rescaling does not change these
ratios. The smaller spread at 1/5 is attractive as a simple Hamiltonian,
but the leading expansion is **not controlled there**. Neither listed
epsilon has a uniform thermodynamic remainder certificate. They are
algebraic realizations, not experimentally or numerically established
three-mode phases. A small emergent scale is distinguished from the
actual microscopic coefficient spread.

## 2. Complete-sector finite check

`escape_check.py` diagonalizes the entire leading Hamiltonian at (1) on
L=12,15,18, dimensions 58,159,442. It does not impose a particle-number
sector. It averages observables over the full numerically degenerate
ground space. N acts as a scalar within that ground space to the stated
rounding precision, making the reported branch response unambiguous.

| L | Mean density | Var(N) | Single-residue probability, conditional on N>=2 | Finite chemical response per site |
|---|---:|---:|---:|---:|
| 12 | 0.08333333 | 0 within rounding | no N>=2 weight | 0 within rounding |
| 15 | 0.11124042 | 0.44882674 | 0.9999999144 | 0.01075196 |
| 18 | 0.09765526 | 0.52959544 | 0.9999997825 | 0.01116098 |

Here a "single-residue" configuration places every occupied site in
one congruence class modulo three. The observable includes superpositions
of the three possible classes; it is not an assertion of a symmetry-broken
finite-ring eigenvector. Conditioning on N>=2 removes the trivial fact
that every one-particle configuration has a single residue. The N>=2
probabilities are 0.5615573 and 0.6005834 at L=15 and 18, respectively.
Residue polarization is distinct from momentum-valley polarization:
a state on one residue is coherent across the three band centers.
This observable does not measure individual valley populations or their
stiffnesses; it tests whether the finite state leaves the known collective
single-residue structure.

The chemical response is the spectral derivative

```
(2/L) sum_(excited e) |<e|N|ground>|²/(E_e-E_ground),         (4)
```

averaged over the ground multiplet. It is a finite susceptibility.
Because N is not conserved, it is **not** a thermodynamic Luttinger-liquid
compressibility or evidence of a positive sound velocity. The full
commutator [H,N] has matrix elements of magnitude 2 in every tested ring;
particle-number fluctuations have not been discarded.

For a low-filling control with the same a,zeta but mu=3/100, all three
sizes have only the two one-particle W ground states at energy -0.04.
They have zero finite response. The exact momentum cancellation in
NORMAL_SPACE.md explains why these W states can remain eigenvectors
despite [H,N]!=0. This control is not evidence for a thermodynamic
one-particle phase: the dilute crossover can exceed the tested rings.

The new escape avoids the leading saturated crystal on the two larger
tested rings, but it does **not** pass a three-component-mixture test.
Its ground support is almost entirely in the single-residue span. These
small rings do not establish a thermodynamic polarized phase either.
No c, stiffness matrix, absence of higher bound clusters, or irrelevant
locking exponent is inferred.

## 3. Decisive next calculation

The next useful analytic calculation is the energy and hopping of a
**mixed-residue interface in the finite-density a!=0 model**, including
all creation/annihilation processes. It should determine whether a
physical complex zeta can make a mixture favorable before an allowed
pair or relative-phase locking vertex grows. The forward and relative
couplings and the locking amplitudes must come from that same (7), not
an independently chosen Luttinger matrix or a fixed-N truncation.

For real zeta the pair of opposite valleys remains momentum matched,
as does `(b†_0)^2 b_+ b_-` if a third valley is present. Thus even a
positive mixture susceptibility would not finish the c>=3 construction.
If the interface cost stays positive under all admissible small normal
directions, that would extend the present obstruction. If it can be
reversed with three positive stiffness eigenvalues and irrelevant matched
locking, it would identify an actual interacting mechanism worth a later
converged calculation. No additional production calculation was launched.
