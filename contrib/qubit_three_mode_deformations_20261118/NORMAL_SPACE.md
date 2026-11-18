# Complete physical normal space at the three-residue scale

Date: 2026-11-18. Status: **[proved local chart and exhaustion]**, using
the previously derived sixth-order local identity; **[numerical full-matrix
verification]** of that identity with all the additional normal controls.
All physical bonds are Hermitian, microscopic nearest-neighbour, uniform
under translation by one site, and independent of system size. No CFT is
inferred from the formal expansion.

## 1. Count physical coordinates before choosing a spin frame

Write `b†=|1><0|`, `n=b†b`, `N=sum n`, `D_r=sum n_j n_(j+r)`,
`T_r=sum b†_(j+r)b_j`, `A=sum b†_j b†_(j+1)`, and `B=sum b†_j`.
Every Hermitian two-qubit bond, modulo a constant and the three
telescoping terms `K tensor I-I tensor K`, is uniquely represented by

```
v D1 + u N + [t T1+p A+f B
             +sum (c_left n_(j-1)+c_right n_(j+1)) b†_j + h.c.].   (1)
```

There are **12 real controls**, equivalently nine exchange Pauli sums
and three fields. Dividing by a common positive energy scale leaves
**11**, not eight. A common onsite spin rotation supplies three additional
unitary equivalences. They are not constants or telescoping terms; we
remove them explicitly only to perform the calculation in a convenient
frame. They can be restored without changing any spectrum or microscopic
range.

Near `D1`, align the isolated eigenaxis of the symmetric exchange matrix
with z. This sets `c_left+c_right=0`. When the remaining antisymmetric
facilitation is nonzero, a common z rotation makes it positive real.
Normalize `v=1` and use it as the analytic parameter `epsilon`:

```
H = D1 + epsilon F + u N + [t T1+p A+f B+h.c.],
F = sum (n_(j-1)-n_(j+1))(b†_j+b_j).                         (2)
```

This leaves eight intrinsic real coordinates: epsilon, u, complex t,
complex p and complex f. The chart assumes the nonzero-facilitation
analytic branch; it does not classify arcs with a different leading
valuation or every qubit bond. The gauge construction and completeness
of (1) were proved in the
[antecedent classification](../qubit_cancellation_classification_20261116/CLASSIFICATION.md).

## 2. An exact five-dimensional normal complement, not only a jet count

For every sufficiently small epsilon and p, let z be the positive local
solution near one of

```
A_parent = sqrt(z) epsilon (z-p*)/(z²-|p|²),
1 = z+|p|²/z-2|A_parent|²,
t_parent = -|A_parent|²,
u_parent = 2(|A_parent|²-|p|²/z).                            (3)
```

The derivative of the scalar equation with respect to z at epsilon=p=0,
z=1 is one. The implicit-function theorem therefore gives a unique real
analytic solution. The rank-one parent with forbidden vector
`(p*/sqrt(z),-A_parent,A_parent,sqrt(z))` has exactly the same v, epsilon
and p as (2), up to the known constant and telescoping bond. Consequently
the following is an **exact local decomposition of the whole chart**:

```
H = H_parent(epsilon,p)
    +delta_u N+[delta_t T1+f B+h.c.].                        (4)
```

Thus the normal complement to the full parent manifold has **five real
controls**: real delta_u, complex delta_t and complex f. A density
interaction independent of these is absent. The two real p directions
are tangent to the full parent manifold, although transverse to its
one-dimensional exact m=3 resonance curve. Together they and (4) give
the seven intrinsic directions transverse to that curve. Uniform spin
rotations are tangents to the corresponding physically equivalent parent
family, so none of its normal directions was lost by choosing (2).

This is the useful new exhaustion argument: varying all twelve physical
bond coordinates does not supply an unexamined sixth-order repulsion
control. It is not an obstruction to higher-order or other-background
mechanisms.

## 3. Preserve the complete lower kernels

The hypotheses here are specific: retain the entire no-11 space at order
zero, the entire no-1,2 space at order two, and the entire no-1,2,4 space
at order four; remove scalar energies and leave the final space scalar
through order five. We do not merely retain a few selected particle
states. These are precisely the lower constraints of the m=3 parent.

In the chart (2), these hypotheses require

```
p = 3 epsilon²+alpha epsilon³+beta epsilon⁴+O(epsilon⁵),
delta_u = -(128/3) mu epsilon⁶+O(epsilon⁷),
delta_t =  (128/3) zeta epsilon⁶+O(epsilon⁷),
f       =  (128/3) eta epsilon⁶+O(epsilon⁷).                (5)
```

Here mu is real; alpha, beta, zeta, eta are complex. To see exhaustion,
suppose the first nonzero coefficient in the last line of (4) occurred
before order six. Its direct compression has three independent
witnesses in every retained space: vacuum-to-one-particle source,
one-particle diagonal energy, and one-step one-particle hop. Virtual
insertions of that coefficient start later. It cannot be scalar on the
retained space. Thus its coefficients must vanish. The parent itself
has `2D2` at order two and `4D4` on the order-two kernel at order four.
The order-five change from alpha is supported on D4 and vanishes after
the third projection.

For a directly implementable bond, (3)--(5) expand to

```
u = 2 epsilon²-34 epsilon⁴-16 Re(alpha) epsilon⁵
    +[196-2|alpha|²-16 Re(beta)-(128/3)mu] epsilon⁶,
t = -epsilon²+8 epsilon⁴+2 Re(alpha) epsilon⁵
    +[-80+2 Re(beta)+(128/3)zeta] epsilon⁶,
f = (128/3)eta epsilon⁶.                                    (6)
```

Unlisted higher coefficients can be zero; this defines a literal
polynomial microscopic family with no L-dependent parameter. The exact
root-of-unity curve has alpha=0 and beta=-24 in this chart. Holding
`p=3 epsilon²` is not the exact resonance curve at finite epsilon, but
beta does not change the leading operator below. We do not confuse the
two curves.

## 4. Complete leading many-body map

Let P exclude all physical separations 1,2,4, and put

```
t3 = 128/3,                 a = alpha*/8,
G8 = sum n_j n_(j+8)(1-n_(j+3))(1-n_(j+5)),
v_a = (0,-1,1,a),           h_a = |v_a><v_a|,
```

where `h_a(j,j+3)` acts on the two indicated sites. The **entire** first
retained operator, including every particle-number sector, is

```
K6/t3 = H(a,zeta,eta,mu)
      = P[sum_j h_a(j,j+3)+(3/16)G8
          +zeta T1+zeta* T1†+eta B+eta* B†-mu N]P.           (7)
```

In particular

```
sum h_a = -(T3+T3†)+(-2+|a|²)D3+2N
          +sum[a(n_(j-3)-n_(j+3))b†_j+h.c.].               (8)
```

This is a seven-real-parameter leading family, with an additional
microscopic epsilon setting the hierarchy of scales. Beta and the fifth
and sixth pair-amplitude jets drop out of (7), after the displayed
counterterms. Five of these seven leading controls are the actual normal
complement (4); a is the tangent-parent displacement that changes the
resonance.

For a=0, (7) follows from the complete undeformed sixth-order calculation.
For general a the local identity is the previously derived
[alpha companion formula](../qubit_resonant_liquids_20261116/DERIVATION.md).
The new content of (3)--(7) is that these terms **exhaust the entire
physical deformation chart under the stated lower-kernel conditions**,
not just a chosen projector path. The direct normal controls in (4)
enter at order six, so virtual corrections involving them start later.

`check.py` constructs literal microscopic bonds for two choices with all
seven controls and arbitrary complex beta, pair5 and pair6, on L=10,12,13.
It eliminates every occupation sector three times and compares every
retained matrix element with a separate occupation-space implementation
of (7). This reuses the antecedent's Schur recursion, so it is not an
independent verification of that recursion. The whole-matrix comparison
and source hashes are saved in `checks.json`.

### Number conservation and low-momentum cancellations

Equation (7) commutes with N **if and only if a=eta=0**, on sufficiently
long rings. Vacuum-to-one creation forces eta=0. The matrix element from
one particle at 0 to particles at 0,3 then forces a=0. There is no
independent repulsion-only direction: its coefficient is |a|² and its
charge-changing companion is a, exactly as in (8).

One must also not overstate this result as immediate relevant locking.
Before the P compression, the charge-increasing cubic term has the
symmetrized Fourier form

```
(-i a/sqrt(L)) sum_(p,k) [sin(3p)+sin(3k)]
                                  b†_p b†_k b_(p+k).        (9)
```

This identity uses `b_j=L^-1/2 sum exp(ikj)b_k`. The coefficient vanishes
at all three root-of-unity centers; indeed it vanishes identically for
a one-particle input with total momentum 0,+/-2pi/3. Equivalently,
creation of a distance-three pair from a one-particle plane wave carries
`a(1-exp(3iK))`, up to a common basis phase. Thus the three one-particle
W states remain exact eigenvectors when eta=0. Derivative vertices and
higher processes still matter. A nonzero local commutator is not by
itself proof of a relevant infrared coupling.

For real a,zeta and eta=0, inversion times product Z is a symmetry.
It permits the opposite-valley pair and the zero-valley pair; at the
exact centers it also permits relative pair transfers such as
`(b†_0)^2 b_+ b_-`. The latter remains momentum matched when the two
nonzero centers are opposite. Their actual matched amplitudes and
dimensions must be calculated, not removed by imposing a fictitious
microscopic U(1). This contribution does not assign those dimensions.

## 5. What the exhaustion does not prove

The formal leading Hamiltonian is not a uniform-in-L convergence theorem
for a fixed microscopic epsilon. A general a,zeta,eta,mu phase diagram is
not solved. In particular the number-changing finite-density class is
not excluded from having extra neutral modes. Nonanalytic approaches,
different lower constraints and other microscopic neighborhoods are
outside (5). Even an exclusion of this entire chart would not give a
universal qubit bound on unitary compact Virasoro c.
