# Sixth order on the further dilute kernels

This extends the complete fourth-order classification; it is **not** a
classification of every complex-q sixth-order interaction or every higher
microscopic jet. The four slices below have complete leading operators,
including all particle-number sectors. The generic complex case is also
implemented by the exact Schur recursion, with small saved diagnostics.

Take the canonical microscopic path

```
H=D1+epsilon F+epsilon²[2N-T1-T1†+q A+q* A†]
  +epsilon⁴[-(4+4 Re q+2|q|²)N+2(1+Re q)(T1+T1†)],           (23)
```

with all unlisted jets zero. Add `L|q|² epsilon⁴` to remove its fourth-order
scalar. Eliminate Q, then P minus P2, then the positive fourth-order
constraint. In the four cases below the entire fifth-order compression
vanishes. The resulting sixth coefficient is the first nonzero operator.

**[Derived finite-order formulas; exact full-matrix verification at L=10;
independent numerical full-matrix checks through L=16]** With projectors
on both sides of each displayed operator:

| q | Further space | Sixth-order operator |
|---|---|---|
| 0 | P4 | `-16N+8(T1+T1†)+(25/3)D5+(25/12)D6+(7/12)D7+D8/8` |
| -1 | P4 | `2L-4N+(256/3)D5+(128/3)D6+(256/11)D7+8D8` |
| 1 | P3 | `2L-50N+24(T1+T1†)-(T4+T4†)-2D4+(10/3)D5+(14/3)D6` |
| 3 | no separations 1,2,4 | `18L-(332/3)N+80(T1+T1†)-(128/3)(T3+T3†)-(256/3)D3+8G8` |

Here

```
G8=sum_j n_j n_(j+8)(1-n_(j+3))(1-n_(j+5))                 (24)
```

counts consecutive occupied particles separated by eight on the last
space. Omitting its conditioning would give incorrect three-particle
diagonal entries. The L=15,16 checks contain three- and four-particle
sectors, not just isolated particles and pairs. Finite rings too short
to separate distances can alias D5,...,D8; the implementation retains the
literal periodic sums, including their multiplicities.

The formulas follow by applying (8) a third time with the fourth-order
diagonal operator as the unperturbed denominator. Exact rational checks
compare the entire matrices, and direct diagonalization of the literal
spin chains checks the leading epsilon^6 scaling for q=1,3. At the
coarsest q=3 point epsilon=.075, the band gap is small (about .07 epsilon^4
on L=10); it grows toward 4 epsilon^4 as epsilon decreases. This is a
finite-L check, not a uniform thermodynamic gap assertion.

## What these slices show

At q=0 and -1 the available sixth-order chemical and nearest-neighbour
hopping counterterms can cancel the displayed motion, leaving positive
penalties at **every separation five through eight**. Thus the next
constraint is P8. These counterterms agree with the sixth-order jets of
the exact rank-one parents (u6=16,t6=-8 at q=0; u6=4,t6=0 at q=-1).
Even the known low-entanglement parents therefore support multiple rounds
of growing exclusion constraints in their perturbative description.
The existence of this hierarchy is not evidence of increasing conformal c.
An all-orders conclusion about its depth is not made.

At q=1, after cancelling the nearest-neighbour hop with t6=-24, a
distance-four hop of coefficient -1 survives. At q=3 the corresponding
distance-three coefficient is -128/3 after t6=-80. A nearest-neighbour
V6 cannot cancel either matrix element. Higher pair/source jets that
preserve the lower scalar/dilute constraints change chemical terms,
nearest-neighbour terms, or derivatives of the fourth-order diagonal
penalties; they do not supply these one-particle matrix elements at this
order. These are exceptions to any argument based on a permanently
degree-two effective one-particle band.

The one-particle bands with the nearest-neighbour coefficient cancelled
are proportional to `-cos(4k)` and `-cos(3k)`, with four and three minima.
This is a **band statement**, not a c=4 or c=3 construction. The displayed
many-body operators have constraints and interactions; extra conservation
laws at leading order can have zero-velocity sectors, and later terms can
couple or gap them. Neither a stable multicomponent liquid nor the absence
of additional neutral sectors follows from the band.

These two exceptional values coincide with limiting transfer-matrix
resonances of the parent family: the two eigenvalues are, to leading order,
`epsilon[-1 +/- sqrt(-q)]`. Their ratios at q=1 and q=3 are fourth and
third roots of unity. This explains why extrapolating the generic
rank-two-ground-state argument through these points is unsafe. It does
not make a detuned Hamiltonian an exact parent.

For arbitrary complex q away from 1 and 3, `dilute_higher` supplies the
complete formal sixth matrix on P4. The two saved complex examples have
vanishing fifth coefficient and no number-changing sixth entries, but
no exhaustive closed form or positivity theorem for all such q is
asserted. The **complete symbolic classification claimed by this
contribution remains the second-order locus and the two fourth-order
operators**, with the stated dense scalar obstruction and these further
dilute slices. More detunings and subsequent constraints remain open.
