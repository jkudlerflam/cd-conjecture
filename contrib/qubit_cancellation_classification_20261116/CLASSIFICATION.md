# Normal form and complete second-order cancellation locus

Throughout, `b†=|1><0|`, `n=b†b`, `D_r=sum n_j n_(j+r)`, `N=sum n_j`,
`T_r=sum b†_(j+r)b_j`, `B=sum b†_j`, and `A=sum b†_j b†_(j+1)`.
Hermitian conjugates of complex-amplitude terms are always included.
We use rings of sufficient length to distinguish local words (L>=7 at
second order); the displayed identities are bulk identities.

## 1. Complete microscopic coordinates [proved]

Modulo constants and telescoping bonds `K tensor I - I tensor K`, every
uniform nearest-neighbour Hermitian qubit Hamiltonian has the unique form

```
H = v D1 + u N + [t T1 + p A + f B
                 + sum_j (a n_(j-1)+b n_(j+1)) b†_j + h.c.].       (1)
```

Here v,u are real and t,p,f,a,b are complex: twelve real parameters.
These span the nine exchange Pauli sums and three uniform fields. The
letter b without a site index in (1) denotes the right-facilitation
coefficient, not the annihilation operator.

For `H=H0+epsilon V1+...`, let P be the **complete** no-`11` subspace.
Then `P V1 P` is scalar iff

```
u1=t1=f1=0.                                                   (2)
```

The scalar is zero after removing constants. Necessity follows by
comparing vacuum and one-particle diagonal elements, the vacuum-to-one
transition, and the one-step one-particle hop. The v,a,b,p terms cannot
have any of those P-to-P elements. Sufficiency follows because every
facilitated creation or pair creation leaves P, while their annihilations
vanish on P. The real kernel dimension is 7.

Write `a=s+g`, `b=s-g`. The complex s direction is precisely the tangent
to common onsite rotations of `D1`; positive scale removes v1. The
stabilizer z rotation changes the phases of g and p with charges one and
two, so their relative phase must be retained. One cannot make both real
in general.

This gauge can be imposed on an entire sufficiently small neighborhood:
choose the isolated eigenaxis near z of the symmetric exchange matrix.
Its two transverse symmetric entries then vanish, giving `a+b=0`.
Divide by the positive D1 coefficient. When g is nonzero, use the residual
common z rotation to make it positive real. These are common onsite
unitaries, with no staggered or site-dependent transformation. On an
analytic arc with nonzero leading g, reparameterize by `g=epsilon`.
Arcs with vanishing leading g require a separate valuation/reindexing;
the unreduced coordinates (1) remain complete there.

## 2. The full second-order operator [proved]

Keep general complex a,b,p in V1, and write u,t,f for the visible
coefficients in V2. Set `R=Q H0^-1 Q`. Then

```
K2 = P V2 P - P V1 R V1 P
   = P { -L |p|²
       + [u+3|p|²-|a|²-|b|²] N
       + [|a-b|²/2-|p|²] D2 - |p|² D3/3
       + [(t-a b*) T1 - |p|² T2/2 + h.c.]
       + sum_j { b†_j [f-p(a+b)*
                    + p(a-b)*(n_(j+2)-n_(j-2))/2] + h.c. }
       } P.                                                   (3)
```

The coefficient convention in the last line is `p times conjugate(a-b)`.
The v1 direction does not contribute to K2, since `D1 P=0`.
Invisible coefficients in V2 also have zero P compression.

Local proof. A virtual pair on an empty bond costs
`1+n_(j-1)+n_(j+2)`. Summing its diagonal contributions on P gives
`-|p|²[L-3N+D2+D3/3]`. Overlapping pair operations give the distance-two
hop `-|p|²/2`. Facilitated creation at an empty site with m occupied
neighbours costs m. The diagonal sum is
`-(|a|²+|b|²)N+|a-b|² D2/2`; exchanging the virtual pair gives `-a b* T1`.
For creation at j, the two mixed paths have amplitudes

```
-p (a+b n_(j+2))* /(1+n_(j+2)),
-p (b+a n_(j-2))* /(1+n_(j-2)).
```

Their sum is the last line of (3). Two raising operations cannot return
to P with Delta N=2. Disjoint operations leave a forbidden pair. This
exhausts the two-step paths. In particular there is no missing correlated
or number-changing term beyond the explicitly conditioned source.

## 3. Cancellation and all classical coexistence points [proved]

The distance-two one-particle matrix element forces p=0 if K2 is to be
diagonal. No nearest-neighbour counterterm supplies that matrix element.
The remaining off-diagonal elements force

```
t=a b*,   f=0,
K2=P(delta N+k D2)P,
delta=u-|a|²-|b|²,   k=|a-b|²/2 >=0.                          (4)
```

K2 is scalar on all of P iff `k=delta=0`. Then a=b is purely rotational:
after removing the common rotation there is no nonzero physical V1.
This case restarts the analysis at the next nontrivial order, rather than
creating a new first-order cancellation mechanism.

For k>0, (4) has two exponentially degenerate coexistence points.

1. **delta=0:** the ground space P2 forbids separations one and two.
2. **delta=-3k:** define `Z3=sum (1-n_j)(1-n_(j+1))(1-n_(j+2))`.
   On P, `Z3=L-3N+D2`, so `K2+kL=k Z3`. Its kernel Pd forbids `11`
   and `000`. Its particle spacings are precisely two or three.

These are complete. For `-3k<delta<0`, write
`delta N+kD2=(delta+3k)N+kZ3-kL`: its minimizing bulk configuration
has spacing three. For delta<-3k the spacing-two crystal minimizes;
for delta>0 the vacuum does. Equivalently minimize the energy per site
over runs between occupied sites: gaps longer than three do not improve
either competing coexistence energy. The limiting densities are 0,1/3,1/2
with crossings at 0 and -3k. Rings incompatible with the crystal period
can contain a bounded number of defects; these do not create an
exponentially large bulk ground manifold away from coexistence.

Other chemical potentials can of course be shifted to a positive
classical Hamiltonian. They select the stated vacuum/crystals rather
than a further extensive constrained space. Positivity without a fixed
choice of energy zero is not an additional cancellation condition.

In the rotation and scale gauge `a=1,b=-1`, the two coexistence paths are

```
V1=F=sum (n_(j-1)-n_(j+1))(b†_j+b_j),
V2=2N-T1-T1†+q A+q* A†          (dilute),
V2=-4N-T1-T1†+q A+q* A†         (dense).                     (5)
```

The complex q is the order-epsilon² pair amplitude. It is not a leading
pair amplitude p1. This distinction explains why (5) survives (3).
