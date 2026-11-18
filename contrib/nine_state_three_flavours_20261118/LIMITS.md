# Error budgets and the exact extent of the result

## The small conducting scale cannot be ignored

The primary point has `v_min/t0=3.93152e-7`. The positive ratio of its
three velocities does not make their absolute scale large. The initial
dimensionless forward and transfer coefficients are of order hundreds
times g²; weak coupling requires those products to be small. The physical
walk, conversion, port, penalty and five controls retain separate scales.

At fixed n=3 the charged gap is `delta=2-sqrt(3)`. The finite-cell
Neumann estimate inherited from the actual resolvent is

```
||K(g)-K0-g²K2-g⁴K4||
 <= (|g|/delta)^6/[delta(1-(g/delta)²)], |g|<delta.     (1)
```

After restoring lambda² alpha⁴ and dividing by t0, the bound is
`epsilon=g⁴/[delta^7(1-(g/delta)²)]`. It must be compared both with
v_min and with the leading interaction, which is another factor g²
smaller than t0.

For the diagnostic value g=10^-6, (1) gives

| Quantity | Size in t0 units |
|---|---:|
| Smallest free velocity | 3.93152e-7 |
| Norm of g² K4 | 1.11956e-13 |
| Conversion tail bound | 1.00840e-20 |

This is a finite conversion diagnostic, not a certified finite tuple
for the complete microscopic model. The script also recomputes roots,
velocities and leading vertices with the explicit K4 quadratic part
included. Its changes are small and saved in `audit.json`; this is not
a self-consistent interacting Fermi-surface calculation.

A useful fixed-band vertex estimate is explicit. Compression to the
vacuum, one- and two-particle blocks and normal ordering give
`||delta W||<=6 epsilon`. Each two-cell Bloch wedge has norm at most 2,
so `|delta Gamma|<=24 epsilon`. After dividing by g² and inserting the
velocity denominators,

```
|delta(D_ij/g²)| <= 48 epsilon/(pi v_min g²).          (2)
```

At g=10^-6 this is .39190, below one tenth of the primary negative-D
margin 7.45124. The same bound gives transfer errors far smaller than
the nonzero bare transfers. This controls the conversion remainder's
constant/quadratic/quartic matching at fixed free band data. At order
g^6, higher fermion operators also occur; (2) is not a replacement for
their filled-sea matching or a volume-uniform RG bound. Smaller g makes
the relative conversion corrections smaller.

Other local errors E need the same small-scale accounting. For example,
the conservative fixed-band target of less than one tenth of that D
margin is

```
E/(t0 g²) < pi v_min * 7.45124 / 480 ~= 1.92e-8.      (3)
```

Here v_min denotes its t0-removed value. Unknown finite-K and higher-port
local bounds must fit an allocated part of this budget, and their
quadratic pieces must also preserve the six crossings. Merely requiring
these errors to be small relative to the walk scale is insufficient.
If bounds of the forms C_K/K and C_port lambda^4 are established, they
must be compared with the right side of (3) times t0 g². The constants
and the necessary volume-uniform bounds have not been established here.
An asymptotic choice lambda/g² -> 0 followed by sufficiently large K can
separate the formal scales at fixed n, but it is not a finite-member
microscopic certificate.

## Assumptions that remain unresolved

1. The band and coefficient calculations use the legal K=infinity code.
   The literal walk and background swaps may leave it at finite K.
   Their induced splittings require an actual local error bound.
2. Marker positions are exactly frozen. No phonon, positional mode or
   mobile-vacancy contribution to c can be added.
3. A zero-port vacuum selector can start from the midpoint choice for
   cell length five: `mu_sel=2-(d_5+d_4)/2`, with
   `a_b=2[1-cos(pi/(b+1))]`, `d_b=(b+1)a_b-b a_(b+1)`.
   This selects the intended legal vacuum length only. With Delta fixed
   at 2-sqrt(3), the nonzero-port interacting ground energy still needs
   comparison against charged, wrong-length, mixed and illegal sectors.
4. A finite two-cell resolvent calculation is not a uniform reduction
   of a gapless infinite chain. Multi-bond port processes and extensive
   error bounds are separate work.
5. The invariant-cone proof concerns the leading one-loop equations
   and numerically matched initial data. It establishes a perturbative
   instability under the usual continuum matching assumptions, not a
   microscopic gap theorem. It does not determine the final Virasoro c.

## What generalizes, and what is special to n=3

The physical alphabet, frozen markers, carrier resolvent, literal five
control profiles and full-Fock extraction of the generated interactions
exist at every finite n. At generic reciprocal momenta, the Cooper
matrix and forward-scattering equations also have the corresponding
intermediate-channel sums. These are structural continuations, not
proofs of more Fermi points or stable liquids.

At n=3, three independent diagonal profiles plus two independent hopping
profiles span all real tridiagonal onsite matrices. This is special:
the space has dimension 2n-1, while the five-control family still has
dimension at most five for n>3. It never supplies a general symmetric M.
The cubic sign-change argument is also a finite n=3 certificate, not an
all-n root construction. The carrier gap shrinks with n, so the same
small-parameter radius cannot be assumed uniform.

The bounded search found no stable generic-momentum three-flavour seed.
That observation and the strict instability of the selected point leave
the full five-control family undecided. The points excluded for additional
or near umklapp would need their extra channels included before any flow
conclusion. No universal d or c bound follows, and no connection to the
unlocated reported d=10 arbitrary-c GitHub example is established.

The useful result is therefore specific: literal band synthesis reaches
six roots without adding states, but the selected generated interaction
destabilizes the three-channel Gaussian theory through coupled transfer.
This is a scaling test of the earlier stable two-flavour seed, not its
extension to an arbitrary-c construction.
