# A literal nine-state bond and its exact code projections

## Physical alphabet and complete local rules

Use matrix-unit notation `|output><input|+h.c.`. The nine orthogonal
physical states are `L0,L1,R0,R1,DA,DB,E,MA,MB`.
Charges are `q(L_a)=q(R_a)=2a`, `q(E)=0`,
`q(DA)=q(MA)=-1`, and `q(DB)=q(MB)=1`. Markers are E,MA,MB.
There is no vacancy or uncounted head state: R0 is also the neutral
head. This is physical d=9, not an encoded or blocked dimension.

The grammar's allowed ordered pairs are:

| Left state | Allowed right states |
|---|---|
| L0 or any marker | Every state except R1 |
| L1 | L0,L1,R0,DA,DB |
| R0,R1,DA,DB | R0,R1,E |

Let Pi_bad be the sum of the other 32 pair projectors. The entire
off-diagonal bond is the following table plus Hermitian conjugates:

| Input | Output | Coefficient |
|---|---|---:|
| R0 R_a | L_a R0 | -1, a=0,1 |
| R0 E | L0 E | -1 |
| D_s R_a | L_a D_s | -z_s(a), a=0,1; s=A,B |
| D_s E | L0 M_s | -1, s=A,B |
| DA R1 | DB R0 | g |
| E R0 | MA DB | -lambda |
| E R0 | MB DA | +lambda |
| L0 L1 | L1 L0 | t0 kappa |
| R0 R1 | R1 R0 | t0 kappa |

Here `z_A(a)=(-1)^a`, `z_B(a)=1`. Set

```
N=|L1><L1|+|R1><R1|,  N_L=|L1><L1|,
P_exc=sum_(DA,DB,MA,MB) |s><s|,
P_M=sum_(E,MA,MB) |s><s|,
o=Delta P_exc + mu_sel P_M + t0(mu N+eta N_L),
h=h_offdiag + K Pi_bad + (o tensor I+I tensor o)/2.       (1)
```

Equation (1) defines a one-site-translation-invariant Hermitian
nearest-neighbour Hamiltonian for every fixed set of its coefficients.
Every rule preserves total additive charge and each physical site's
marker projector, including on illegal words. It is therefore **frozen**.
`literal.py` implements these matrix units independently of the scalar
resolvent. No interaction is supplied by hand in logical space.

## Legal code and parameter convention

A cell has p=n+1 nonmarker sites and one marker, total b=n+2. For
n bits a_0,...,a_(n-1), its legal words with head position j=0,...,n are

```
L_(a_0)...L_(a_(j-1)) T R_(a_j)...R_(a_(n-1)) E,
T in {R0,DA,DB}.
```

The three terminal words are `L_(a_0)...L_(a_(n-1)) L0 m`,
with m=E,MA,MB. Thus the legal cell dimension is `3(n+2)2^n`,
and its neutral vacuum code has dimension `2^n`. The fixed L0 at the
terminal costs one flavour; omitting it would miscount the code.

Put s=n+2. At g=lambda=t0=0, the neutral clock is a length-s path,

```
E_v=-2 cos(pi/(s+1)),
v_j=sqrt(2/(s+1)) sin(pi(j+1)/(s+1)), j=0,...,s-1,
G|a>=sum_j v_j |j;a>,  alpha=v_0=v_(s-1).              (2)
```

The charged A path is gauged by the parity of crossed payload bits;
B needs no such gauge. Choose `Delta=2+E_v` at the intended n.
Their common positive scalar Hamiltonian above the selected vacuum is

```
H_car=2I-adjacency(path_s),   delta_car=2-2cos(pi/(s+1)).
```

Only positions 0,...,n-1 convert. Position n and the terminal marker
remain inactive carrier vertices; dropping either changes the answer.
Conversion is `|B><A| f_i + |A><B| f_i^dagger` at position i.
In global fermions the port source, after removing alpha², is
`|B,A>-|A,B> P_left`, with position source `(s-1,0)`.

The primary seed has n=2, s=b=4, `Delta=2-2cos(pi/5)` and
`alpha=sqrt(2/5)sin(pi/5)`. Define `t0=lambda² alpha⁴ g²`.
All scalar coefficients in the subsequent files use t0 as the
quadratic energy unit. The microscopic coefficients remain fixed when
comparing other marker sectors; Delta is not retuned cell by cell.

The vacuum-clock length selector can start from the inherited legal,
zero-port choice `mu_sel=2-(d_4+d_3)/2`, where
`a_b=2[1-cos(pi/(b+1))]`, `d_b=(b+1)a_b-b a_(b+1)`.
It selects b=4 among those vacuum cells. This is not a nonzero-port
or finite-K ground-sector certificate; see LIMITS.md.

## The added physical terms really give the required onsite matrix

Projection of the total physical N is exactly logical total number.
Projection of N_L is diagonal, with coefficient

```
p_i=sum_(j>i) v_j².
```

The two background swaps together project to adjacent logical hopping

```
B_(i,i+1)=B_(i+1,i)=1-v_(i+1)².
```

For the first formula a bit is L precisely when the clock has passed
it. For the second, adjacent payload bits are physically contiguous
except when the head lies between them. Swaps involving the shared
R0 head or the fixed terminal bit leave the legal code and vanish in
the leading projection; their finite-K effects are not silently discarded
from a claimed finite-K theorem.

Consequently the added leading onsite matrix is

```
M=mu I+eta diag(p_i)+kappa B.                           (3)
```

For n=2, `p=(.861803398875,.5)` and `B_01=.638196601125`.
These three physical controls span every real symmetric 2 by 2 M.
The full Fock-space projections are independently checked, at n=2
and 3, rather than inferred from one-particle matrix elements alone.

There are nine matched Hermitian walk transitions, one conversion,
two port transitions and two background swaps: 28 directed off-diagonal
matrix units in the displayed basis. The grammar and onsite diagonal
projectors are additional operators. Matching several coefficients does
not reduce this literal count. Allowing the listed grouped coefficients
to vary gives nine real couplings including the common walk scale,
eight after scale fixing; the seed constrains several of them by the
explicit formulas above. A complete finite-member coefficient spread
is not supplied because g,lambda,K and the dressed selector have not
been certified at a chosen finite tuple. The required scale separation
is a material simplicity cost, not merely an emergent small velocity.

Even inside the three added weak controls, the physical onsite
coefficients are `mu+eta` on L1 and `mu` on R1. Their common t0-removed
values, together with the two equal swap coefficients, are approximately
`.00236015,.01851689,.00791520,.00791520`, with spread about 7.85.
Reporting only the grouped mu,eta,kappa ratios would hide this difference.
No full microscopic spread is inferred from this restricted comparison.
