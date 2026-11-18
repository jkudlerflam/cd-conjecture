# The literal five-control extension, still physical d=9

## Alphabet, grammar and local Hamiltonian

The physical states are `L0,L1,R0,R1,DA,DB,E,MA,MB`, with charges
`(0,2,0,2,-1,1,0,-1,1)`. Markers are E,MA,MB. R0 doubles as the
neutral head; there is no uncounted vacancy or boundary state.

The allowed ordered pairs are inherited unchanged:

| Left state | Allowed right states |
|---|---|
| L0 or E,MA,MB | Every state except R1 |
| L1 | L0,L1,R0,DA,DB |
| R0,R1,DA,DB | R0,R1,E |

Pi_bad is the sum of the other 32 ordered-pair projectors. The
off-diagonal local rules, each with its Hermitian conjugate, are:

| Input | Output | Coefficient |
|---|---|---:|
| R0 R_a | L_a R0 | -1, a=0,1 |
| R0 E | L0 E | -1 |
| D_s R_a | L_a D_s | -z_s(a), a=0,1; s=A,B |
| D_s E | L0 M_s | -1, s=A,B |
| DA R1 | DB R0 | g |
| E R0 | MA DB | -lambda |
| E R0 | MB DA | +lambda |
| L0 L1 | L1 L0 | t0 kappa_L |
| R0 R1 | R1 R0 | t0 kappa_R |

Here z_A(a)=(-1)^a and z_B(a)=1. Define

```
N=P_L1+P_R1,  N_L=P_L1,
P_exc=P_DA+P_DB+P_MA+P_MB,  P_M=P_E+P_MA+P_MB,
o=Delta P_exc+mu_sel P_M+t0(mu N+eta N_L),
h=h_off+K Pi_bad+(o tensor I+I tensor o)/2
      +t0 b P_M tensor P_L1.                         (1)
```

The last term is an actual nearest-neighbour boundary density, not an
addressed logical field. The same h is used on every physical bond.
Every term conserves total charge and both sites' marker projectors
on the full 81-dimensional bond space. Frozen marker positions therefore
remain exact, even on illegal words. `physical.py` checks these identities
and implements (1).

There remain 14 Hermitian transition pairs, or 28 directed off-diagonal
matrix units. Independently weighting LL and RR adds a parameter, not
another transition. The boundary density adds three diagonal matrix units,
one for each marker. Grammar and onsite diagonal projectors are additional.
There are eleven grouped real coefficients including the common walk
scale, or ten after fixing that scale. Matching coefficients does not
reduce the literal operator count.

## Five-site cell and physical projections

For three payload bits, the legal cell has four nonmarker positions
and one marker. Head positions j=0,1,2,3 have words

```
L_(a_0)...L_(a_(j-1)) T R_(a_j)...R_(a_2) E,
T in {R0,DA,DB};
```

the three terminal words are `L_(a_0)L_(a_1)L_(a_2)L0 m`, with
m=E,MA,MB. There are 120 legal words, and eight neutral code states.
The terminal L0 is retained. The neutral path has five vertices, with

```
E_v=-sqrt(3),
v_j=sqrt(1/3) sin(pi(j+1)/6), j=0,...,4,
v_j²=(1/12,1/4,1/3,1/4,1/12), alpha²=1/12.
```

Set Delta=2-sqrt(3) at this intended cell size. The gauged charged
carrier above the vacuum is `H_car=2I-adjacency(path_5)`. Its conversion
vertices are j=0,1,2; the last nonmarker and terminal vertices remain
inactive. Put `t0=lambda² alpha⁴ g²=lambda² g²/144`.

For arbitrary n the five projections are

```
M_n=mu I+eta diag(p_i)+b p_0 |0><0|+kappa_L B_L+kappa_R B_R,
p_i=sum_(j>i) v_j²,
(B_L)_(i,i+1)=sum_(j>=i+2) v_j²,
(B_R)_(i,i+1)=sum_(j<=i) v_j².                         (2)
```

At n=3 this becomes

```
p=(11/12,2/3,1/3),
B_L adjacent entries=(2/3,1/3),
B_R adjacent entries=(1/12,1/3).
```

The three diagonal profiles are independent. The determinant of the
two hopping profiles is 7/36, so (2) spans every real tridiagonal 3 by 3
matrix. It supplies no (0,2) entry. All calculations stay in this space.

`verify_projections` assembles the LL and RR matrices on actual legal
words and projects on the entire eight-dimensional code. It projects
the literal boundary density on the entire 64-dimensional pair code,
including the left cell's physical marker. The results agree with the
fermionic quadratic operators from (2), not merely with one-particle
blocks. The largest discrepancy is 4.5e-16. Swaps involving the shared
head or fixed terminal bit can leave the legal code; their finite-K
effects remain a separate issue.

## The selected rational coefficients

The five coefficients in (2) are

```
(mu,eta,b,kappa_L,kappa_R)=(-3385,9387,-7486,-2560,11999)/500000,

M=[ -19709/6000000, -2827/2000000,             0
     -2827/2000000,   2873/500000, 9439/1500000
                0, 9439/1500000,      -8/15625 ].      (3)
```

These rational numbers define the point, not rounded printed eigenvalues.
After combining the overlapping onsite fields, the five added physical
coefficients on L1, R1, boundary density, LL swap and RR swap are
`.012004,-.006770,-.014972,-.005120,.023998`, respectively, with the
common t0 removed. Their absolute-value spread is 4.687109375. This is
only the added weak-control subset; a full microscopic spread cannot be
reported without choosing and certifying g,lambda,K and the selector.
The order-one walk, large penalty and very small t0 remain a substantial
hierarchy. The narrow independent control box in BANDS.md is another
simplicity cost, separate from the small emergent velocities.
