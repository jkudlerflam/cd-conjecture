# Exact cubic determinant and six simple roots

## Carrier coefficient and admissible onsite matrix

All matrices here use t0 as the quadratic energy unit. At n=3 set

```
R=(H_car tensor I+I tensor H_car)^(-1), x=R|4,0>,
U_i=(|i><i| tensor I)x, V_i=(I tensor |i><i|)x, i=0,1,2,
T=-2 U^T R V.
```

The combined generated onsite coefficient cancels at this order,
as checked on the full pair Fock space. The only added onsite M is
the literal five-control projection in BOND.md. Numerically,

```
T=[ -.000413852, -.000526713, -.000381879
    -.002120031, -.002472602, -.001657905
    -.007234894, -.007009032, -.004048489 ].
```

The defining rational entries, obtained by an exact 25 by 25 inverse,
are saved in `exact.json`. `exact.py` verifies that inverse by exact
multiplication. No fitted hopping matrix is used.

## Cubic certificate

Let `S=(T+T^T)/2`, `A=(T-T^T)/2`, and
`a=(A_12,-A_02,A_01)`. For x=cos(k), the reciprocal symbol has

```
h(k)=M+2xS+2i sin(k) A,
p(x)=det h(k)
    =det(M+2xS)-4(1-x²) a^T(M+2xS)a.                 (1)
```

Equation (1) is a real cubic with leading coefficient `8 det(T)`.
This follows directly by expanding a 3 by 3 determinant; no eigenvector
congruence or nonorthonormal basis is used. The entire polynomial has
rational coefficients at (3) of BOND.md and is saved exactly.

**[Proved finite algebra.]** Its signs at the following rational endpoints
are certified exactly:

| Interval for x | Sign at left | Sign at right |
|---|---:|---:|
| (-7/10,-13/20) | - | + |
| (9/50,23/100) | + | - |
| (9/20,1/2) | - | + |

A nonzero cubic with these three disjoint sign changes has exactly three
distinct, simple roots in (-1,1). Hence the symbol has exactly six simple
Fermi roots on the momentum circle. A multiple zero eigenvalue would
make the determinant derivative vanish; it is therefore excluded here.

Exact interval arithmetic preserves all six endpoint signs when each
of the five physical coefficients varies independently by at most 10^-9.
This certifies six roots throughout that box. It does not certify a broad
parameter region, and it is not a box certificate for the interacting flow.
The 32 corner flow evaluations reported elsewhere are numerical samples.

## Velocities and basis

The three cosine roots are approximately
`(-.671731949099,.206477694842,.469499752875)`.
The sorted roots and actual eigenvalue slopes are:

| Root label | k | slope / t0 |
|---:|---:|---:|
| 0 | -2.30734060562 | -7.42972037e-7 |
| 1 | -1.36282262103 | +3.93151597e-7 |
| 2 | -1.08207220871 | -5.71505210e-7 |
| 3 | +1.08207220871 | +5.71505210e-7 |
| 4 | +1.36282262103 | -3.93151597e-7 |
| 5 | +2.30734060562 | +7.42972037e-7 |

The Dirac channels (R,L) are (5,0), (1,4), (3,2). All velocities and
vertices use normalized eigenvectors of the original Hermitian h(k).
The velocity ratio is about 1.89; their very small absolute scale matters
for matching errors and weak coupling.

There is also an exact conservative velocity bound. On each rational
root bracket, interval arithmetic gives a nonzero lower bound m for
|p'(x)|. Let

```
H_bound=||M||_infinity+||T||_1+||T||_infinity.
```

Then `||h||_2<=H_bound`. At a simple root with other eigenvalues E1,E2,
`|v|=sqrt(1-x²)|p'(x)|/|E1 E2|`. Using sqrt(1-x²)>1/2 gives
`|v|>=min(m)/(2 H_bound²)>1.20214e-8` at the primary point.
This bound is weaker than the numerical velocities but strictly positive.

## Momentum channels and search scope

All 21 unordered root pairs with replacement are compared modulo 2pi.
The only collisions found are among (0,5), (1,4), (2,3), all of sum zero.
The minimum other chord separation is .11256365. Thus the primary point
retains three Cooper-pair transfers; no additional marginal umklapp
appears in the complete finite momentum audit. Repeated-fermion local
quartics vanish without derivatives. Six-fermion operators have free
scaling dimension three and are irrelevant sufficiently near zero coupling.

The optional search solved three determinant-coefficient constraints
within five tridiagonal variables, using twelve initial guesses per
root design and fixed evaluation budgets. All 108 solves gave six-root
symbols. Forty-eight had extra or nearly resonant momentum channels and
were excluded from the generic flow equations; the remaining sixty all
grew under those equations. A finite denominator ladder then selected
the rational primary point. None of these searches proves that stable
six-root points are impossible elsewhere in the permitted control space.
