# Exact four-root symbol, with reciprocal pair resonance retained

All matrices below use the energy unit `t0=lambda² alpha⁴ g²`.
Constants and telescoping onsite contributions have been removed.

## Coefficient from the common carrier

For n=2 the carrier has four vertices. Set

```
H_car=2I-adjacency(path_4),
R=(H_car tensor I+I tensor H_car)^(-1), x=R|3,0>,
U_i=(|i><i| tensor I)x, V_i=(I tensor |i><i|)x, i=0,1,
T_ij=-2 U_i^T R V_j.
```

The exact rational result is

```
T = [ -29114923/28749600000,  -9838709/9583200000
      -53812309/9583200000,  -69800369/14374800000 ].       (1)
```

Its Hermitian part S=(T+T^T)/2 is indefinite. For A=(T-T^T)/2,

```
det S = -1010260494073711/165307900032000000000,
det T = -140105620441711/165307900032000000000,
rho²=A_01²/(-det S)=870154873632000/1010260494073711 < 1.   (2)
```

There is no appeal to the failed positive-Hermitian-hopping test.
`exact_symbol.py` certifies the inverse identity and rational signs.

## Explicit onsite matrix and root identity

Since S has one eigenvalue of each sign,

```
|S|=[(tr S) S-2(det S) I]/sqrt[(tr S)²-4 det S].
```

Choose c=1/5 and r=2/3. Define

```
x0=(1-rho²)c,
epsilon=2 sqrt[rho²+(1-rho²)(r²-rho² c²)],
M=epsilon |S|-2x0 S,
h(k)=M+T exp(ik)+T^T exp(-ik).                           (3)
```

The three physical controls of BOND.md implement M exactly. Their
dimensionless values are approximately

```
mu=.018516885003506695,
eta=-.016156732294428884,
kappa=.007915197829885345.
```

Equations (1)--(3), followed by (3) of BOND.md, define the point;
the printed decimals are not its definition.

**[Proved finite algebra.]** With x=cos k,

```
det h(k)=4(-det S)(1-rho²)[r²-(x-c)²].                  (4)
```

Congruence by |S|^(-1/2) sends S to an involution of signature (1,1),
and sends A to an antisymmetric matrix with squared magnitude rho².
The transformed determinant is
`epsilon²-4[(x-x0)²+rho²(1-x²)]`. Multiply by det|S|=-det S
and substitute the definitions. Congruence preserves zeroes, not
physical eigenvalues, velocities or interaction vertices.

The four and only four roots are
`k=+/- arccos(-7/15), +/- arccos(13/15)`.
Both cosine roots lie strictly inside (-1,1), and (4) has nonzero
derivative there. The trace is positive for all k: use
`epsilon²[(tr S)²-4det S]>4(tr S)²`, tr S<0 and x0>0.
The other eigenvalue is nonzero, so the four crossings are simple.
All four occur in the **lower** Bloch band; the upper band stays
empty and gapped. There are two Dirac channels, not two separately
crossed orbital bands.

Actual eigenvectors and slopes are calculated in the original
orthonormal code basis:

| k | Eigenvalue slope |
|---:|---:|
| -2.05631444909 | +.000194934017244 |
| -.522314821806 | -.000464058241425 |
| +.522314821806 | +.000464058241425 |
| +2.05631444909 | -.000194934017244 |

The two velocities have ratio about 2.38. Their small absolute values
relative to the walk and T are emergent; they do not remove the
independent microscopic hierarchy in t0.

## Exactly the resonance under audit

Label the four sorted roots 0,1,2,3. All ten unordered pairs with
replacement are checked. The only collision of distinct pairs is
`k0+k3 = k1+k2 = 0 mod 2pi`. Other pair sums have minimum numerical
chord separation .405247. Repeated-fermion quartics vanish without
derivatives, but including their sums checks C84's stronger convention.
There is no other marginal umklapp or momentum-allowed inter-root
bilinear here. Physical U(1) excludes number-nonconserving terms.

Cooper-pair transfer between the two Dirac channels is therefore allowed.
It is retained, and its nonzero generated amplitude is used in
INTERACTIONS.md. Strict pair-sum separation is unnecessary for that
perturbative stability argument.

The simple roots and positive RG margin are strict conditions, allowing
an open neighbourhood of real onsite matrices by continuity. Nine
samples with `c in {3/20,1/5,1/4}` and `r in {3/5,2/3,7/10}` retain
four roots and positive leading RG margins. No certified rectangular
interval in all microscopic coefficients is claimed. Exploratory radii
1/3 and 1/2 also passed; increasing r improved the margins. The chosen
point does not tune pair transfer to zero.
