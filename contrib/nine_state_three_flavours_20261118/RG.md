# The complete leading interaction and the coupled growth cone

## Generated coefficient, not an independently chosen interaction

Let R0 be the inverse of the two free charged-carrier Hamiltonians,
C the sum of their literal conversions without g, and S0 the literal
port source without lambda alpha². In gauged species order AA,AB,BA,BB,
S0 has components `(0,-P_left,I,0)` at carrier positions (4,0). Then

```
K(g)=-S0^dagger (R0^(-1)+g C)^(-1) S0
    =K0+g² K2+g⁴ K4+O(g⁶),
K_(2r)=-S0^dagger R0 (C R0)^(2r) S0.                  (1)
```

K0 is scalar, K2 is quadratic on the entire 64-dimensional pair code,
and K4 is exactly reconstructible as constant, quadratic and quartic
pieces. Odd powers vanish by species return. The literal excited-pair
matrix has dimension 6400. Independently assembling it from physical
words gives agreement through K4 within 5.6e-17. The added weak physical
controls contribute M directly; their insertions inside (1) start at
higher port order, not at the leading lambda² g⁴ quartic.

Normal order using the six canonical fermions of the two cells:

```
K4_quartic=sum_(i<j,k<l) W_(ij),(kl) f_i^dagger f_j^dagger f_l f_k.
```

Subtract the vacuum constant and full one-particle quadratic lift from
the two-particle block to obtain W. Reconstructing every particle sector
checks that no higher-degree part of K4 was omitted. All 225 real entries
of W are saved, along with its pair ordering and the K4 quadratic part.

For a normalized physical-band eigenvector u_a at k_a let
`w_a=(u_a,exp(i k_a)u_a)`. The full antisymmetrized coefficient is

```
Gamma_(ab;cd)=(w_a wedge w_b)^dagger W (w_c wedge w_d). (2)
```

This includes direct and exchange contributions. No nonunitary change
of basis enters the matching. All coefficients below multiply
`lambda² alpha⁴ g⁴=t0 g²`; velocities multiply t0.

## Every marginal channel at the selected point

In the channel order of BANDS.md, define G_i as the R_i,L_i density
coupling, F_ij as the R_i,L_j density coupling for i!=j, and J_ij as
the coefficient of
`R_i^dagger L_i^dagger L_j R_j+h.c.`. Reciprocity gives F_ij=F_ji.

| Channel | Coefficient |
|---|---:|
| G1 | +.0000563980448694 |
| G2 | -.0004920800841751 |
| G3 | -.0005912988497160 |
| F12 | +.0000154277137697 |
| F13 | -.0002996705879118 |
| F23 | -.0000013765687524 |
| J12 | +.0001828408039589 |
| J13 | -.0002899336076042 |
| J23 | +.0005629133396652 |
| Same-chiral density, channels 1,2 | -.0001674130901892 |
| Same-chiral density, channels 1,3 | -.0000097369803076 |
| Same-chiral density, channels 2,3 | -.0005642899084176 |

The same-chiral values apply to both R,R and L,L. Altogether there are
15 diagonal density vertices and six directed pair-transfer vertices,
including Hermitian conjugates. They are all saved in `audit.json`.
The complete pair-sum audit finds no other momentum-allowed quartic at
this point. Physical charge excludes number-nonconserving terms.

The signs of individual J depend on field phases; the real triangle
product `J12 J13 J23` is invariant under real channel rephasings.
It is negative here. No independent channel U(1) prohibits these transfers.

## Coupled one-loop equations

In this subsection G,F,J and v are physical continuum couplings and
positive velocities. With ell the logarithmic increase in length, the
leading equations are

```
dG_i/dell = -sum_(k!=i) |J_ik|²/(2pi v_k),
dF_ij/dell = |J_ij|²/[pi(v_i+v_j)],
dJ_ij/dell = -D_ij J_ij
             -sum_(k!=i,j) J_ik J_kj/(2pi v_k),
D_ij = G_i/(2pi v_i)+G_j/(2pi v_j)
       -2F_ij/[pi(v_i+v_j)].                         (3)
```

The last sum is essential at three flavours. A Cooper contraction
through channel k gives the shell integral 1/(2pi v_k), the matrix
product J_ik J_kj, and the minus sign of the second cumulant.
The two mixed-channel particle-hole contractions yield the forward
term in D. Self-channel Cooper and particle-hole terms cancel for the
pure density interactions. These contractions can equivalently be
organized with the spinless chiral currents R_i^dagger R_j and
L_i^dagger L_j. Same-chiral densities shift velocities at first order
and do not contribute a logarithm to (3).

For normalization and sign checks, the two-channel restriction agrees
with the n=1 specialization of Szirmai--Solyom,
[Eqs. (21), (23)--(25), (27)](https://arxiv.org/pdf/cond-mat/0605300),
using their parallel coupling equal to minus our density coupling.
The intermediate-channel contraction has the Cooper matrix structure
derived by Lin--Balents--Fisher in
[Section III and Appendix A](https://arxiv.org/pdf/cond-mat/9703055).
That paper treats spinful ladders; (3) is the spinless contraction written
here, not an identification of its Hubbard model with this bond.

Normalize

```
b_i=G_i/(2pi v_i), f_ij=F_ij/[pi(v_i+v_j)],
c_ij=J_ij/(2pi sqrt(v_i v_j)),
a_ij=4v_i v_j/(v_i+v_j)²,  D_ij=b_i+b_j-2f_ij.
```

Then b_i'=-sum c_ik² and f_ij'=a_ij c_ij². For three flavours,
in edge order 12,13,23, (3) closes as

```
D_ij'=-2(1+a_ij)c_ij²-c_ik²-c_jk²,
c_ij'=-D_ij c_ij-c_ik c_kj,                          (4)
```

where k is the third channel. `audit.py` compares (4) to an independent
matrix implementation retaining all forward terms. Removing the third
channel recovers the earlier two-flavour KT equations, but that
restriction is not used to decide this point.

## A strict invariant cone, not a solver divergence claim

At the primary point, after removing the common g² from the dimensionless
initial data,

```
D=(-195.766578055,-7.451235561,-362.961693129),
c=(53.842720798,-70.814527202,189.004242770),
a=(.905193268,.982984133,.965816372).
```

Set tau=g² ell and evolve these coefficient vectors by (4). Since the
triangle product is negative, changing channel pair signs by (1,-1,1)
makes every c negative. Put x_ij=-c_ij>0. Equation (4) becomes

```
x_ij'=-D_ij x_ij+x_ik x_kj,
D_ij'<=0.
```

**[Proved property of the leading ODE.]** The cone x>0, D<0 is invariant.
Every transfer grows, including its coupling through the third channel;
all forward dimensions become more relevant. With
`kappa=min(-D(0))=7.45123556` and `q=min(x_ij)`, the lower Dini derivative
satisfies `q'>=kappa q+q²`. Thus the formal one-loop solution cannot stay
bounded beyond

```
tau <= log(1+kappa/q(0))/kappa = .0173949712.           (5)
```

This proves a leading weak-coupling instability from strict initial
inequalities; it does not rely on extrapolating a numerical integrator
to a pole. The numerical trace stops at tenfold coefficient growth,
tau=.00232580153. At sufficiently small g all physical couplings can
still be weak there. One-loop equations cease to justify microscopic
predictions once the physical couplings become strong.

The numerical matching of the initial inequalities is distinct from the
ODE proof. Their margins are large compared with floating errors. All
32 corners of the band-certified control box also have negative D and
negative triangle product; the smallest sampled negative-D margin is
1.3117. This corner check is not an interval proof for the nonlinear
band matching throughout the box. Strict inequalities at the central
point do imply an open instability neighbourhood by continuity.

The resulting conclusion is an unstable three-channel Gaussian theory
under the stated matching/RG assumptions. We do not assign the eventual
central charge, prove relative-mode gaps, or assert instability of the
whole five-control family.
