# Resonant leading models, scattering, and admissible detunings

Date: 2026-11-16. All microscopic Taylor coefficients below are constants
independent of L. Fix epsilon before taking any thermodynamic limit.
Finite-order elimination at fixed L is distinct from uniform infrared control.

Use the sibling's occupation notation: n=b†b, N=sum n,
D_r=sum n_j n_(j+r), T_r=sum b†_(j+r)b_j,
A=sum b†_j b†_(j+1), and
F=sum (n_(j-1)-n_(j+1))(b†_j+b_j).
We write **P_{not 1,2,4}** for the q=3 projector to avoid the ambiguous
P134 label. Every displayed leading operator is compressed on both sides.

## 1. Exact leading-model result

For q=1 use m=4,t=1 and P3 (separations 1,2,3 excluded). For q=3 use
m=3,t=128/3 and P_{not 1,2,4}. After canceling the scalar, chemical term
and nearest-neighbour hop in the sixth coefficient, define

```
K4 = -(T4+T4†)-2D4+(10/3)D5+(14/3)D6,
K3 = -t(T3+T3†)-2t D3+8G8,
G8 = sum n_j n_(j+8)(1-n_(j+3))(1-n_(j+5)).
```

Subscripts on K denote the hopping range, **not perturbation order**.
Both appear multiplied by epsilon^6. The actual chemical perturbation is
-mu6 N; the threshold is mu6=-2t. Write H_m=K_m+2tN.

**[Proved for these effective operators.]** On compatible rings L=mM,
large enough to distinguish the local distances, a hop by m cannot pass
another particle. The full cyclic ordered residue word is conserved,
not merely the m sublattice occupations. Boundary wrapping preserves
residues only for L divisible by m; we do not apply this assertion to
incompatible rings used for local coefficient checks.

Let A_graph count the allowed directed moves with multiplicity and let
L_graph=diag(degree)-A_graph. A gap blocks both approaches precisely at
4,5,6,7 for m=4 and at 3,5,7 for m=3. The other gap increases and remains
admissible. Therefore

```
H4 = L_graph +(16/3)D5+(20/3)D6+2D7,
H3 = t L_graph+2t(D5+D7)+8G8.                         (1)
```

Here D5,D6,D7 in the first line and D5,D7 in the second count consecutive
gaps. G8 also counts consecutive gaps, because its conditioning removes
the 3+5 and 5+3 decompositions. These are positive sums.

In every mixed residue word, reduce a gap between unlike neighbours by m
until it reaches its minimum: 5/6/7 for m=4, 5/7 for m=3. Each preceding
move is legal. The connected component thus reaches a positive-potential
vertex. A Laplacian zero vector is constant on a component, so no such
component belongs to the kernel of (1). On a single residue class the
coarse chain is exactly

```
t sum (I-SWAP) = -t(T1+T1†)-2t D1+2t N.              (2)
```

Its fixed-N ground vector is the uniform hard-core occupation state.
The full kernel is one vacuum plus m states for each N=1,...,M:
**dim ker H_m=L+1**. At each fixed N>0 the leading ground states select
single-residue sectors (or superpositions of those sectors), not arbitrary
degenerate mixed-color words. A mixed-color interface can have vanishing
energy as density tends to zero; no uniform color gap is asserted.

The grand ground energy of K_m-mu6 N is minimized at N=0 if mu6<-2t and
at N=M if mu6>-2t, since H_m>=0 and both endpoints attain zero. At mu6=-2t
every N is degenerate. Thus the naive leading construction sits at a
density jump/ferromagnetic threshold, not a stable compressible liquid.
The one-particle energy above threshold is 2t(1-cos k) on the coarse chain,
quadratic near k=0. Counting four/three minima as four/three positive-velocity
conformal modes is invalid.

These identities and kernel counts were independently checked by the
coordinating chat with integer arithmetic in
[its review](../research_coordination_20261116/README.md). Our `leading.py`
independently constructs the projected hopping and diagonal operators;
`check.py` additionally compares the full matrices to the sibling's generic
microscopic Schur recursion. That last recursion is shared code, not an
independent verification of the recursion itself.

## 2. Two-body channels at the undeformed point

**[Derived exactly.]** Use total coarse momentum Q and ordered separation
r=b+m n, n>=0. Both particles have hopping -t. Relative hopping has magnitude
J(Q)=2t cos(Q/2) for -pi<Q<pi. In a same-residue channel b=m, the only
potential is V_0=-2t. The half-line equations are

```
E f_0=V_0 f_0-J f_1,
E f_n=-J(f_(n-1)+f_(n+1))  (n>=1).
```

At Q=0, |V_0|=J: the zero-momentum bound state is exactly at threshold,
with a non-normalizable constant relative wavefunction. At Q != 0 a
normalizable same-residue pair exists, with

```
z=J/(2t)=cos(Q/2),
E_bound=-2t-J²/(2t),
E_continuum_min=-2J,
E_bound-E_continuum_min=-2t[1-cos(Q/2)]².             (3)
```

Its global energy never lies below -4t. Binding near Q=0 is quartic in Q;
the center-of-mass branch is quadratic. Thus “no bound pair” would be too
strong: there are finite-momentum ferromagnetic strings, but no negative
binding below the global two-free-particle minimum at the threshold.

Different ordered residues give the following half-lines. All their
potentials are nonnegative, so none has a bound state below its continuum.
At Q=0 write the asymptotic zero-energy wave as f proportional to r-a.

| q | Ordered residue difference | Minimal b | Potential | Scattering length a, physical sites |
|---|---:|---:|---|---:|
| 1 | 1 | 5 | V(b)=10/3 | 7/2 |
| 1 | 2 | 6 | V(b)=14/3 | 24/5 |
| 1 | 3 | 7 | none | 3 |
| 3 | 1 | 7 | none | 4 |
| 3 | 2 | 5 | V(8)=8 | 56/19 |

For an endpoint potential, a=b-m J/(J+V_0). For q=3's interior barrier
v=8/J=3/32, f_n is proportional to n+(1-2v)/(1+2v) beyond the barrier,
which gives 56/19. These are oriented infinite-line channels. On a ring
the two ends connect the two orientations; finite-ring spectra must retain
both. `scattering.py` verifies (3) and the absence of mixed-channel bound
states on a 2000-site relative half-line. No entropy fit is used.

## 3. The alpha detuning contains an essential off-diagonal term

Set p(epsilon)=q epsilon²+alpha epsilon³+beta epsilon⁴, with q=1 or 3.
Maintain the lower scalar/chemical/hopping cancellations. Through order six
the required microscopic visible coefficients, before adding -mu6 N, are

```
u2=2, t2=-1,
u4=-(4+4q+2q²),                  t4=2+2q,
u5=-(4+4q) Re(alpha),            t5=2 Re(alpha),
u6=cN-2|alpha|²-(4+4q)Re(beta),  t6=cT+2 Re(beta),
(cN,cT)=(50,-24) at q=1; (332/3,-80) at q=3.
```

The scalar addition is
L[q² epsilon⁴+2q Re(alpha)epsilon⁵+
(|alpha|²+2q Re(beta)-cL)epsilon⁶], with cL=2 or 18.
The scalar only fixes energy zero. Source jets through order five are zero;
the sixth-order source is set to zero for this test. All coefficients are
fixed numbers; beta and alpha are not tuned with L.

**[Derived finite-order identity, full-matrix checked through L=13.]** The
fifth coefficient vanishes. The new sixth coefficient is

```
K_m(alpha)=K_m(0)+d |alpha|² D_m
  + kappa sum [(n_(j-m)-n_(j+m)) alpha* b†_j+h.c.],    (4)

(m,t,d,kappa)=(4,1,1,1) or (3,128/3,2/3,16/3).
```

Beta cancels from (4). The extra diagonal term follows from the exceptional
fourth penalty: |q-1|²D4 or (2/3)|q-3|²D3. **It is not the full answer.**
The linear-alpha, number-changing term is equally leading. In particular
<0,4|K4(alpha)-K4(0)|0>=alpha* and
<0,3|K3(alpha)-K3(0)|0>=(16/3)alpha*.
The sign reverses for creation on the other side. Restricting the microscopic
elimination to fixed N would wrongly turn (4) into a simple XXZ detuning.

There is an exact consequence of these checked coefficients. Put
g=kappa alpha*. They satisfy d|alpha|²=|g|²/t. On each coarse bond,

```
v=(0,-sqrt(t),sqrt(t),g/sqrt(t)),
h_alpha=|v><v|,
sum h_alpha=-t(T+T†)+(-2t+|g|²/t)D1+2tN
            +sum[g(n_(j-1)-n_(j+1))b†_j+h.c.].       (5)
```

**[Proved for (4).]** Globally K_m(alpha)+2tN is the compression of the
positive sum (5) on the m residue chains, plus the original nonnegative
D5,D6 or G8 penalties. It is therefore positive. For alpha !=0, each
coarse periodic chain has precisely the vacuum and uniform one-particle W
state in its kernel. To see exhaustion, its forbidden-vector transfer matrix
is a nontrivial Jordan matrix with a nonzero repeated eigenvalue. Its M-th
power is not scalar, so Bravyi--Gosset Proposition 3 gives dimension two;
the two stated vectors visibly span it. This is an application of their
[periodic ground-space theorem](https://arxiv.org/html/1503.04035#S2.SS2),
not of a conformal classification.

The unconstrained m-chain kernel is a tensor product of those two-dimensional
spaces. Any product containing W on two different residues has nonzero
amplitudes on forbidden nearby particle pairs. Different occupied-residue
subsets have disjoint occupation support and cannot cancel them. Intersecting
with the physical constrained space leaves exactly

```
vacuum, W_0, ..., W_(m-1);    dim ker = m+1.           (6)
```

Thus this alpha detuning **does not** create a finite-density multicomponent
liquid at the shifted threshold. It does not even conserve N. The reduced
cyclic residue word, with consecutive repeated labels erased, is still
conserved: (4) only inserts/deletes a duplicate neighbour's label. It does
not create a new residue color or exchange different ordered labels.
At other chemical potentials the parent argument no longer decides the
phase; no exclusion of all such detunings is claimed.

## 4. Actual label changes and their scale

At order epsilon^6, a residual microscopic t6=t6_required+tau gives
tau T1+tau* T1† directly in the constrained model. This genuinely changes
the residue of a particle. A source f6 gives a uniform creation term and
also mixes colors/particle numbers. Neither damages the order-two or
order-four blockades; both change the leading sixth dynamics rather than
canceling it. The coefficients are freely chosen microscopic nearest-
neighbour jets. Generic tau splits the original m minima, so this is not
automatically a route to c=m.

Even with tau=f6=alpha=beta=0, the next static Schur coefficient has
residue-changing, number-changing matrix elements. In our fixed occupation
basis and the counterterm convention above, the L>=14/15 checks give

```
<0,5|K7|0> = 13/3       (q=1),
<0,5|K7|0> = -128/3     (q=3),                         (7)
```

with the opposite sign for creation five sites to the left. These are
epsilon^7 processes: O(epsilon) relative to the sixth-order kinetic scale.
Short rings alias the local processes; all L-specific raw coefficients
are retained, including those that differ from (7).
The seventh matrix is Hermitian in these checks. Effective normalization
first contributes at relative order epsilon² here, so it does not alter
the displayed seventh coefficient. A uniform f7 cannot cancel both signs;
t7 cannot cancel a number-changing matrix element.

**Important limit:** a nonzero off-block matrix element is not a new
gapless color velocity. On the checked compatible rings, projection of K7
onto the entire L+1-dimensional zero manifold of H_m vanishes to rounding.
It couples that manifold to excited configurations. An induced splitting
inside it requires a higher-order calculation with the energy-dependent
metric and degenerate elimination. Also, effective basis changes can move
off-block terms between orders; (7) identifies the process in this specified
occupation-basis convention, not a proof that no locally dressed conserved
labels exist. That remaining step is not claimed complete.

## 5. Scope of the obstruction

The undeformed m=3,4 leading models are exactly ferromagnetic thresholds;
the alpha-only cure retains a positive parent structure and has only
vacuum/one-particle zero states at the shifted threshold. These are clean
obstructions to the two naive multivalley arguments. They are not a
no-go theorem for every fixed finite-epsilon microscopic detuning. In
particular, tau, chemical filling, later non-number-conserving processes,
and their thermodynamic control remain substantive infrared questions.
No microscopic c>1 or universal upper bound is established here.
