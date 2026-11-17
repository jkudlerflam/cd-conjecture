# Further projection, parent identities, and the remaining obstruction

The two classical branches and notation are defined in
[CLASSIFICATION.md](CLASSIFICATION.md). All expansions here are formal
finite-order expansions at fixed ring length. A claim about the complete
thermodynamic infrared theory is deliberately separate.

## 1. Which effective operator is being computed

Let V1=F, V2=W, `A=Q F P`, `B=Q W P`, `C=Q F Q`, `D=Q W Q`,
`R=(Q H0 Q)^-1`. With the zero-order space P identified with its occupation
basis, the static Schur complement has coefficients

```
S2 = P W P - A† R A,
S3 = P V3 P + A† R C R A - A† R B - B† R A,
S4 = P V4 P - B† R B
     + A† R C R B + B† R C R A + A† R D R A
     - A† R C R C R A
     - A† R Q V3 P - P V3 Q R A.                              (6)
```

First shift the scalar energy at every lower order. If I is the kernel
of S2 and J its complement within P, the next coefficients are

```
K3 = I S3 I,
K4 = I S4 I - I S3 J (J S2 J)^-1 J S3 I.                     (7)
```

When K3=0, (7) is also the leading Hermitian effective Hamiltonian.
Metric/normalization terms have zero compression because S2 I=0.
Keeping only `I S4 I` gives the wrong answer. The same reasoning applies
to the leading fifth/sixth coefficient after all lower coefficients are
scalar and removed. We do **not** interpret arbitrary higher static Schur
coefficients as a Hermitian effective Hamiltonian when a lower nonzero
operator survives.

For reproduction at any order, solve `D(epsilon)X(epsilon)=B(epsilon)`
in the eliminated space. The recursion is

```
X_n = R [B_n - sum_(k=1)^(n-1) D_k X_(n-k)],
S_n = A_n - sum_(k=1)^(n-1) B_k† X_(n-k).                    (8)
```

Apply (8) once to H0 and again to the order-two constraint. The code uses
the entire constrained bases, without a fixed-N truncation.

## 2. Dilute branch through fourth order

Use `W=2N-T1-T1†+q A+q* A†`. On P2, `K3=P2 V3 P2`.
The visible third-order coefficients must therefore obey
`u3=t3=f3=0` to keep this whole space scalar. The arbitrary pair p3 has
no effect on K4. After the rotations/normalization of the normal form,
these are all the third-order possibilities.

**[Derived local algebra; exact polynomial matrix check]** Write q=x+iy.
Before adding the visible fourth-order counterterms, the complete answer is

```
K4 = P2 { -L |q|² + (4+4x+2|q|²)N
          -2(1+x)(T1+T1†)
          +(2/3)|q-3|² D3 + |q-1|² D4 } P2.                 (9)
```

In particular there is no number-changing term at this order. This is an
operator identity on all sectors, not a one-particle compression.

One way to check the local path calculation is to separate the two pieces
of (7). Put `C2=P sum b†_j(n_(j-2)-n_(j+2)) P2`. Then

```
S3 P2 = (1-q) C2,
C2† C2 = P2 [2N-2D3-2D4-T2-T2†] P2,
```

and all states reached by C2 have S2 energy 2. The direct S4 compression
has the following coefficients, found by the four-step paths in (6):

| Operator | I S4 I |
|---|---:|
| L | `-|q|²` |
| N | `5+2x+3|q|²` |
| T1+T1† | `-2(1+x)` |
| T2+T2† | `-|q-1|²/2` |
| D3 | `5-2x-|q|²/3` |
| D4 | 0 |

The subtraction `-|q-1|² C2†C2/2` gives (9), including the cancellation
of the entire distance-two hop. These rules also show why discarding
number-changing intermediate states would fail.

The conditions for the fourth operator to be diagonal are

```
t4=2(1+x) (real),  f4=0.                                    (10)
```

Choosing `u4=-(4+4x+2|q|²)` removes its chemical term. Its positive
constraint is then

```
(2/3)|q-3|² D3 + |q-1|² D4.                                 (11)
```

It is never scalar on all of P2: q=3 and q=1 cannot hold simultaneously.
The next kernels are different:

* q different from 1,3: no separations 1,2,3,4 (P4).
* q=1: no separations 1,2,3 (P3).
* q=3: no separations 1,2,4; separation three remains allowed.

Other u4 values select different densities. Formula (9) is complete, but
it does not classify all the subsequent dynamics on these kernels.

## 3. A general exact rank-one parent, and its exceptions

**[Proved elementary identity]** Let a,c be complex and b real, and set

```
psi=(c,-a,a,b) in the ordered basis 00,01,10,11,
v=b²+|c|²-2|a|²,       u=2(|a|²-|c|²),
t=-|a|²,               p=b c*,
g=b a*+a c*,           f=0,
G=[[0,c a*],[a c*,0]].
```

The literal nearest-neighbour bond with facilitations `(g,-g)` satisfies

```
h + |c|² I + G tensor I - I tensor G = |psi><psi|.             (12)
```

The difference telescopes on rings. Thus this is a uniform rank-one
frustration-free parent up to a constant; no size-dependent bond is used.
For v=1, g=epsilon real and p=q epsilon², put z=b² and solve

```
1 = z + |p|²/z - 2 z epsilon² |z-p|²/(z²-|p|²)²,
a = sqrt(z) epsilon (z-p*)/(z²-|p|²),  c=p*/sqrt(z).           (13)
```

The derivative of the first equation with respect to z at
`epsilon=0,z=1` is nonzero; the real implicit-function theorem gives
a unique local solution. Expanding it yields

```
t=-epsilon²+2(1+Re q)epsilon⁴+O(epsilon⁶),
u=2epsilon²-(4+4Re q+2|q|²)epsilon⁴+O(epsilon⁶).             (14)
```

Thus every cancellation in (10)--(11), with its chemical tuning, is
tangent through fourth order to (12). Higher analytic pair amplitudes
p(epsilon) can also be inserted in (13). This proves tangency, not that
arbitrary higher jets are on the exact manifold.

**[Read, with application]** Bravyi--Gosset's transfer matrix here is

```
Tpsi=[[-a*,b],[-c*,-a*]].                                    (15)
```

When it is invertible, their Proposition 3 says that the periodic ground
space has dimension two unless `Tpsi^L` is proportional to I, when it has
dimension L+1. For the nonresonant diagonalizable case it is spanned by
two product states. The defective case is obtained from vacuum and W
under a product of invertible onsite maps, hence also has Schmidt rank at
most two for every ground state. These conclusions use the construction
in their Propositions 2--3, not just the ground-space dimension.

Product forbidden vectors (`bc+a²=0`) require separate treatment, since
Tpsi is singular. For nonparallel factors, the invertible two-state basis
argument forbids one ordered pair; on a ring only the two constant words
survive. This recovers the earlier q=-1 continuation. The q=0 continuation
has a nontrivial Jordan block and recovers the earlier vacuum/W result.

**Caveat:** the general family also includes roots-of-unity eigenvalue
ratios in (15). The enlarged ground spaces on compatible rings do not
inherit the preceding rank-two argument. There is no claim that every
ground state of every rank-one parent has uniformly bounded entanglement.
Nor does their Theorem 1 identify a conformal universality class or bound
all nearby detunings. This is why (12) alone cannot close local boundedness.

Primary source: [Bravyi--Gosset, Sections 2.1--2.2 and Theorem 1](https://arxiv.org/html/1503.04035).

## 4. Dense branch: a different fourth-order cancellation locus

Use `W=-4N-T1-T1†+q A+q* A†`, and add `2L epsilon²` to the energy.
Then S2=2 Z3 on P, with kernel Pd (spacings two or three).

Here `Pd B Pd=0`: adding a particle is impossible and removing one
creates a run of three zeros. The entire third-order scalar condition is
`u3=t3=0`; a complex uniform source **f3=f is free**. The arbitrary p3
again does not affect K4. In this section f denotes this third-order
source, not the first-order f in the microscopic normal form.

**[Derived local algebra; exact polynomial matrix check]** Let
`Q=|q|²`, `F=|f|²`, `x=Re q`. The fourth-order operator is

```
K4=Pd [C L + B N + A4 D4 + tau T1 + tau* T1†] Pd,             (16)
C = 1/2-x+(Q+F)/6,
B = -7/2+3x-5Q/6-2F/3,
A4=|q-1|²/2-F/6,
tau=(|q-1|²-F)/4 + (i/2) Im[(3+q*)f].
```

To verify (16), insert the dense W and source in (6), shift the energy
before eliminating Q, and use `R2=(2 Z3)^-1` on P minus Pd. Classify the
allowed words by their neighboring gaps 22,23,32,33. The only surviving
off-diagonal move is 23<->32 (T1); all Delta N not equal to zero cancel.
The diagonal possibilities reduce to L,N and the number D4 of neighboring
22 gaps. Equation (16) gives their complete polynomial coefficients.
`symbolic_check.py` checks every matrix entry with four free real symbols;
larger rings check arbitrary spectators and longer tile strings.

All off-diagonal dynamics at this order can be removed with `t4=-tau`.
For the **entire Pd space** to remain scalar, the other conditions are

```
u4=-B,    |f|²=3|q-1|².                                    (17)
```

The residual scalar is CL. f4 and p4 are invisible at this order.
Equation (17) is a genuine branch of further cancellation. It cannot be
tangent to a positive rank-one parent approaching H0: such a parent's
second jet has u2=+2 after the first cancellation, whereas this has u2=-4.
Adding a scalar or telescoping bond does not change this difference.

For A4 nonzero, (16) still gives classical fourth-order subproblems after
cancelling the hop. Their lower-dimensional kernels and chemical
coexistences are **not** covered by the scalar-Pd obstruction below.

## 5. Dense scalar branch at orders five and six

First take the representative with only the jets specified in (17),
and subtract the scalar CL at order four. Set all remaining microscopic
coefficients to zero. **[Derived local coefficient formulas; exact
symbolic finite-ring verification, plus larger-ring checks]** Its entire
fifth-order operator is a tile conversion

```
33 <-> 222,
1001001 -> 1010101,
gamma=(q+3)f*/6+(q*-1)f/2,                                  (18)
```

with gamma for the displayed direction and gamma* for its reverse.
The occupied endpoints are retained. The formula sums over all positions
with arbitrary legal spectators. This operator changes particle number
by one and would be absent from a fixed-N calculation.

Together with (17), gamma=0 gives exactly

```
(q,f)=(1,0),
or |q|²=3 Re q,  |f|²=3|q-1|²,
   (q+3)f*+3(q*-1)f=0.                                    (19)
```

The second set is the circle `|q-3/2|=3/2` with two permitted source
phases. Write q=x+iy and f=r+is. On it,

```
0<=x<=3, y²=x(3-x), r²=(3-x)(x+1), s²=x(x+1),
xr+ys=0, (3-x)s+yr=0.                                     (20)
```

At order six consider the correlated move

```
223 -> 322,
10101001 -> 10010101.                                      (21)
```

It moves two interior particles; a microscopic nearest-neighbour bond
cannot have this direct matrix element. After the two required
eliminations, its coefficient is

```
-1/2                                      at (q,f)=(1,0),
(15x²+4x-3)/16 + i(5x-1)s/4                on the circle.    (22)
```

The latter never vanishes: its imaginary part can vanish only at x=0 or
x=1/5, where its real part is -3/16 or -1/10, respectively. At x=3 the
imaginary part is nonzero. Equations (18), (22) are obtained using the
path recursion (8). Exact rational/polynomial checks, including reduction
by the algebraic relations (20), are saved separately from floating-point
checks. Finite rings L=9,10,12,13 reproduce the local coefficient.

**Effect of remaining microscopic jets.** In the canonical gauge write
`p=epsilon² q(epsilon)`, `source=epsilon³ f(epsilon)` and allow analytic
u,t. Requiring the whole Pd space to remain scalar through fifth order
forces its earlier u,t coefficients as above. Variations in the next
q,f coefficients differentiate (16), and hence contribute only L,N,D4,T1
at the next orders. Differentiating (18) contributes only the conversion.
The direct u6,t6 counterterms contribute N,T1; source6 and pair6 have
zero Pd compression. None contributes the four-bit change (21).
Thus these unused jets cannot cancel the witness (22). This reasoning
concerns scalar cancellation on the *whole Pd space*, not an arbitrary
selection of a smaller ground manifold by earlier diagonal terms.

The local derivation and exact finite-ring checks establish a concrete
finite-order obstruction to that scalar-cancellation strategy. The
checks are not a proof of uniform convergence of the expansion with L,
nor of a bound on the complexity of the sixth-order interacting model.

## 6. What remains for a finite upper bound

There is no claim of an all-orders obstruction to every nesting route.
In particular, the following remain:

1. Dynamics on P4, P3, and the `no 1,2,4` dilute kernel, including jets
   away from the exact parent family and its resonances.
2. Dense classical fourth-order branches with A4 nonzero, and their
   additional chemical coexistences. They can select a smaller space
   before the scalar-Pd obstruction applies.
3. Analytic arcs with different leading orders, nonanalytic parameter
   approaches, and bonds outside this neighborhood.
4. A uniform thermodynamic construction of the constrained effective
   theories and their complete infrared degrees of freedom, including
   composite liquids, neutral critical sectors, and very small velocities.

Even a proof that every branch develops motion by a bounded order would
not itself bound interacting c. A range or valley count is not a count
of all conformal sectors.

The exact missing local lemma for compactness is: **there exist an open
neighborhood U of the normalized blockade bond and a finite B such that
every fixed bond in U with the stipulated complete unitary conformal
bulk limit has total c<=B.** Equivalent neighborhood statements are
needed around every normalized bond, including nonconformal limits.

A sufficient, stronger route would construct a finite collection of
uniformly controlled continuum descriptions covering this neighborhood,
prove that they include all low-energy sectors, and bound the total c
in each. Another sufficient route is the one-sided finite-interval
entropy estimate in the compactness note, with common finite constants.
This contribution proves neither premise. No conclusion c<=2, or even
finite universal B(2), follows from the present perturbative classification.
