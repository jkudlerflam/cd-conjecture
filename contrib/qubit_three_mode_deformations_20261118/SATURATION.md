# A local certificate of crystal preemption

Status: **[proved, with an exact rational arithmetic certificate]** for
the displayed complete leading Hamiltonian. This is an all-length
operator inequality, not a fit of finite-ring ground energies. It is
not yet a microscopic finite-epsilon theorem.

Use (7) of [NORMAL_SPACE.md](NORMAL_SPACE.md), set a=eta=0, and measure
energy in units t3=128/3. Thus

```
H = P[sum h0(j,j+3)+(3/16)G8+zeta T1+zeta* T1†-mu N]P,
h0 = |01-10><01-10|,
gamma = |Re(zeta)|+|Im(zeta)|.                               (1)
```

P forbids separations 1,2,4. Work on periodic rings L divisible by three,
L>=12. N<=M=L/3. The three fully occupied residue classes are called
the period-three crystals; each has energy -mu L/3.

## 1. Exact saturation and gap theorem

**Theorem.** If

```
0 < gamma <= 1/10,              mu >= (7/10) gamma,          (2)
```

the entire ground space of (1) is precisely the span of the three
period-three crystals. For strict inequality in the second condition,

```
gap_L >= mu-(7/10)gamma.                                    (3)
```

The density is fixed to 1/3 throughout this region; its compressibility
vanishes. Bound clusters, a polarized lower-density gas, or any mixed
state cannot undercut this energy, because the certificate bounds the
whole Hilbert space. The ground degeneracy is discrete translation
breaking, not three conformal modes.

### The twelve-site block

Let P_open forbid physical separations 1,2,4 within sites 0,...,11, with
no wrapping. Define

```
Q12(zeta,mu) = P_open [
    (1/9) sum_(j=0)^8 h0(j,j+3)
  + (3/64) sum_(j=0)^3 n_j n_(j+8)(1-n_(j+3))(1-n_(j+5))
  + (1/11) sum_(j=0)^10 [zeta b†_(j+1)b_j+h.c.]
  + mu(1/3-N12/12)] P_open.                                (4)
```

When summed over all translated blocks on a constrained ring, each
distance-three bond is counted nine times, each distance-eight term four
times, each nearest-neighbour hop eleven times, and each onsite number
twelve times. Compressing once more to the full ring constraint gives
the exact identity

```
sum_x P Q12(x) P = H+mu L/3.                                (5)
```

Local states or transitions forbidden by particles outside a block are
removed by the outer P; positivity survives this further compression.
Blocked-hop diagonal entries in h0 are retained. Replacing h0 with only
the allowed hopping graph would invalidate (4).

At (zeta,mu)=(0,0), (4) is visibly positive. Exact rational LDL
elimination certifies positivity at the four additional vertices

```
(zeta,mu) = (+1/10,7/100), (-1/10,7/100),
            (+i/10,7/100), (-i/10,7/100).                   (6)
```

The negative imaginary vertex is the complex conjugate of the positive
one. The five number sectors have dimensions **1,12,37,30,6**. At each
nonzero vertex, all N<4 blocks are strictly positive; the N=4 block has
rank three, with kernel exactly the three local crystal words. For the
imaginary vertex, realification doubles those dimensions and ranks.
The largest rational matrix tested has dimension 74.

`block_certificate.py` constructs these rational matrices directly from
occupation words. `check.py` performs ordinary Fraction arithmetic:
every LDL pivot is strictly positive except the declared crystal zeros,
and each zero pivot has an exactly zero remaining column. This is a
complete PSD test, including singular cases, with no numerical tolerance.
The input hashes, ranks and smallest positive pivots are in `checks.json`.
The proof relies on this small reproducible exact computation; it is
not presented as a hand diagonalization or independent human review.

For general zeta in (2), put x=Re(zeta), y=Im(zeta). Use weights
`10|x|`, `10|y|`, `1-10 gamma` on the correspondingly signed real
vertex, imaginary vertex and origin. Their convex combination is exactly
`Q12(zeta,(7/10)gamma)`. It is positive. At least one nonzero vertex has
positive weight, so its kernel is exactly the three local crystal words.
Increasing mu adds the positive operator `(mu-(7/10)gamma)(4-N12)/12`
and leaves that kernel unchanged.

Equation (5) now proves the energy bound. A global zero vector must
restrict to a local crystal word on every twelve-site window. Overlapping
windows fix the same residue, so their intersection is exactly the three
global crystals. This proves the ground-space assertion, including the
boundary mu=(7/10)gamma.

Finally, N is conserved. At mu0=(7/10)gamma,
`H(mu0)+mu0 M>=0`; consequently

```
H(mu)+mu M >= (mu-mu0)(M-N).                                (7)
```

The N=M space already consists only of the three crystals. Every other
number sector has M-N>=1, proving (3). No volume-independent gap was
assumed in obtaining this bound.

As a separate numerical implementation check, the literal block sum was
embedded in complete periodic spaces L=15,18 and compared with the
independent leading-model code. The largest discrepancy was 1.78e-15.
The all-length identity itself follows from the multiplicities in (5).

## 2. Three equal dilute thresholds require zeta=0

The one-particle compression of (1), before its chemical shift, is

```
e(k) = 2-2 cos(3k)+2 Re(zeta exp(-ik)).                      (8)
```

**Lemma.** Three distinct equal global minima of (8) require zeta=0.

Proof. Set z=exp(ik), and let E be their common value after subtracting
the constant two. The degree-six polynomial

```
z³[e(k)-2-E] = -z⁶+zeta* z⁴-E z³+zeta z²-1                (9)
```

has a double zero at each minimum. It must therefore equal
`-(z³-s1 z²+s2 z-s3)²`. The absent z⁵ coefficient gives s1=0;
the absent z coefficient gives s2=0 because s3 is nonzero. The z⁴ and
z² coefficients then force zeta=0. At zeta=0 the exact leading parent
has the known ferromagnetic threshold and density jump, not a stable
three-component dilute gas. The lemma concerns equal **global** minima;
there can certainly be three unequal local minima.

## 3. Filling the higher third valley is preempted in the certified region

One proposed escape is to keep the unequal valleys and raise the chemical
potential until all three one-particle wells are occupied. The crystal
theorem rules out this region of the complete interacting model, not
just a two-body approximation.

Let q_a=2pi a/3, a=0,1,2, and minimize (8) on each interval
`|k-q_a|<=pi/3`. For gamma<=1/10, each minimum lies within
`|k-q_a|<pi/6`: outside that smaller interval the parent part is at least
2 and the perturbation at least -2gamma, whereas the center has energy
at most 2gamma. There are three distinct local wells. Let mu3 be the
largest of these three well minima. Write rho=|zeta|<=gamma. On the
smaller interval, with delta=k-q_a,

```
2(1-cos(3delta)) >= (72/pi²)delta²,
2 Re[zeta(exp(-ik)-exp(-iq_a))] >= -2rho |delta|.
```

Completing the square gives

```
minimum in well a >= 2 Re(zeta exp(-iq_a))-(pi²/72)rho².
```

If x=Re(zeta), y=Im(zeta), the largest center value is
`max(2x,-x+sqrt(3)|y|)`. Its minimum at fixed |x|+|y|=gamma is
`(sqrt(3)-1)gamma`. Hence

```
mu3 >= (sqrt(3)-1)gamma-(pi²/72)gamma²
     > (71/100)gamma > (7/10)gamma.                         (10)
```

The weak rational bound 71/100 suffices; the displayed analytic bound is
stronger. Therefore for **every** chemical potential mu>=mu3, the exact
leading ground is a period-three crystal and its gap exceeds gamma/100.
This covers arbitrary phases of zeta in the certified diamond, including
chiral detunings intended to avoid commensurate locking. It resolves this
finite-filling escape without assuming positive compressibility,
relative stiffness or absence of clusters.

The scope is deliberately precise. Equation (10) compares with **bare
one-particle** well thresholds. It does not prove that interactions cannot
produce three dressed low-energy sectors at mu<mu3, nor does it classify
larger hopping, a!=0, or eta!=0. A general interacting mode count cannot
be deduced from (8). The higher-order microscopic terms also require a
uniform control estimate before (3) can be transferred to a fixed qubit
Hamiltonian. No general c<=2 claim follows.
