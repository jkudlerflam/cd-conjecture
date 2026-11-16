# A revised central-charge conjecture for one-site translation-invariant qubit chains

Date: 2026-11-16. AI assistance: OpenAI Codex, GPT-6 family. This is a new
working conjecture and a set of supporting calculations, not a proof of the
universal bound or an independent confirmation of the existing c≈2 candidate.
No independent agent or human review was performed.

## Claim and status

**[Heuristic / conjecture.]** Let

\[
 H_L=\sum_{j=1}^{L}h_{j,j+1},\qquad
 h=h^\dagger\in\operatorname{End}(\mathbb C^2\otimes\mathbb C^2),
 \qquad L+1\equiv1.
\]

The same fixed, finite bond acts on every link, independently of L. In
particular, **translation by one physical lattice site is required**.
Uniform on-site fields are allowed and can be split between adjacent bonds.
There is no imposed spin symmetry, reflection symmetry, time-reversal
symmetry, integrability, or restriction on coupling ratios.

If the complete low-energy bulk theory is a compact unitary conformal field
theory, or a finite collection of such sectors with finite nonzero limiting
velocities, then

\[
 \boxed{c_{\rm tot}=\sum_a c_a\leq 2.}
\]

Here c_a is the usual nonchiral Virasoro central charge: a compact boson
counts as 1 and an Ising CFT as 1/2. It is not c_L+c_R. Different sectors may
have different velocities. The statement includes isolated multicritical
points as well as extended critical phases. With spontaneous symmetry
breaking, choose a bulk ground-state sector; do not add central charges of
distinct vacua as if they were independent gapless fields.

**Separate sharpness conjecture:** c_max(2)=2, with a fixed finite-coupling
chain attaining c=2. The blockade candidate is the proposed realization.
Its finite-size data support sharpness but do not establish it.

Take L→∞ at fixed couplings before any singular-coupling limit. Energetically
generated constraints at finite coupling are allowed. A projected constrained
model alone, an L-dependent bond, a longer-range Hamiltonian, or an explicit
multi-site unit cell is not an example in this class. A nonconformal gapless
point with z≠1 does not acquire a c merely by fitting its entropy.

**Updated assessment (2026-11-16):** c≤2 remains a tentative upper-bound
guess, with substantially weaker support than the c=2 attainability
conjecture. The structural checks below concern particular mechanisms;
they do not establish that F(2) is finite. A sequence of fixed,
one-site-invariant qubit bonds with c tending to infinity remains open,
including sequences with increasingly small positive velocities and long
crossover scales. Confirming the candidate would establish a lower bound,
not this upper bound. The [small-d continuation](../small_d_lifting_20261116/README.md)
therefore develops constructive attainability conjectures and keeps all
proposed maxima separate.

## Evidence as of November 16

The live GitHub evidence used here is STATUS.md dated November 12 and
REPORT §9.205 at commit
[33248f8482dbc4c0b0f8b6782e5b156c47c48cc0](https://github.com/jkudlerflam/cd-conjecture/commit/33248f8482dbc4c0b0f8b6782e5b156c47c48cc0),
checked on November 16. The local checkout predates these results; this note
does not silently treat its older STATUS.md as current.

**[Numerical, existing work; not rerun here.]** The two leading microscopic
bonds read c=2.00 from 128→256 sites after bond-dimension convergence. From
256→512 they read 1.95–1.98 at χ=2048; a c=2 control reads approximately 1.95
at these sizes. The microscopic 512-site entropy increment still slightly
exceeds the registered convergence tolerance. Infinite-chain convergence
of the microscopic model remains unresolved.

The effective model reads c≈2 on the infinite chain. Its latest finite-chain
sequence is 1.82, 1.88, 1.93 over 64→128→256→512. Its 512→1024 value is 1.60
at χ=2048, with an entropy increment of 0.026; this fails convergence and is
not evidence for an infrared c=1.60. These are three distinct pieces of
evidence: microscopic finite chains, effective finite chains, and effective
infinite chains.

A simple microscopic representative, with Pauli normalization, is

\[
 H=\sum_j\left[\tfrac12(X_jX_{j+1}-Y_jY_{j+1})
 +\tfrac1{160}(X_jY_{j+1}-Y_jX_{j+1})+10Z_jZ_{j+1}\right]
 -19.976047\sum_j Z_j.
\]

Existing broad searches, including about 103,000 sampled bonds, found no
other convincing high-c region. This is weak evidence against narrow or
fine-tuned exceptions: those searches initially missed the blockade corner
itself. Earlier claims of a universal single-species description cannot be
used to prove the revised conjecture.

## A firm restriction: exact microscopic U(1) still gives c≤1

**[Proved algebraic reduction; established XXZ physics for the CFT value.]**
Suppose the chain conserves a nontrivial uniform additive spin component
Q=Σ_j n·σ_j. A common on-site rotation makes Q=Σ_j Z_j. Orthogonality of
Pauli strings in [H,Q] then gives

\[
 H=\sum_j[a(XX+YY)+b(XY-YX)+uZZ+hZ]_j+\text{constant}.
\]

A site-dependent rotation about z removes the hopping phase, giving XXZ
with a longitudinal field and a boundary twist. In the thermodynamic bulk
its conformal sectors have c≤1. Saturation/ferromagnetic endpoints with
quadratic dispersion are excluded by the CFT hypothesis. If a=b=0 the
Hamiltonian is classical and has no additional conformal liquid.

The site-dependent rotation is a solution method; the original Hamiltonian
still has exactly one-site translation symmetry. An *emergent* U(1) at low
energies does not satisfy this microscopic hypothesis. Thus the c≈2
candidate is consistent with this result.

As an algebra check, `check_structure.py` takes the commutator of the full
12-parameter uniform qubit Hamiltonian with ΣZ on rings of four and five
sites. Its kernel is exactly the four-dimensional span above, to less than
5e-16 in the kernel projector. This numerical check supplements the Pauli
algebra argument; it does not prove an infinite-chain spectrum.

## Why two is a natural ceiling in the known mechanism

Write n=(1−Z)/2, b†=|1⟩⟨0| and start from VΣn_j n_(j+1), V>0. Its low-energy
space P has no adjacent occupied sites. Consider adjacent-pair creation
amplitude p and the first-order-invisible transverse term

\[
 \sum_j\big[g(n_{j-1}-n_{j+1})b_j^\dagger+\mathrm{h.c.}\big].
\]

**[Derived by second-order perturbation theory; finite-ring matrix checks.]**
The isolated-particle hopping coefficients induced by these terms are

\[
 t_1^{\rm ind}=|g|^2/V,\qquad
 t_2^{\rm ind}=-|p|^2/(2V),\qquad t_{r\geq3}^{(2)}=0.
\]

The bare nearest-neighbour hop adds freely to t_1. Mixed p,g terms also
change particle number; they must be retained in the many-body effective
Hamiltonian. The displayed one-particle compression is not itself an
invariant sector of that full Hamiltonian.

The second-neighbour amplitude has one squared-magnitude source: transverse
facilitation cannot cancel it at this order. With p≠0, this gives a
dispersion of the form

\[
 \epsilon(k)=\epsilon_0+2\operatorname{Re}(t_1e^{ik}+t_2e^{2ik}).
\]

**[Proved band lemma.]** A nonconstant real trigonometric polynomial of
degree two has at most two nondegenerate minima. Multiply ε′(k) by z²,
z=e^(ik), to obtain a polynomial of degree at most four. Its simple roots
on the circle alternate between minima and maxima. At a regular Fermi
level the same polynomial-degree argument gives at most two occupied
pockets. Away from degenerate stationary points, sufficiently small smooth
corrections preserve that count.

Within the ordinary dilute-gas/two-component Luttinger description, those
two pockets supply at most two bosons. Locking a density or phase
combination removes a mode; it does not add an oscillator to that theory.
This explains the proposed c=2 phase and makes two a natural candidate
ceiling in its controlled regime.

**Scope:** the band lemma is not a theorem about interacting central charge.
It does not exclude independent composite fluids, an additional neutral
critical sector, or an irreducible strongly interacting CFT. Nor does
microscopic range one imply that the exact effective Hamiltonian stops at
range two. These are precisely the gaps between this argument and a
universal bound. A common-velocity Lorentz-invariant c-theorem is not being
assumed for all lattice flows.

`check_structure.py` independently constructs the complete second-order
matrix on P for 60 complex-amplitude cases on rings L=8,9,10. The displayed
coefficients agree within 7.2e-15. It also counts stationary points for
10,000 range-one and 10,000 range-two random bands. None exceeds its range
in number of minima; a deliberately out-of-scope range-three control has
three. This tests the calculation, not 20,000 interacting quantum chains.

## New attack: can repeated tuning manufacture a larger effective range?

Set p=0 and define the uniform nearest-neighbour operators

\[
 F=\sum_j(n_jX_{j+1}-X_jn_{j+1}),\quad
 T=\sum_j(b^\dagger_{j+1}b_j+\mathrm{h.c.}),\quad N=\sum_j n_j.
\]

Study, with fixed small ε,

\[
 H(\epsilon)=\sum_jn_jn_{j+1}+\epsilon F+
 \epsilon^2(-T+2N).
\]

All three terms have the same bond at every site. Although PFP=0, the
second-order virtual term is P(T−2N+2Σn_j n_(j+2))P. Consequently the
chosen counterterms cancel the leading motion and leave

\[
 H_{\rm eff}^{(2)}=2\epsilon^2P\sum_jn_jn_{j+2}P.
\]

**[Exact second-order identity.]** There really is a second blockade:
the next low-energy space P2 excludes separations one and two. Therefore
an argument that dismisses longer effective constraints just because the
microscopic Hamiltonian is nearest-neighbour is invalid.

**[Numerical diagnostic, not a c measurement.]** `attack_nested.py` keeps
all particle-number sectors and extracts the low-energy Hamiltonian on P2
by a polar-decomposition/direct-rotation construction. It checks rings
L=8,9,10 at ε=0.1,0.075,0.05,0.035. The second-order identity is exact in
these matrices and the P2 compression of the third-order term vanishes.
The fourth-order one-particle coefficients approach an on-site term 4ε⁴
and a nearest-neighbour hop −2ε⁴. Longer one-particle hops decrease relative
to ε⁴. At L=10, ε=0.035, the scaled hops at distances 1,2,3 are
−1.98965, −0.001216, −0.001411. This particular path did not expose a
three-valley leading band. Interactions and correlated processes remain;
no central charge is inferred from these small rings.

### An exact continuation of the cancellation is not conformal

Consider instead the finite-coupling family

\[
 H=V\sum_j n_jn_{j+1}+\lambda F+wT-\mu N,
 \quad w=\frac{V-\sqrt{V^2+8\lambda^2}}4,\quad \mu=2w,
 \quad V>0,\ \lambda>0.
\]

It agrees with the preceding tuning through order λ²/V, and its next
counterterms cancel the measured isolated-particle motion at fourth order:
w=−λ²/V+2λ⁴/V³+O(λ⁶/V⁵).

**[Proved reduction, using a published ground-space theorem.]** Put
a=√(−w), b=λ/a. The literal bond, including half of the field at each end,
factorizes exactly:

\[
 h=|\psi\rangle\langle\psi|,\qquad
 |\psi\rangle=-a|01\rangle+a|10\rangle+b|11\rangle,
 \qquad b^2=V-2w.
\]

The transfer matrix of Bravyi–Gosset is

\[
 T_\psi=\begin{pmatrix}-a&b\\0&-a\end{pmatrix}.
\]

It is invertible and no positive power is proportional to the identity.
Their Proposition 3 therefore gives exactly two periodic-chain ground
states. They are visibly the vacuum and the uniform one-particle W state:

\[
 |0\rangle^{\otimes L},\qquad
 |W_L\rangle=L^{-1/2}\sum_j b_j^\dagger|0\rangle^{\otimes L}.
\]

Every linear combination has Schmidt rank at most two across any spatial
cut, hence S≤log 2 independently of L. This exact family does not have a
positive-c conformal ground-state limit. The increasingly small gaps of
this family must not be interpreted as many conformal modes.

The factorization was checked for 200 parameter choices, with relative
error below 6.3e-16. Rings L=4,…,10 have exactly the two predicted ground
states. The finite-ring checks are supplementary; the all-L ground-space
statement uses the published theorem and the displayed transfer matrix.

This maps one natural continuation of the nested cancellation to a known
frustration-free class. It does **not** show that every such continuation
must follow this curve, nor that detuning it cannot produce a high-c phase.
An all-orders exclusion of hierarchical constructions remains open.

### Following the cancellation one stage further

**[Numerical extension.]** Add an adjacent-pair term αε² Pair, where
Pair=Σ(b†_j b†_(j+1)+h.c.), to the preceding H(ε). It leaves the second-order
blockade unchanged. The 24 complete-matrix calculations at α=0.5,1,2,
L=9,10, and four ε values again show vanishing next-neighbour hopping
relative to ε⁴ as ε decreases. It is essential to eliminate the states
violating the *second* blockade, not just to inspect the first effective
Hamiltonian's next-neighbour hop.

The extension also exposes another exceptional tuning: α=−1 suppresses
the fourth-order nearest-neighbour hop. Four further L=10 calculations
record this, rather than discarding it as an outlier. With an additional
chemical-potential adjustment, the literal microscopic family becomes

\[
 H=D_1+\epsilon F+\epsilon^2(-T+2N-\mathrm{Pair})-2\epsilon^4 N.
\]

`attack_deeper.py` extracts its lowest band on states excluding separations
one through four, on L=10,11,12 at three ε values. Band gaps and projection
overlaps support this finite-ring extraction, but particle-number mixing
remains appreciable at the tested ε. The long-distance hopping matrix
elements decrease relative to ε⁶ as ε decreases. **These data do not
establish a closed single-particle model or its central charge.**

There is a second exact continuation, now including pairing:

\[
 w=p=-\lambda^2/V,\qquad
 \mu=V-\sqrt{V^2+4\lambda^2}
 =-2\lambda^2/V+2\lambda^4/V^3-4\lambda^6/V^5+\cdots.
\]

**[Proved elementary factorization.]** Parameterize it by 0<r<1, with
V=(1−r²)², λ=r(1−r²), w=p=−r², μ=−2r²(1−r²). Up to an overall scale this
covers arbitrary λ/V>0. The bond obeys

\[
 h+r^4 I+r^3(I\otimes X-X\otimes I)
   =|\psi_r\rangle\langle\psi_r|,\qquad
 |\psi_r\rangle=(-r|0\rangle+|1\rangle)
                 \otimes(r|0\rangle+|1\rangle).
\]

The added single-site difference telescopes on the periodic chain and
does not alter its one-site translation symmetry. The remaining shift is
the scalar Lr⁴. Let |a⟩∝|0⟩+r|1⟩ and |b⟩∝|0⟩−r|1⟩. In this invertible
local basis, the projector constraint forbids the ordered pair ba and
allows aa, ab and bb. A periodic binary string without ba must be all a
or all b. Thus the entire periodic ground space is exactly

\[
 \operatorname{span}\{|a\rangle^{\otimes L},|b\rangle^{\otimes L}\}.
\]

This reasoning applies to arbitrary superpositions because the local
basis is invertible and the constraint kills each forbidden pair's
coefficient separately. Every ground state again has Schmidt rank at
most two. This second exact continuation is therefore also not a
positive-c conformal example. The identity was checked at 200 r values
to 2.3e-16, and the predicted full ground space on rings L=4,…,10.

The two exact continuations strengthen the case against these particular
attempts to obtain more modes by repeated cancellation. Neither rules out
all detunings, other singular limits, or additional strongly interacting
sectors. This distinction is part of the conjecture's remaining risk.

## Adversarial assessment and what would falsify the conjecture

| Possible exception | Assessment |
| --- | --- |
| Three independent microscopic subchains or a higher-range free chain | Requires interactions beyond the stipulated nearest-neighbour qubit bond. Arbitrarily high c at longer range is compatible with this conjecture. |
| More valleys in the ordinary blockade reduction | Excluded by the leading band lemma when p≠0, within its regular perturbative regime. |
| Repeated cancellation of the leading dynamics | A genuine loophole. Calculations reveal further constraints and two exact continuations with bounded-entanglement ground spaces. Other detunings remain open. |
| A third composite fluid | Not excluded universally. Earlier bound-pair searches do not constitute an exclusion of all higher composites. |
| An extra Ising sector on top of the c=2 liquid | Would give c=5/2 and refute the bound. The two-boson description cannot exclude a separate microscopic sector becoming soft. |
| An irreducible interacting CFT with c>2 | Not excluded by counting valleys. For example, an emergent SU(2)_5 theory has c=15/7; microscopic on-site SU(2) is not required for emergence. |
| A small-system entropy estimate above two | Insufficient: the repository already contains severe finite-size and unconverged overshoots, including in c=2 controls. |

Old potential exceptions such as Ising×boson (3/2), Potts×boson (9/5), or
Z_k parafermions (2(k−1)/(k+2)<2) do not refute this revised inequality.
This is a benefit of the conjecture's scope, not evidence that those
universality classes have been realized by an in-scope bond.

A serious falsification attempt should keep a single fixed bond, show a
stable estimate strictly above two with controlled size and bond-dimension
errors, and independently identify its conformal modes or spectrum.
Multiple velocities must be resolved rather than fitting Σc_a v_a as though
it were a single c. Exact constructions would be stronger. A finite scan
cannot rule out the isolated points included in the conjecture.

The strongest remaining target is an additional neutral/composite sector
near a singularly tuned bond, rather than merely increasing a small-ring
entropy estimate. No exhaustive classification of those sectors is claimed.

## Reproducing the new calculations

Python with numpy is sufficient; the recorded version is numpy 2.3.5.

```bash
python3 contrib/d2_updated_conjecture_20261116/check_structure.py
python3 contrib/d2_updated_conjecture_20261116/attack_nested.py
python3 contrib/d2_updated_conjecture_20261116/attack_nested.py --sizes 9 10 --eps 0.05 0.035 0.025 0.0175 --pair 0.5 1 2 --output contrib/d2_updated_conjecture_20261116/nested_pair_results.json
python3 contrib/d2_updated_conjecture_20261116/attack_nested.py --sizes 10 --eps 0.075 0.05 0.035 0.025 --pair -1 --output contrib/d2_updated_conjecture_20261116/nested_cancel_results.json
python3 contrib/d2_updated_conjecture_20261116/attack_deeper.py
python3 contrib/d2_updated_conjecture_20261116/check_product_parent.py
```

The scripts write the correspondingly named `*_results.json` files beside
themselves. They retain all parameter-grid outcomes and error measures,
with seeds for randomized checks and timing for the numerical extractions.
The numerical library was limited to two threads. No DMRG run, remote compute job, or
independent human verification is implied by these calculations.

## Primary sources and relation to existing work

- [Latorre–Sierra, 2403.17242v2](https://arxiv.org/html/2403.17242v2):
  original nearest-neighbour, translation-invariant conjecture and
  established examples. This note replaces its d=2 upper bound provisionally.
- [Verresen–Vishwanath–Pollmann, 1903.09179](https://arxiv.org/abs/1903.09179):
  constraints can suppress pairing and yield an emergent U(1) liquid.
- [Chepiga–Mila, 2206.11754](https://arxiv.org/html/2206.11754):
  nearest-neighbour interacting Kitaev/XYZ criticality; its real-hopping
  family is not a survey of all complex chiral bonds.
- [Mishra–Greschner–Santos, 1411.7270](https://arxiv.org/abs/1411.7270):
  two-component and paired liquids in frustrated triangular ladders,
  relevant to the effective mechanism rather than an in-scope microscopic bond.
- [Verresen–Jones–Pollmann, 1709.03508](https://arxiv.org/html/1709.03508)
  and [Lahtinen–Ardonne, 1504.07044](https://arxiv.org/abs/1504.07044):
  free critical chains and higher-c cluster models; longer spin support
  cannot be omitted when applying a local-dimension conjecture.
- [Bravyi–Gosset, 1503.04035, §2.2, Proposition 3](https://arxiv.org/html/1503.04035):
  the periodic rank-one frustration-free ground-space result used above.
- [Masaoka–Soejima–Watanabe, 2406.06414v2](https://arxiv.org/html/2406.06414v2):
  quadratic-or-softer gapless excitations in nearest-neighbour frustration-free
  qubit systems, consistent with excluding the exact parent from a wholly
  conformal low-energy description.

Repository inputs include C169 (structured searches), C173/C180 (two-valley
mechanism), C185/C186 (general cusp perturbations), the
[complete facilitated reduction](../cusp_validation_20261109/FACILITATED_REPAIR.md),
and live REPORT §§9.192–9.205. They were assessed as research records with
their stated limitations, not as a completed proof of the global claim.
