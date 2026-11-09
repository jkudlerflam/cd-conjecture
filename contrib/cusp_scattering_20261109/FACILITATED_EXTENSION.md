# A microscopic deformation that weakens the pair binding at small amplitude

AI assistance: Codex (OpenAI GPT-6 family), 2026-11-09. This is a
perturbative calculation on a newly proposed nearest-neighbour qubit
family. It is not a counterexample or a finite-density phase result.

## 1. Allowed microscopic deformation and its full leading model

Add to the original microscopic cusp Hamiltonian

\[
\lambda F,\qquad
F=\sum_j(n_jX_{j+1}-X_jn_{j+1}),\qquad n_j=(1-Z_j)/2.
\]

On a ring the one-site terms telescope, so this is exactly the allowed,
one-site invariant, nearest-neighbour bond
`(lambda/2)(X tensor Z-Z tensor X)`. It breaks the original parity.
Its compression to the blockaded space is zero, `PFP=0`.

[proved by local path enumeration; exact finite-matrix checks] Writing
`alpha=lambda/Delta` and `tau=Delta^2/(2V)`, the second-order contribution
from `F^2` is

\[
\frac{\lambda^2}{V}P\sum_j
\left[b^\dagger_{j+1}b_j+\mathrm{h.c.}-2n_j+2n_jn_{j+2}\right]P.
\]

Here `b_j` is the local spin-lowering operator, with no Jordan-Wigner
string. The mixed pairing-facilitation contribution is

\[
\frac{\Delta\lambda}{V}P\sum_j
(n_{j+2}-n_{j-2})X_jP. \tag{1}
\]

The vacuum-to-single-particle terms cancel on a ring, but (1) survives.
It cannot be discarded when increasing the repulsion. These formulas
were independently enumerated in this calculation on a nine-site ring
and verified entry by entry on rings 7 through 12 by the
[parallel reduction audit](../cusp_validation_20261109/FACILITATED_REPAIR.md).
That folder contains the reproducible exact matrix comparison.

Tune the bare real nearest hop by `-lambda^2/V` to cancel the generated
real hop, and absorb the diagonal number shift in the chemical potential.
In units `tau`, the complete leading deformation of C173's model is

\[
H(\alpha)=H(0)+
4\alpha^2\sum_j n_jn_{j+2}
+2\alpha\,P\sum_j(n_{j+2}-n_{j-2})X_jP. \tag{2}
\]

Projectors around the density term in (2) are implicit. Thus simply
inserting `v=4alpha^2` into the binding criterion would not describe this
microscopic deformation. A correct calculation must include the
particle-number-changing term.

## 2. The pair-to-trimer source

The following uses the original C173 gauge, `w=i`, `t2=-1`, and total
pair momentum `K=pi`. The exact two-body energy is
`E_b=-4.519526534666399`. Let `g_r` be the normalized real wavefunction
from [DERIVATION.md](DERIVATION.md), and let
`f_r=i^r g_r` be the leftmost-coordinate wavefunction in the original
gauge.

[proved, matrix elements] The term linear in `alpha` in (2) couples a
one-particle state at momentum `K` to the separation-two pair with
centered-coordinate amplitude `4i sin K`. Therefore the old bound pair
at `K=pi` has zero matrix element to the one-particle sector exactly.
It does couple to three particles.

For three ordered particles at positions `0,g,g+h`, `g,h>=2`, removing
each of the three particles gives the exact source for that coupling:

\[
B(g,h)=2\left[
\delta_{g,2}\{f_h-f_{g+h}\}
+\delta_{h,2}\{f_{g+h}-f_g\}
\right]. \tag{3}
\]

It is antisymmetric under `g <-> h`; the script checks this independently.
There is no continuum approximation in (3).

## 3. Pair energy shift at fixed chemical potential

[derived, perturbation theory with a spectral-isolation assumption]
Let `H3(K=pi)` be the number-conserving three-body Hamiltonian without a
chemical potential. At fixed `mu`, the unperturbed pair grand energy is
`E_b-2mu`. Provided the relevant resolvent exists and this branch can be
continued as an isolated eigenvalue, the coefficient of `alpha^2` is

\[
\Sigma_{\rm pair}(\mu)=4|f_2|^2
-\langle B|(H_3-\mu-E_b)^{-1}|B\rangle. \tag{4}
\]

The first term is the extra repulsion, and the second term includes the
entire three-body sector, including its interactions. Replacing that
sector by noninteracting atoms would not implement (2).

[numerical] `facilitated_pair_shift.py` uses C175's saved `nbody_box.py`
for `H3`, with independently constructed exact pair data and source (3).
Each gap is cut off at `R`, through `R=240`, a dimension of 57,121.
Linear-solve residuals are below `1e-10`. The direct term is
`4|f_2|^2=0.605277302659`.

| `mu` | Pair shift coefficient at `R=240` | Change from `R=160` |
| --- | ---: | ---: |
| -2.40 | -0.668887530300 | 2.22e-7 |
| -2.30 | -0.999614714156 | 3.22e-7 |
| -2.27 | -1.163523967762 | 3.80e-7 |
| -2.26 | -1.230870306648 | 4.05e-7 |

These coefficients lower the pair energy. To decide the binding change,
one must also lower the two-atom threshold; the direct term alone gives
the wrong comparison.

## 4. Both dressed valleys and the binding comparison

[numerical; independent cross-check] The coordinating agent computed
the atom coefficients with a separate two-body resolvent implementation,
saved in
[facilitated_atom_results.json](../cusp_binding_certificate_20261109/facilitated_atom_results.json).
This agent independently reproduced both coefficients using a banded
linear solve through 512 relative sites.

At `mu=-2.27`, the two valleys have

\[
k_a=-\arcsin(1/4),\quad k_b=\pi+\arcsin(1/4),
\]

\[
\Sigma_a=-0.825224871396,\qquad
\Sigma_b=-2.663400053460.
\]

The unequal shifts are consequential: the facilitated term breaks the
symmetry which exchanged the valleys. Comparing twice the first atom
shift to the pair shift would be incorrect.

An additional real nearest hop `c alpha^2`, in units `tau`, balances the
two valley energies through this order when

\[
c=-\frac{\Sigma_a-\Sigma_b}{4\cos k_a}
=-0.474614791170. \tag{5}
\]

Their common shift is `(Sigma_a+Sigma_b)/2`. The pair at `K=pi` has zero
linear response to a real nearest hop, so this retuning does not change
(4) at order `alpha^2`. At this fixed chemical potential the binding
therefore has the expansion

\[
\begin{split}
E_{\rm bind}(\alpha)
&=0.0195265346664+
(\Sigma_a+\Sigma_b-\Sigma_{\rm pair})\alpha^2+O(\alpha^4)\\
&=0.0195265346664-2.32510095709\,\alpha^2+O(\alpha^4).
\end{split} \tag{6}
\]

[heuristic, not a prediction with a controlled error] Setting the displayed
quadratic expression to zero gives `alpha approximately 0.09164`.
This is a promising scale for a direct calculation, but it is not an
established unbinding point. The omitted higher orders, valley-momentum
shifts, and branch isolation have not been bounded there. All coefficients
depend on the held-fixed chemical potential because particle number is
broken; they cannot be transported unchanged to finite density.

## 5. What this adds, and what remains

The exact second-order Hamiltonian supplies a microscopic way to move the
pair threshold using only allowed qubit nearest-neighbour terms. The
resolvent calculation includes its unavoidable charge-one term and finds
that the balanced-valley binding decreases at small amplitude. That is
a constructive lead beyond finding another artificial effective model.

The next decisive calculation is finite-amplitude spectroscopy of (2),
retuning the real hop to equalize the **dressed** valley minima and tracking
the pair pole against the dressed continuum. If binding is removed,
the remaining tests are absence of competing larger clusters, repulsive
low-energy scattering compatible with two components, stability against
number-changing resonances, and a comparison to the finite-coupling
microscopic qubit chain. No value of `c` follows from (6).

Reproduce the pair coefficient from the repository root:

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
python contrib/cusp_scattering_20261109/facilitated_pair_shift.py
```

The script takes approximately ten seconds and writes
[facilitated_pair_shift.json](facilitated_pair_shift.json).

## 6. Exploratory small-u limit

[numerical, exploratory] To ask whether reducing the original binding
can help, the derivative calculation was repeated at
`mu=epsilon_min-0.1u^4`. This is a different chemical-potential path from
the `mu=-2.27` calculation above. Both atom valleys were evaluated with
relative cutoff 8192, and the pair-to-trimer resolvent with increasing
cutoffs in both gaps.

| `u` | Largest gap cutoff | Binding derivative divided by `u^2` | Formal quadratic zero `alpha/u` |
| --- | ---: | ---: | ---: |
| 1 | 240 | -1.65707 | 0.10855 |
| 0.5 | 600 | -2.44526 | 0.10575 |
| 0.25 | 1000 | -2.66448 | 0.10643 |

The `u=0.5` pair coefficient changes by `4.44e-5` between cutoffs 360
and 600. At `u=0.25`, the remaining tail is estimated at about `4e-4`
from the known pair decay and the observed cutoff sequence; this last
row is **not** converged to the printed precision. Saved intermediate
boxes make that limitation visible. These scans suggest a useful path
`alpha` proportional to `u`, with smaller mixing weights as `u` falls.
They do not prove an asymptotic coefficient or a controlled phase.

There are two additional obstacles to such an inference. First, the
average atom shift is of order `alpha^2 u`, whereas the starting vacuum
gap is only of order `u^4`. Along `alpha` proportional to `u`, the chemical
potential must also be retuned by
`[(Sigma_a+Sigma_b)/2] alpha^2` to maintain the intended vacuum gap.
This retuning leaves the displayed binding derivative unchanged at
order `alpha^2`, but affects higher orders.

Second, a small admixture weight does not establish a small error on the
binding scale. Individual self-energies are of order `u^3` along this
path, while the original binding is of order `u^4`. An order-`u` relative
correction to an individual self-energy can therefore be as large as
the binding itself. One must prove cancellation between atom and pair
corrections or bound the remainder relative to the binding. Neither
has been done here. The small-u path is a lead, not a controlled
unbinding theorem.

Reproduce these scans with

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
python contrib/cusp_scattering_20261109/small_u_scaling.py
```

The sparse calculation takes about ten seconds on the development
machine and retains only the compact
[small_u_scaling.json](small_u_scaling.json), not the large vectors.
