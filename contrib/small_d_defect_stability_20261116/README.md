# Positive dressed defects in the frozen d=20 precursor

Date: 2026-11-16. Author: OpenAI Codex, GPT-6 family, with AI assistance
throughout. No independent agent or human review is claimed.

## Claim and status

**[Proved, awaiting independent review.]** For the literal frozen
20-state nearest-neighbour Hermitian bond of the
[new smaller-d family](../smaller_d_unbounded_20261116/README.md), a
sufficiently small **strictly nonzero** port coupling preserves the
selected length-n ground sectors on every ring divisible by n+1.
The bond is repeated with one-site translation symmetry and all couplings
are fixed independently of chain length. The result includes arbitrary
wrong-cell lengths, finite grammar penalty, domain boundaries,
marker-free sectors, and the small positive readout required by the
leading-band construction.

Write mu=-(n+2)E_n, gamma=E_(n+1), and
B=sqrt(2+2 eta)|lambda|. Under the parameter hypotheses in
[PROOF §1](PROOF.md), it suffices that

    gamma > 5B,             K > 2|mu| + 7B,             f >= 0.

If S is the exact all-n ring ground energy, a sector with N_W wrong cells
and N_I forbidden bonds has

    E_sector - S >= N_W (gamma-5B) + N_I (K-2|mu|-7B).

For illegal sectors, N_W counts the complete cells between forbidden
bonds; the remaining end fragments are already included in the bound.
The legal marker-free sector has energy >=0>S. Thus every competing
finite-ring ground sector is excluded, including one with a single
defect. This is a sector gap, not a gap within the interacting payload.

No chemical-potential correction or scalar counterterm is needed.
The bulk reference is the exact interacting selected-chain energy
density. The proof never substitutes a many-body perturbation error
proportional to the full volume.

## An explicit sufficient schedule

For every fixed n>=1 and 0<q<=1/4, set

    gamma_low=(1-q^2)^2 q^(2n+4),
    0<lambda<=gamma_low/[10 sqrt(2+2 eta)],   K=1.

This gives cost at least gamma_low/2 per wrong cell and cost greater
than 0.976 per forbidden bond. The readout can be precisely the parent's
positive field

    f=lambda^2 g^2 |v_0 v_n|^2 D^(-2n-2) nu/q^2,   nu>0.

All earlier parameters may be chosen in the parent's stated band order;
lambda can then be reduced further as needed. Constants and energy
scales may deteriorate arbitrarily with n. This schedule certifies
defect stability, not the parent's additional finite-band tolerances
or any microscopic conformal limit.

## Why the comparison closes

The selector and conserved projectors are rederived in the proof.
Positive readout leaves the all-zero-bit vacuum path unchanged, so every
isolated cell still has minimum mu+E_k at any f>=0. Ports have the
length-independent norm B. For an open run of j selected cells, exact
ground energies obey |A_j-j u|<=B, where u is their thermodynamic energy
per cell. This follows from cutting and joining one bounded port.

Cutting around N_W wrong cells costs at most 2B N_W; the surviving
selected runs cost at most another B N_W. A periodic selected reference
has at most a 2B boundary remainder. Because there is at least one defect,
the stronger gamma>5B condition absorbs this entire remainder. The
negative marker energy costs at most 2|mu| per illegal cut. These are
exact operator/variational estimates and require no global payload gap.

The method is inherited from
[the earlier route-A comparison](../counterexample_campaign/route_a/README.md)
and the [bulk-subtracted frozen analysis](../decisive_attack_20261028/frozen/DEFECT_PRESSURE.md).
The new work checks all its hypotheses for this 20-state bond, handles
its negative chemical potential, and supplies an all-n nonzero schedule.

## Checks and reproduction

**[Numerical corroboration, not interval certificates.]**
[checks.py](checks.py) independently reconstructs the entire bond from
the transition table and assembles matrices directly on physical letter
words. It imports the parent model only for an entrywise bond comparison.
The saved [checks.json](checks.json) records:

- Exact floating-point entrywise agreement with the parent bond at
  three diagnostic parameter points; individual site marker conservation,
  neighboring forbidden-projector conservation, and the exact port norm.
- Selector inequalities for k=0,...,18, n=1,...,12, and three q values,
  at 110 decimal digits. The proof, not this finite scan, covers the tail.
- Complete isolated-cell spectra through k=5 with nonzero positive
  field; grammar dimensions and minima agree with the formulas.
- Every conserved sector of full physical rings L=2 and L=3, respectively
  n=1 and n=2: all 400+8000 basis states and 80 sectors are included,
  including illegal and marker-free sectors.
- Thirteen mixed-length comparisons through L=6 and n=3, using the
  strictly nonzero certified port and the parent's field scaling. The
  largest selected space has dimension 21,952. All bounds pass.
- Joining open runs of one, two and three cells at a larger diagnostic
  lambda=0.08, where dressing is numerically resolved. This last check
  tests the boundary lemma only and is outside the defect schedule.

The very small readout and dressing shifts at the certified schedule
can lie below double precision. The checks validate the much larger
comparison margins; they do not measure those tiny shifts or a CFT.
The largest recorded eigenpair residual is below 1.6e-12. Sparse
eigensolver output is evidence, not a rigorous eigenvalue enclosure.

From the repository root, with Python, NumPy, SciPy and mpmath:

```sh
OPENBLAS_NUM_THREADS=1 python3 -B contrib/small_d_defect_stability_20261116/checks.py
```

The recorded run took about 22 seconds. It overwrites only the adjacent
small JSON result. The supplied `requirements.txt` lists dependencies.

## What this settles and leaves open

This closes the frozen d=20 nonzero-port ground-sector selection input
at the displayed sufficient schedule. The proposed c=n remains
conditional on the **full interacting** selected-chain infrared theory;
c denotes total nonchiral Virasoro central charge, counting one compact
boson as one. No thermodynamic central charge was measured here.

Compatibility with L divisible by n+1 is essential. The existing frozen
[incompatible-ring obstruction](../decisive_attack_20261028/frozen/ALL_RING_OBSTRUCTION.md)
persists, as do the translated, translation-breaking marker vacua.
These facts cannot be counted as extra conformal modes.

For the d=22 mobile model, fast-vacancy dynamical matching, positive
dressed defects at fixed nonzero J/t, and the interacting multichannel
CFT identification are still separate tasks. The
[later detailed route audit](../route_audit_20261109/README.md) explicitly
retains such model-specific hypotheses despite stronger older headlines
in STATUS.md. This theorem supplies one previously missing frozen input;
it does not turn the leading-band calculation into an unbounded-c
counterexample.

The top-level README and REPORT link this result and preserve these
qualifications. No external mathematical or physics premise was needed
beyond the local bond and the elementary estimates proved here.
