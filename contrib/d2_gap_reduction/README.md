# Proven qubit subclasses and a reduced staggered field problem

## Claim and scope

For the c–d conjecture at d=2, consider one fixed real exchange matrix J and
uniform real field B,

\[
H_L=\sum_j[\boldsymbol\sigma_j^T J\boldsymbol\sigma_{j+1}
 +B\cdot\boldsymbol\sigma_j],
\]

on the full periodic qubit tensor product. Pauli matrices have eigenvalues
±1. The conjectured bound is on total unitary Virasoro central charge
\(c_T\leq1\), with nonzero limiting mode velocities.

This contribution completely closes the **symmetric ferromagnetic rank-one
exchange** subclass, giving the stronger bound \(c_T\leq1/2\) for arbitrary
uniform fields. It also exhausts the microscopic additive-spin-symmetry
route and reduces its only additional interacting case to a two-parameter
staggered-field problem. The general d=2 conjecture remains open.

## Results and status

1. **[Proved using a primary-source gap theorem.]** For
   \(J=-\kappa nn^T\), \(\kappa>0\), \(|n|=1\), every nonzero longitudinal
   field gives \(\Delta_L\geq2|B\cdot n|\). With zero longitudinal field,
   this is exactly transverse-field Ising. Its only conformal locus is
   \(|B|=\kappa\), with \(c_T=1/2\). The proof includes the ordered
   doublet and decoupled-spin endpoints. See [rank_one.md](rank_one.md).
2. **[Proved using the same source.]** The larger criterion
   \[
   Jn=J^Tn=\lambda n,\quad
   \lambda\leq-\|(I-nn^T)J(I-nn^T)\|_{\rm op},\quad B\cdot n\ne0
   \]
   gives the same size-independent gap. Transverse fields are unrestricted;
   this excludes interacting points outside the previously controlled
   free, uniform-U(1), and zero-field-XYZ classes.
3. **[Proved algebraically; standard XXZ input for the CFT conclusion.]**
   Classify every exact charge \(Q=\sum_j n_j\cdot\sigma_j\), allowing
   arbitrary site dependence, amplitudes and zero coefficients. At B=0,
   all such Hamiltonians reduce to XXZ or commuting Ising, hence satisfy
   \(c_T\leq1\) with the established XXZ input. At B≠0, the only extra
   class is \(a(XX-YY)+uZZ+bZ\), equivalent to XXZ with a staggered field.
   See [additive_charge.md](additive_charge.md).
4. **[Proved exact spectrum on even rings.]** For that remaining class set
   \(r=\sqrt{a^2+b^2/4}\). If u>r, there are exactly two Néel vacua and
   the gap above them is \(4(u-r)\). At u=r, a≠0, the ground degeneracy
   is L+1 and a quadratic branch excludes an ordinary conformal limit.
   If u≤−|a| and b≠0, the gap theorem gives a unique ground state and
   gap at least 2|b|. See [staggered_slice.md](staggered_slice.md).

The literature input for results 1–2 is Bravyi, Gosset, Liu and Wong,
[arXiv:2607.14401v1, Definition 1 and Theorem 1](https://arxiv.org/html/2607.14401v1),
read directly. Their theorem, rather than an unproved phase-diagram
assumption, supplies the nonintegrable gap. Result 4 has an elementary
positive-operator and one-defect proof here; the corresponding exact line
and gapped side also appear in
[Zhang–Song, arXiv:2609.29728v1, §§III–IV](https://arxiv.org/html/2609.29728v1).
Their broader gapless phase identification uses finite-size DMRG.

## The precise next proof problem

A uniform rotation and positive rescaling set a=1, b>0 in the unresolved
staggered-charge class. The exact exclusions reduce its candidate conformal
region to

\[
\boxed{-1<u<\sqrt{1+b^2/4},\qquad b>0,\qquad u\ne0.}
\]

The line u=0 is quadratic; b=0 is XYZ/XXZ. The box is a **remaining search
region**, not a claim that every point in it is critical. An exhaustive
one-component infrared description inside it would finish the conjecture
for all qubit chains with an exact additive on-site continuous symmetry.
That would still leave chains whose currents are emergent, multispin, or
absent microscopically.

The boundary u=r is the factorization line already analyzed in
[C15 §§2.3–2.4](../../codex/work_C15/RESULT.md). Its gap vanishes quadratically
with momentum, so a relativistic c-theorem does not control its entire
neighborhood. This contribution proves the adjoining gapped wing directly;
it does not assume that all perturbations of the boundary have one boson.

Two further exact closure results prevent redundant searches:

- Arbitrary site-dependent product rotations followed by standard
  parity-preserving Jordan–Wigner cannot enlarge the known quadratic locus.
- A zero-field uniform bond can be rotated into a homogeneous XYZ bond in
  the bulk exactly when \(JJ^T=J^TJ\). This is already zero-field XYZ or
  axial XXZ plus DM. Inhomogeneous target bonds and entangling circuits
  are outside this statement.

## Reproduction and checks

Run from the repository root with Python 3 and NumPy:

```sh
python3 contrib/d2_gap_reduction/checks.py
```

The script writes [checks.json](checks.json). The recorded run took under
one second and passed:

- 15 complete finite-ring spectra at L=4,6,8 checking the source-based gap,
  including arbitrary transverse fields, a chiral exchange block, a global
  rotation, equality in the dominance condition, and an exactly saturated
  pure-field bound;
- 9 parameter/size checks of the exact staggered-family gap, with 9 boundary
  checks of the L+1 ground degeneracy;
- 14 additive-charge kernel cases and their rotated copies, including odd
  and even rings, proper and improper twists, and ranks zero through three;
- full spin-matrix commutators for all nullspace basis vectors, with maximum
  Frobenius residual below \(1.2\times10^{-13}\); eight independent random
  checks of the Pauli-coefficient norm identity agree within
  \(1.8\times10^{-15}\);
- critical-Ising spin matrices against the exact finite-ring ground energy.
  Its exact Casimir estimator is 0.5000022 at L=256, approaching 1/2.

These are convention and algebra checks. No numerical phase scan or new
numerical determination of an unknown central charge is claimed. The axis
checker uses floating-point tolerances; the displayed analytic conditions,
not its Boolean output near a tolerance threshold, certify the theorem.

## Provenance and relation to existing work

Written 2026-10-28 by Codex (GPT-6-based assistant) with three parallel
workers. Two workers independently derived and cross-checked the charge
classification and rotation closure. The coordinator and another worker
independently read the gap theorem and checked its signs and assumptions.
The coordinator ran the saved finite checks. No human review is claimed.

This extends [C6](../../codex/work_C6_d2/c6.tex) and the
[imported proof extensions](../chat_followup/proof_extensions.md).
The Ising solution, uniform-U(1) solution, and C15 factorization line are
inherited, not claimed as new. New to this repository are the complete
mixed-field ferromagnetic rank-one corollary, the larger rigorous gapped
region, the exhaustive additive-charge and rotation-closure arguments,
and their combined reduction to the displayed staggered-field region.
