# Mechanisms and the comparable-coupling continuation

OpenAI Codex, GPT-6 family, 2026-11-10. This is rationale for bounded searches,
not an exclusion theorem or a counterexample claim. The completed Stage A
motivates the separately registered follow-ups; this note does not analyze its
numerical files. The physical target remains one-site translation-invariant,
nearest-neighbour qubits. Longer-range effective Hamiltonians are diagnostic
tools, not substitutes for that target.

## What improving the existing mechanism actually costs

Use the conventions of [C190 models](../../codex/work_C190/models.py):

\[
\Delta=J_x-J_y=1,\quad q=J_x+J_y,\quad
r=V/\Delta=4J_z,\quad w=q+2iD,\quad |t_2|=\frac1{2r}.
\]

On the symmetric line, put \(u=|w|/|t_2|\). Then

\[
J_x=\tfrac12,\quad J_y=-\tfrac12,\quad
J_z=\frac r4,\quad D=\frac{u}{4r},\quad
h_z=-\frac r2+\frac{6+\mu_{\rm eff}}{4r}.
\]

The last formula is a **leading-order effective-model field prescription**,
not a reliable density prediction at small r. Fields must be retargeted in the
full microscopic chain. At u=1, r>=2 and mu_eff approximately -2.17, the
exchange coefficient spread is r^2 and the total spread is approximately
2r^2-3.83. Thus r=8 still gives total spread about124, r=4 about28, and r=2
about4.17. An overall rescaling changes none of these ratios.

Reaching the preferred spread of at most3 therefore requires leaving the
controlled large-r blockade expansion, increasing u, or changing mechanism.
The choice r=2,u=3/2 is a useful rational seed:

\[
H=\sum_j\left[\frac{X_jX_{j+1}-Y_jY_{j+1}+Z_jZ_{j+1}}2
 +\frac3{16}(X_jY_{j+1}-Y_jX_{j+1})+h_zZ_j\right].
\]

At the simple starting field h_z=-1/2 its coefficient spread is8/3, with six
expanded Pauli monomials or four grouped terms. Changing D to1/4 gives spread2.
Neither is asserted to be critical. Their fields and their finite-size behavior
require direct measurement. Increasing u is not cost-free:
[C182](../../codex/work_C182/RESULT.md) finds pairing/clustering and generally
c approximately1 homogeneous liquids at effective-model u approximately2–3.
Those findings motivate including intermediate u=3/2; they do not exclude the
strongly dressed microscopic r=2 chain.

## Invariant scales and cancellations

Write a bond exchange as \(\sum_{ab}J_{ab}\sigma_j^a\sigma_{j+1}^b\), with

\[
J=\begin{pmatrix}(1+q)/2&D&0\\-D&(q-1)/2&0\\0&0&r/4\end{pmatrix}.
\]

Its symmetric part has eigenvalues {(1+q)/2,(q-1)/2,r/4}; the DM-vector
magnitude is |D|, its antisymmetric Frobenius norm is sqrt(2)|D|, and the field
norm is |h_z|. A common on-site spin rotation preserves these quantities.
The interaction singular values are

\[
\left\{\frac{|1+\sqrt{q^2+4D^2}|}{2},
\frac{|1-\sqrt{q^2+4D^2}|}{2},\frac{|r|}{4}\right\}.
\]

At the first r=2 seed they are {11/16,5/16,1/2}, giving spread11/5.
At D=1/4,q=0 they are {3/4,1/4,1/2}, giving spread3. Singular values alone
understate the old candidate's very weak DM: its symmetric and antisymmetric
parts must also be reported. Exact zero scales are listed separately; a nearly
zero singular value is an additional hierarchy to disclose.

Track |q|/(|Jx|+|Jy|) and |4Jz+2hz|/(4|Jz|+2|hz|) explicitly. An exact
symmetry line q=0 and ordinary adjustment of a field to criticality are
different from relying on a tiny nonzero cancellation. Comparable coefficients
do not automatically imply comparable generated dynamics or a broad critical
window. A tuning resolution is not a measured tolerance: report both the step
used and any actual parameter interval over which behavior survives.

## Microscopic charge locking: a warning with substantial uncertainty

[C191, section8](../../codex/work_C191/PROOF.md) discusses the symmetry-allowed
microscopic two-pair operator
\(O_4=(\Psi_a\Psi_b)^2\propto e^{2i\sqrt2\theta_c}\), absent from the
number-conserving second-order effective Hamiltonian. In that treatment it is
resonant on q=0 and has dimension x4=2/Kc approximately1.45–1.55. Combining
this with the estimate y4 approximately c4(16/9)r^-2 gives a locking length of
order10^6 at r=40 and order10^3 at r=8, assuming an order-one unknown c4.

These are **heuristic length estimates, not established microscopic gaps**.
The nonuniversal coefficient c4 was not computed there. Kc is inferred using
the effective model and approximate scattering/continuum matching; the charge
parameters and velocities of the full chain need not equal those inputs.
The power-law estimate assumes a weak perturbation of the stated Gaussian
fixed point. It cannot be quantitatively extrapolated to r approximately2,
where the blockade expansion itself is uncontrolled. Nor does symmetry alone
prove a nonzero coefficient at every microscopic point. C191 explicitly
limits its microscopic-chain analysis to this estimate.

Nevertheless, a finite-size entropy near2 on q=0 does not by itself rule out
eventual charge locking. Reducing r increases this estimated bare perturbation
in t2 units. Nonzero q can detune the resonance by shifting the two valley
momenta, but whether two occupied valleys survive that shift is a separate
microscopic question. This is a reason to compare both lines, not a proof
that either has c>1 or that the symmetric line must be gapped.

## A bounded off-symmetric rational follow-up

Proposed next grid, **not registered as already run**: keep r=2, Delta=1,
take D in {3/16,1/4}, and q in {0,+/-1/8,+/-1/4,+/-1/2}. This is14 exchange
bonds, with q=0 controls. For an initial bounded field bracket use
h_z in {-3/4,-5/8,-1/2,-3/8,-1/4}, then refine only a bracket whose microscopic
diagnostics justify it. This is70 points before symmetry reductions. Record
field steps as1/8; any later refinement must have a new registration.

Here Jx=(1+q)/2 and Jy=(q-1)/2, so |Jx|+|Jy|=1 and the nonzero cancellation
fractions are1/8,1/4,1/2. These are explicit moderate offsets, not hidden
arbitrarily small differences. At D=1/4 the displayed coefficient spread is
at most3 throughout this field bracket. At D=3/16 it is at most4: classify each
point separately as preferred or exploratory. Also report singular values:
for example q=1/2,D=1/4 has singular-value spread about5.83 despite raw spread3.

This grid can test whether a comparable-coupling two-component region has
appreciable width in q. It is not a sensitive search for an extremely narrow
partially polarized window. In the leading-order band, at fixed D a small q
splits the formerly symmetric valley minima by approximately
4|q|sqrt(1-u0^2/16), where u0=4rD. At r=2, q=1/8 this is already about1.7–1.9
times |t2| for these D values. That can empty one valley before it usefully
detunes charge locking at low density. This estimate is qualitative at r=2;
measure occupations/correlations rather than vetoing points from it. A null
result here would mean that these moderate rational offsets did not help,
not that all smaller offsets or all comparable-coupling mechanisms fail.

## Different mechanisms and the microscopic range test

The search should admit any c>1, including c=3/2 from an Ising mode coexisting
with a liquid. Two primary examples motivate that target but do not satisfy
the requested microscopic range:

- [Rahmani et al., interacting Majorana chain](https://arxiv.org/abs/1505.03966)
  finds an Ising-plus-Luttinger regime. The explicit spin representation in
  [Karcher et al., Eq.(17)](https://arxiv.org/pdf/1905.12743) is
  \(\sum_j[t_1X_j-t_2Z_jZ_{j+1}-g_1X_jX_{j+1}-g_2Z_jZ_{j+2}]\).
  Its last term is a microscopic next-nearest-neighbour exchange. Counting
  Majoranas as physical qubit sites would also misstate the local Hilbert space.
- [Gotta et al., pairing in spin chains](https://arxiv.org/abs/2006.07128)
  finds c=3/2 transitions between ordinary and cluster liquids in models with
  nearest- and next-nearest-neighbour interactions. It is a mechanism example,
  not an admissible nearest-neighbour qubit counterexample.

[Yi et al.](https://arxiv.org/abs/1907.05303) studies the actual XYZ+DM+field
family and its antiferromagnetic, paramagnetic and chiral phases; it supplies
no c>1 construction. The local
[multicritical search](../d2_multicritical/README.md) found tested band-edge
Ising/floating junctions separated by a gapped strip or ending at
nonconformal Lifshitz points. This is a warning against assuming that an
Ising line and a liquid simply add their central charges, not a complete
classification of strong symmetry-breaking or isolated multicritical bonds.
Its free chiral c=1 control also produced ring entropy fits up to2.4 at
L<=20. Thus a fresh mechanism must survive controls, sector selection,
size/bond-dimension checks, and tests for order, shell effects, clustering and
nonlinear dispersion before being interpreted as a conformal c>1 state.
