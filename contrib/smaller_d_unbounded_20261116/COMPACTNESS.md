# Compactness and the missing local central charge bound

Date: 2026-11-16. AI assistance: OpenAI Codex, GPT-6 family. This records
the compactness argument considered before the smaller-d construction.
The topological reduction and entropy inequality are proved below; the
local dynamical hypothesis needed to finish an upper bound remains open.
No independent agent or human review has been performed.

## Scope and normalization

Consider a single finite Hermitian nearest-neighbour bond h on
C^d⊗C^d, repeated at every physical site independently of chain length.
The central charge is the sum of the usual nonchiral Virasoro charges of
the complete unitary conformal bulk theory, allowing a finite collection
of sectors with finite nonzero velocities. Distinct symmetry-breaking
vacua do not contribute additive central charges.

Remove additive constants and telescoping terms A⊗I−I⊗A, and fix a
positive overall energy scale. For qubits a representative is

\[
 H(J,b)=\sum_j\left[\sum_{a,b=x,y,z}J_{ab}\sigma_j^a\sigma_{j+1}^b
 +\mathbf b\cdot\boldsymbol\sigma_j\right],
 \qquad \|J\|_F^2+|\mathbf b|^2=1.
\]

**[Proved.]** This normalized space K_2 is compact. For general finite d,
the same construction uses a finite-dimensional real space of physical
bond coefficients and a unit sphere, giving a compact K_d. Positive
rescaling preserves the ground-state problem and c; negative rescaling
is not identified. Quotienting by common onsite rotations is optional.

Let C_d⊂K_d denote the bonds with the stipulated conformal infrared limit.
There is no assumption that C_d is closed or that c is continuous on it.

## The exact finite cover reduction

**[Proved.]** The following two statements are equivalent:

1. There is B(d)<∞ such that c(h)≤B(d) for every h∈C_d.
2. Every h_*∈K_d has an open neighborhood U_* and a finite B_* such
   that c(h)≤B_* for every h∈U_*∩C_d.

The first statement immediately implies the second. Conversely, the U_*
cover K_d. A finite subcover U_1,...,U_N gives

\[
 c(h)\le\max_{1\le i\le N}B_i<\infty.
\]

Equivalently, if c(h_n)→∞, a subsequence converges to a bond h_* at which
c is unbounded in every neighborhood. This limiting bond may be
nonconformal, so proving local bounds only at conformal points is
insufficient. A finite-valued upper-semicontinuous extension of c to K_d
would suffice, but such an extension has not been established.

This reduction does not give a numerical bound such as c≤2. It permits an
isolated finite-c multicritical point of any particular central charge,
provided it has a locally bounded neighborhood. Consequently, establishing
some finite B(2) is weaker than identifying c_max(2).

## A sufficient entropy condition

Suppose a bulk ground state has the usual conformal interval entropy
S_h(ell)=(c(h)/3)log ell+s_h+o(1). Define

\[
 D_h(\ell)=\frac{3}{\log2}\,[S_h(2\ell)-S_h(\ell)].
\]

Then D_h(ell)→c(h). If a chosen symmetry-breaking sector has endpoint
oscillations, use a consistent translation average of interval entropies;
the entropy inequality below survives that average. The conformal entropy
law identifies the intended infrared observable and is not itself inferred
from finite local dimension.

**[Proved.]** Subadditivity for the two adjoining ell-site blocks gives

\[
 S_h(2\ell)-S_h(\ell)\le \ell\log d,
 \qquad D_h(\ell)\le3\ell\log_2d.
\]

**[Additional hypothesis; open in the microscopic class.]** Near every
h_* there is a finite integer ell_*, a finite A_*≥0, and a neighborhood
U_* such that all h∈U_*∩C_d satisfy

\[
 c(h)\le D_h(\ell_*)+A_*.
\]

This would imply the local bound B_*=3ell_*log_2d+A_*, and compactness
would finish the proof. For qubits this is 3ell_*+A_*.

The hypothesis requires a common finite-scale upper estimate for eventual
c in a neighborhood. It does not require all those models to have already
reached their final CFT at that scale, nor is a uniformly bounded crossover
length a necessary condition for bounded c. Nevertheless, this one-sided
control is substantive and is not supplied by compactness or subadditivity.

The entropic c-theorem provides decreasing entropy differences in unitary
Lorentz-invariant continuum theories. Applying it here would require
justified matching of a neighborhood of lattice Hamiltonians to a common
continuum description, including every low-energy sector. The theorem
does not automatically hold from the lattice spacing onward, nor for an
arbitrary coupled flow with different velocities. See
[Casini and Huerta](https://arxiv.org/html/cond-mat/0610375v2), §3.

Continuity of finite-chain matrices alone does not control the
thermodynamic limit. In particular, taking ell→∞ in the finite-scale
inequality D_h(ell)≤3ell log_2d gives no finite upper bound on c.

## Why the singular qubit neighborhoods matter

Stable gapped phases and regular continuum neighborhoods with controlled
matching are natural places to establish local bounds. These observations
do not exhaust all qubit bonds.

**[Proved example.]** Consider

\[
 H_0=\sum_j n_jn_{j+1},\qquad n_j=(1-Z_j)/2.
\]

On periodic chains with L≥3 it has a gap of one above its zero-energy
space. Its ground words are the binary strings without adjacent ones,
whose number is

\[
 \operatorname{tr}\begin{pmatrix}1&1\\1&0\end{pmatrix}^{L}
 \sim\varphi^L.
\]

Thus a gap above an exponentially degenerate ground space does not make
this a stable gapped phase. Arbitrarily small perturbations can create
nontrivial dynamics inside that space. Bounding nearby c requires control
of that dynamics, rather than a classification of H_0 alone.

The proposed qubit two-valley mechanism and repeated blockade calculations
are relevant to this issue. Their leading band descriptions involve
finitely many modes, but they do not exclude every composite sector,
neutral sector, or hierarchy of further effective constraints. The
[qubit note](../d2_updated_conjecture_20261116/README.md) retains these gaps.

**[Proved restricted obstruction.]** If a nontrivial rank-one onsite
qubit projector M_j is conserved separately at every site, a basis rotation
makes every M_j diagonal. Commuting with all M_j then makes H diagonal in
the product basis. Such an exact stationary one-site marker cannot leave
an internal quantum degree of freedom on its one-dimensional complement.
For larger d the complement can have dimension at least two. This is an
obstruction to a particular qubit cell mechanism, not a theorem excluding
emergent markers or all qubit encodings.

## Higher dimensions and the smaller candidate

The finite-cover theorem and entropy criterion apply to every fixed d.
Their dynamical premise is the unresolved part. One consistent failure
mechanism is a sequence of bonds selecting longer cells with more logical
degrees of freedom, whose velocities and energy scales approach zero.
Each individual Hamiltonian may still have finite nonzero conformal
velocities after its own thermodynamic limit is taken.

The [new d=20 frozen and d=22 mobile candidates](README.md) pursue that
mechanism with explicit finite-cell and leading-band ingredients for
arbitrary n. Their positive zero-port selection margin E_(n+1) shrinks
with n, as does the payload scale. Thus compactness does not contradict
their proposed unbounded c. The mobile defect comparison and conformal
infrared conclusion remain conditional.

Neither a universal finite B(2) nor a universal finite B(d) has been proved
in this discussion. Conversely, the new candidate does not establish an
unconditional violation of such a bound. For a qubit upper-bound attempt,
a useful target is a neighborhood theorem at degenerate bonds with any
finite constant, rather than first trying to prove the sharp value two.

The topological reduction overlaps the earlier
[C16 analysis](../../codex/work_C16/RESULT.md). The entropy criterion and
its connection to the new small-alphabet family are recorded here to keep
the conclusions of this discussion together.
