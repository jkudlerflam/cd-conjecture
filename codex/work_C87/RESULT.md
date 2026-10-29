# C87 result

[Result] No <=10-letter bond with a certified member n>=a was obtained.
The fallback is four obstructions beyond C86's single-token count, plus
an exact example showing why general interacting ports remain open.
There is no universal eleven-letter lower bound here.

[Proved] Arbitrarily many tokens on an open path cannot generate a
configuration-space flux through homogeneous nearest-neighbor exchanges
alone. For a swap st -> ts with phase chi_st, a diagonal gauge built
from the number of inverted st pairs removes every phase. This holds
with interacting diagonal energies and with cycles in configuration
space. For two tokens, forbidding exchange conserves their order;
allowing homogeneous exchange connects orders but still gives no flux.
If the common-carrier reduction retains point sources and diagonal
conversion projectors in this gauge, its spectrum is reciprocal.
C84's pair-sum margin then vanishes for n>=2. Species-changing conversion
or contextual exchange can escape these hypotheses.

[Proved] A fixed number k of tokens whose positions store all the data
has code dimension at most O(b^k), so it cannot retain n=b-O(1) independent
bits. Ten bits require at least 33 sites for two distinguishable tokens
on one chain, or 46 for identical tokens. A fixed multiset of charged
letters also makes the additive charge scalar on the code, preventing
C84's logical U(1) identification. Connected hopping components have
unique ground states after the phase gauge. These are separate
obstructions, not an exclusion of every longer or fragmented encoding.

[Proved] A diagonal phase dressing U of a connected quadratic hopping
operator preserves Gaussianity on the full Fock space if and only if

    U=exp(i F(N)) exp(i sum_j theta_j n_j).

The F(N) factor commutes with the quadratic operator; the other factor
only rephases modes. Thus this class of occupation-dependent port phases
cannot produce new band data while remaining Gaussian. An explicit
conditional phase fails by operator-norm distance at least one from
every Gaussian, and a product amplitude dressing fails by at least 9/4,
in the dimensionless second-order coefficients tested.

[Proved] A frozen marker storing k charged binary bits, with injective
payload-preserving source maps of charge shifts +/-delta/2, needs at least

    2^(k+1)+binomial(k,floor(k/2))

marker states, even allowing the two source images to overlap. One bit
costs five marker letters. Keeping four L/R backgrounds and two carrier
states therefore still costs at least eleven letters. Contextual markers
and ports that move the marker payload elsewhere are outside this count.

[Proved positive algebraic example] Four orthogonal equal-gap channels
with sources P L_1, (1-P)L_1, Q L_2, (1-Q)L_2 give

    K_2=-L_1^dagger L_1-L_2^dagger L_2,

although each individual channel contributes interactions. Here
L_1=f_0+f_1, L_2=f_2+f_3, P=n_2, Q=n_0. This identity holds on the full
16-dimensional code. It rules out extending the dressing obstruction
to all data-dependent ports. It has no microscopic alphabet, cell
selector, or n>=a band certificate. Its fourth-order coefficient is
interacting.

[Proved, computer-assisted] Integer gauge checks include a graph with
210 configurations and 241 independent cycles. A 36-dimensional exact
product-resolvent inverse checks the multi-token scalar coefficient.
All 256 diagonal sign dressings on three modes satisfy the proved
classification, with 64 Gaussian-preserving cases. The port identities
and interaction witnesses use exact integer matrices.

[Numerical] Two-token carrier ED uses dimensions 6 through 72. Port
gadget ED uses dimensions 32,64,80, with all 16 code levels checked at
lambda=0.1,0.05,0.025. Second-order energy errors scale as lambda^4;
port eigenpair residuals are below 3.3e-15. These are finite gadgets,
not ED of a new nearest-neighbor letter bond.

[Rejected checker, repaired] The first rational inverse failed because
NumPy fixed-width integers entered Fraction arithmetic. Conversion to
Python integers repairs it, and the exact inverse identity now passes.
The failure is retained in rejected_checks.json.

[Completed] Derivations and hypotheses are in [PROOF.md](PROOF.md).
Run `python3 checks.py` for all new exact checks and small ED, and
`python3 input_audit.py` for the input inventory and inherited rational
sign/budget audits. Outputs are checks.json and input_inventory.json.
No new root certificate, finite-K schedule, selector, or central-charge
claim is made. No independent referee or subagent was used. All new
files are inside work_C87.

[Open] A construction below eleven must evade these hypotheses, for
example by coordinated interacting channels, conversion-dependent
configuration phases, or a new marker and charge encoding. The channel
cancellation above is a surviving algebraic ingredient, not a completed
design.
