# Exact-encoding shortcuts and the depth required to use them

## Claim and status

[proved, AI derivation; not independently reviewed] An exact finite-depth
circuit dressing does not remove the separated-boundary obstruction to
encoding a free multichannel hopping chain, unless the circuit reaches a
fixed fraction of the cell length. More generally, an exactly invariant
product code for an exponentially local transformed Hamiltonian requires
interaction tails at the middle of a cell of at least the logical hopping
scale. These are obstructions to proposed shortcuts, not bounds on the
central charge of general nearest-neighbour chains.

The precise theorems and proofs are in [PROOF.md](PROOF.md). The finite-depth
result extends the product-code theorem already proved in
[`codex/work_C8/RUN8_RESULT.md`, section 10](../../../codex/work_C8/RUN8_RESULT.md).
The old theorem is credited rather than presented as a new result.

No new counterexample, microscopic conformal limit, or general proof of
c <= d-1 has been obtained. The strongest current constructive route still
requires a controlled low-energy reduction. An exactly invariant encoding
cannot simply be substituted for that step without checking the locality
conditions below.

## Main statements

For a microscopic nearest-neighbour chain, a circuit of depth D with gates
on at most r consecutive sites has propagation radius R = D(r-1). Suppose
cells contain b sites, there are at least four cells, and the entire
circuit-dressed tensor-product code is exactly invariant. If

    b >= 4 D(r-1) + 2,

then all adjacent *connected* logical bonds commute. They cannot equal the
standard free hopping bond

    K = J sum_(a=1)^(2n) (i Gamma_a P) tensor Gamma_a,  J != 0.

Thus a nontrivial exact realization of this target needs

    D(r-1) > (b-2)/4.

For an exponentially local transformed interaction Phi with

    M_mu = sup_x sum_(X contains x) exp(mu diam(X)) ||Phi(X)||,

exact invariance and that same exact hopping target imply, with
R = floor(b/2),

    |J| <= (sqrt(3)+1) b M_mu exp(-mu R).

This last inequality is compatible with the extremely small hopping scales
in the repository's proposed constructions. It is not a refutation of them.
It makes explicit why a sufficiently strong exact decoupling claim must be
compared with the logical bandwidth, rather than only with a microscopic
energy unit.

## Verification and provenance

Author: a Codex GPT-6-family subagent, 2026-10-28, in the campaign coordinated
by the chat "Identify a Conjecture Proof Strategy". The proofs are independent
extensions of the cited C8 argument. No human review is claimed.

`python3 contrib/decisive_attack_20261028/exact/checks.py` checks the Clifford
identities and exact bond-commutator norm at n=1,2,3, and the collar arithmetic
for a finite set of circuit parameters. It uses numpy and completes in a
few seconds. These checks do not replace the proofs or test all circuits.
