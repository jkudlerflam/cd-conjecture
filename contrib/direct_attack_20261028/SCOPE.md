# Frozen sectors, ground-state uniqueness, and what must actually be proved

## Source statement [read]

The original conjecture, [Latorre and Sierra, arXiv:2403.17242v2,
Section 2](https://arxiv.org/html/2403.17242v2#S0.SS2), concerns a single
Hermitian two-site bond, periodic boundaries, one-site invariance of the
Hamiltonian, and a low-energy spectrum described by unitary CFTs. The text
does not impose a unique ground state or explicitly forbid spontaneous
translation breaking. It allows distinct velocities. The paper does not
resolve how its CFT hypothesis should treat immobile defects forced by
incompatible ring lengths.

Consequently, a unique ground state on every ring is a useful validation
target for the mobile construction, not an additional axiom that should be
attributed to the paper. Conversely, a logarithmic entropy slope along one
subsequence does not by itself settle the stated spectral hypothesis.

## A fixed number of sectors does not multiply c [proved]

Suppose a chain has b extremal states, where b is fixed as the system grows.
Write their interval reduced states as rho_a(I), and assume each has

    S(rho_a(I)) = (c/3) log |I| + O(1).

For any fixed probabilities p_a, entropy concavity and the classical
mixture upper bound give

    sum_a p_a S(rho_a(I))
        <= S(sum_a p_a rho_a(I))
        <= sum_a p_a S(rho_a(I)) + H(p).

For completeness, attach a classical register A containing the label a.
The joint entropy is H(p)+sum p_a S(rho_a). Nonnegative classical-quantum
conditional entropy gives the upper bound after discarding A; entropy
concavity gives the lower bound. Since 0<=H(p)<=log b is independent of
|I|, the mixture has the same logarithmic coefficient c/3.

For frozen marker patterns the translated sectors are distinguishable on
the complement of a finite interval in a sufficiently large ring. Tracing
out that complement kills off-diagonal terms between these sectors, so
the reduced state of their coherent translation eigenstate is precisely
such a mixture. This last step requires the local distinguishability of
the frozen patterns; it is not asserted for arbitrary orthogonal states.

There is also a direct blocking observation. If a logical chain is embedded
by isometries acting independently on b-site cells, an interval consisting
of complete cells has exactly the logical entropy. Moving either endpoint
by at most b sites changes entropy by at most 2b log d, by the entropy
triangle inequality. Therefore

    S_physical(ell) = (c/3) log(ell/b) + O(1)
                    = (c/3) log ell + O(1).

This is an exact statement about the embedded reference. Transferring it
to the finite-coupling ground state still needs control of the effective
theory; a norm bound on energy density alone does not provide it.

## Consequence for the present program [research judgment]

The b-fold degeneracy on compatible rings is not by itself an entropy
obstruction. Nor should b identical direct-sum sectors be counted as b
independent tensor-factor CFTs with central charge bc.

The stronger obstruction recorded in the frozen-cell audit is different:
on incompatible lengths, immobile defects give a number of ground states
growing with L and additional low-energy structure. The finite-b entropy
argument above does not handle that. It also does not establish the
payload's finite-coupling n-boson limit. The mobile d=24 and d0=4464 routes
continue to be the cleanest way to avoid this interpretation issue.

Accordingly, neither an assertion that symmetry breaking alone refutes the
candidate nor an assertion that compatible-ring entropy alone solves the
conjecture is justified by this note.

Provenance: coordinating Codex GPT-6-family agent, 2026-10-28. The source was
opened directly; the entropy argument is derived here. No human review.
