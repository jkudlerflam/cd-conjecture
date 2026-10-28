# Progress

The task is closed book. No pre-existing files have been read or listed. All work uses the self-contained setting and files created in this directory.

Plan: prove the orbit span, examine the positive cone modulo the free terms, construct finite gate recipes where possible, test the proposed uniformity, and verify small instances numerically. Every unproved assumption will be recorded.

## Analytic findings

1. A Pauli-character average of conjugated Bell projectors isolates any one nonidentity port Pauli coupled to any nonidentity Pauli on s,c. Clifford conjugation suffices. This proves the stated span and gives its dimension (q^2-1)(4q^2-1)+1.
2. Replacing signed character weights by 1 plus or minus the character makes every weight nonnegative and adds only an identity. Therefore positivity is not an obstruction modulo the free identity.
3. The target minus its pure-(s,c) part has q^2-2 Pauli terms, each with coefficient magnitude 2/q^2. The resulting explicit construction has q^2(q^2-2)/2 nonzero frames, each weight 4/q^2, and identity offset 2-4/q^2.
4. The claimed Q-uniformity under conjugation on s is impossible for nonidentity Q, even modulo both free classes. Partial trace over s is invariant under that conjugation but contains the distinct port-c Pauli coupling labeled by Q.
5. Gate order matters at m=1: with Pi_W=W^dagger Pi W, W=Pauli_c SWAP gives no Bell-frame variation. W=SWAP Pauli_c does.

The gate construction uses entangling Clifford gates CNOT(s_l,c). A clarification about this interpretation of the allowed Cliffords has been requested. Numerical verification and a complete written derivation remain underway.

## Completed direct numerical verification

All 21 targets were checked: one at m=1, four at m=2, and sixteen at m=3. The weighted sums were built from the actual gate matrices and rank-two conjugated Bell projectors. The largest entrywise residual was 4.45e-16; the largest Frobenius residual was 2.96e-15. Runtime was about 8 seconds. The m=1 Bell-frame formula also passed, and the incorrectly ordered circuit was confirmed to produce identical frames.

## Repair of uniformity

Let V_Q be controlled-Q from the first qubit to the remaining qubits, and K_Q=V_Q* on s. Although conjugating the synthesized operator by a Pauli on s cannot work, conjugating every gate product W_i by K_Q does work. The proof uses the invariance of Omega under V_Q on p and V_Q* on s. V_Q maps X_1,Y_1 to X_1 Q,Y_1 Q. This construction uses the same positive weights and only Clifford gates. It is a correction to the requested Pauli-frame uniformity, not a proof of that false statement.

The complete gate lists for m=1,2 are being exported. Verification of the corrected uniform construction replays every primitive gate for all Q through m=3.

## Final result

RESULT.md contains the complete proofs, explicit finite indexing of every gate word and weight, time-fraction normalization, the m=1 Bell-frame formula and gate-order correction, the obstruction to Pauli-frame uniformity, and its controlled-Clifford repair. It also proves impossibility under the alternative interpretation allowing only single-qubit Cliffords and SWAP(s_1,c).

The corrected uniform construction passed for all 21 targets through m=3. Its largest Frobenius residual was 9.12e-15 and its largest entrywise residual was 8.89e-16. The full replay took about 8.6 seconds. The direct and repaired constructions were both checked independently against the target definition. Every exported weight is strictly positive. The final consistency check confirmed all 42 verification cases pass a 1e-12 Frobenius tolerance and all frame counts match the proved formula.

No non-Clifford gate or modified Bell-check term is needed under the entangling-Clifford interpretation. No claim is made about optimal frame counts, optimal gate counts, or finite-pulse physical implementation errors. The literal ordinary-Pauli-frame uniformity requirement is impossible, rather than an unfinished derivation.
