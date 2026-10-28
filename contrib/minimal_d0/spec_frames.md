# Positive-weight synthesis of a port coupling from one Bell-check term
Closed book: do not read or list any files you did not create. Work only inside the current working directory.
Python (numpy/scipy) allowed; keep total runtime under about 60 minutes. Be honest: a precisely located obstruction
is as valuable as a construction. Flag every step you assume rather than prove. Write PROGRESS.md as you go and
RESULT.md at the end.

## Setting (self-contained)
Registers: a port p and a partner s, each of dimension q = 2^m, viewed as m qubits p_1..p_m and s_1..s_m; and one
qubit c. Let |Omega> = tensor_l |Phi+>_{p_l s_l} with |Phi+> = (|00>+|11>)/sqrt2, Pi = |Omega><Omega| (x) 1_c.
Gates may act on s (x) c only, never on p. For a unitary W on s (x) c write Pi_W = W^dagger Pi W.

Grading on p: Gamma = Z_{p_1}. Odd port operators are the Pauli strings F = A (x) Q with A in {X, Y} on p_1 and Q a
Pauli string on p_2..p_m (r = q^2/2 of them). For a fixed Q define the target coupling
    T_Q = sum_{A in {X,Y}} [ (A (x) Q)_p Pi_Omega (x) A_c  +  h.c. ]        (operator on p (x) s (x) c).
(At m = 1, T = (X_p X_c + Y_p Y_c)/2 + (X_s X_c - Y_s Y_c)/2 up to normalization conventions; check this.)

Two classes of terms are "free":
 (F1) any Hermitian operator acting on s (x) c only (it can be added or cancelled separately, with either sign);
 (F2) any operator of the form D (x) 1_c with D diagonal in the Bell basis of p (x) s, i.e. a real combination of the
      projectors onto (1 (x) P)|Omega>, P a Pauli string on s.

## Questions
1. Span. Prove (or find the exact statement that holds) that the real span of {Pi_W : W unitary on s (x) c} equals the
   set of Hermitian H on p (x) s (x) c whose partial traces tr_p H and tr_{s,c} H are both multiples of the identity.
   (Numerically the ranks are 46 = 64 - 15 - 3 at q = 2 and 946 = 1024 - 63 - 15 at q = 4.) Deduce that T_Q minus its
   pure-(s,c) part lies in the span.
2. Positive cone. The physical construction can only produce sum_i w_i Pi_{W_i} with ALL w_i >= 0 (time fractions),
   plus free terms (F1), (F2) of any sign. For m = 1 show that
       sum_i w_i Pi_{W_i} = alpha*1 + g * T_port + (F1) + (F2),   g != 0,
   is achievable with W_i = SWAP(s,c) followed by Paulis on c (I believe weights (f1, f2, f, f) on the four Bell
   frames of (p,c), with 2f = f1 + f2, give exactly -(f1-f2)/4 (X_p X_c + Y_p Y_c); verify). Then for general m and each
   Q, either give an explicit finite list of unitaries W_i on s (x) c, each a short product of gates from a SMALL FIXED
   set (Cliffords on s_l and c, SWAP(s_1,c), and at most one fixed non-Clifford rotation if needed), with explicit
   positive weights, achieving g T_Q + (F1) + (F2) + alpha*1 with g != 0 and NO other junk; or prove that no such
   positive combination exists and identify the smallest extra resource that fixes it (for example letting the
   Bell-check term carry a factor Y_c, i.e. using Pi_Omega (x) Y_c instead of Pi_Omega (x) 1_c).
3. Uniformity. The different Q (there are 4^(m-1) of them) should be reachable by conjugating ONE fixed recipe with a
   Pauli frame on s_2..s_m (which is itself of type (F2) when applied alone). Check that conjugating the m-qubit recipe
   for Q = identity by the Pauli frame Q^T on s gives the recipe for T_Q, with all intermediate frames still of the
   allowed types, and count the gates per Q.
4. Verify numerically for m = 1 and m = 2 (and m = 3 if cheap): build the weighted sum, subtract the free parts, and
   report the residual against g T_Q.

Deliverable: RESULT.md with the statements, proofs or obstructions, the explicit recipes (gate lists and weights),
numerical residuals, and a short list of what remains assumed.
