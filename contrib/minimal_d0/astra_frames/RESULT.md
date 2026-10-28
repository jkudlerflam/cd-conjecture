# Positive-weight synthesis from one Bell check

All statements below use the normalized maximally entangled state in the question. The work is closed book. No pre-existing files were read or listed, and all generated files are in the working directory.

The span statement is true. Positive synthesis is possible for every Q using Clifford gates, with no non-Clifford rotation and no modified Bell check. The proposed uniformity under an ordinary Pauli frame on s is false. A controlled-Q Clifford conjugation of the gate words repairs uniformity, with the same weights for every Q.

There is one interpretation of the gate allowance to make explicit. The construction uses entangling Clifford gates CNOT(s_l,c), in either direction. If “Cliffords on s_l and c” means only single-qubit Cliffords, together with SWAP(s_1,c), synthesis is impossible for every Q when m >= 2. Section 8 proves this alternative obstruction. The mathematical span statement concerns all unitaries and is unaffected by this distinction.

## 1. Conventions and the result

Write B=s tensor c, q=2^m, and let P_m denote the q^2 Hermitian Pauli strings without phases. Define

\[
\eta(P)=(-1)^{\#Y(P)},\qquad P^T=\eta(P)P,
\]

and

\[
\chi_R(U)=\begin{cases}+1,&UR=RU,\\-1,&UR=-RU.\end{cases}
\]

All identity operators below have the dimensions implied by their subscripts. An unlabelled identity acts on the whole system. Define the pure-(s,c) part of the target by

\[
 L_Q=I_p\otimes\frac{2}{q^2}
 \sum_{A=X,Y}(A\otimes Q)_s^T\otimes A_c.
 \tag{1}
\]

The explicit construction in Section 4 gives

\[
 \sum_{i=1}^{M} \frac{4}{q^2}\Pi_{W_i}
 =\left(2-\frac{4}{q^2}\right)I+T_Q-L_Q,
 \qquad M=\frac{q^2(q^2-2)}{2}.
 \tag{2}
\]

Thus g=1, the F1 term is -L_Q, and the F2 term can be zero. Every listed weight is strictly positive. Repetitions, if any, can be merged by adding weights. This is a construction, not a minimum-frame claim.

If the weights must sum to one, divide (2) by 2(q^2-2). Each frame then has weight 1/M and

\[
 \frac1M\sum_i\Pi_{W_i}
 =\frac{I}{q^2}+\frac{T_Q-L_Q}{2(q^2-2)}.
 \tag{3}
\]

In particular, the normalized coupling remains nonzero for every m >= 1. Both signs of the coupling are possible by reversing the character selection below.

## 2. Exact orbit span

The normalized Bell projector has the Pauli expansion

\[
 \Pi=\frac{1}{q^2}\sum_{P\in P_m}
 P_p\otimes P_s^T\otimes I_c.
 \tag{4}
\]

For completeness, this follows by expanding in the orthogonal Pauli basis. The identity

\[
 \langle\Omega|A\otimes B|\Omega\rangle
 =\frac1q\operatorname{tr}(A^TB)
\]

makes the coefficient of P tensor R vanish unless R=P, with coefficient eta(P)/q^2 when R=P.

For every unitary W on B,

\[
 \operatorname{tr}_p\Pi_W=\frac{I_B}{q},\qquad
 \operatorname{tr}_B\Pi_W=\frac{2I_p}{q}.
 \tag{5}
\]

The first equality follows from tr_p Pi=I_B/q and conjugation by W. The second follows from invariance of the trace over B. This proves the necessary condition for the real span.

For sufficiency, a Hermitian operator with both marginals scalar has precisely the form

\[
 H=aI+\sum_{R\in P_m\setminus\{I\}}
       \sum_{D\in P_{m+1}\setminus\{I\}}
       h_{RD}\,R_p\otimes D_B,
 \qquad a,h_{RD}\in\mathbb R.
 \tag{6}
\]

Indeed, its coefficients of I_p tensor D_B and R_p tensor I_B vanish by the two partial-trace conditions. The scalar marginal coefficients are already related by the total trace, so they are not independent constraints on a.

Here is a finite construction of every basis element in (6). Fix R != I and any Clifford C on B. Equation (4) gives

\[
 \Pi_{U_s C}
 =\frac1{q^2}\sum_{P\in P_m}
 \chi_P(U)\,P_p\otimes C^\dagger(P_s^T\otimes I_c)C.
 \tag{7}
\]

The Pauli characters satisfy

\[
 \frac1{q^2}\sum_{U\in P_m}\chi_R(U)\chi_P(U)
 =\delta_{RP}.
 \tag{8}
\]

To prove orthogonality when R != P, their character product is a nontrivial Pauli character. Choose a Pauli anticommuting with its label and pair terms under multiplication by that Pauli, ignoring phases. The paired signs cancel. For R=P every summand is one.

Consequently,

\[
 \frac1{q^2}\sum_U\chi_R(U)\Pi_{U_s C}
 =\frac1{q^2}R_p\otimes C^\dagger(R_s^T\otimes I_c)C,
 \tag{9}
\]

while the same average without the character is I/q^2. Section 5 constructs a Clifford C with

\[
 C^\dagger(R_s\otimes I_c)C=D_B
 \tag{10}
\]

for every nonidentity D. Multiplying (9) by q^2 eta(R) gives R_p tensor D_B. This proves sufficiency using only a finite subset of the Clifford orbit. No irreducibility assumption or numerical rank argument is needed.

The dimension is

\[
 1+(q^2-1)(4q^2-1)=4q^4-5q^2+2.
 \tag{11}
\]

It is 46 at q=2 and 946 at q=4, agreeing with the ranks in the question.

There is also a useful general positivity consequence. For t=+1 or -1, replacing the signed weights in (9) by 1+t chi_R(U) gives

\[
 \frac1{q^2}\sum_U[1+t\chi_R(U)]\Pi_{U_s C}
 =\frac{I}{q^2}+\frac{t\eta(R)}{q^2}R_p\otimes D_B.
 \tag{12}
\]

Every weight is nonnegative. Applying this to each term of (6) proves that every element of the orbit span becomes a positive combination after adding a sufficiently large scalar identity. Thus positivity cannot obstruct the present target once the identity is free.

## 3. Target decomposition and the one-qubit check

Put F_A=A tensor Q. Equations (4) and the definition of T_Q give

\[
 T_Q=\frac1{q^2}\sum_{A=X,Y}\sum_{P\in P_m}
 \{F_A,P\}_p\otimes P_s^T\otimes A_c
 =\frac2{q^2}\sum_{A=X,Y}\sum_{[P,F_A]=0}
 (F_AP)_p\otimes P_s^T\otimes A_c.
 \tag{13}
\]

The product F_AP is Hermitian in the second sum because the two factors commute. The only terms with identity on p are P=F_A. They give (1). Equivalently,

\[
 L_Q=I_p\otimes\frac1q\operatorname{tr}_p T_Q.
\]

Every remaining term has nonidentity Pauli factors on both p and B. Therefore T_Q-L_Q has zero partial trace over either factor and lies in the span proved above.

A nonidentity Pauli commutes with q^2/2 strings. Removing P=F_A from each of the two A sums leaves

\[
 N=q^2-2
\]

distinct Pauli products, each with coefficient magnitude 2/q^2. They are distinct within each A sum because P determines both factors, and between the two sums because their c factors differ.

At m=1, the commuting strings are {I,X} for A=X and {I,Y} for A=Y. Hence, with the normalization in the question,

\[
 T=\frac12(X_pX_c+Y_pY_c)
   +\frac12(X_sX_c-Y_sY_c).
 \tag{14}
\]

This is an equality, with identities on omitted registers. Define T_port=(X_pX_c+Y_pY_c)/2.

Let S_sc be SWAP(s,c), and use the four words

\[
 W_I=S_{sc}I_c,\quad W_Z=S_{sc}Z_c,\quad
 W_X=S_{sc}X_c,\quad W_Y=S_{sc}Y_c.
 \tag{15}
\]

Their projectors are the Bell projectors on (p,c), tensored with I_s, in the order Phi+, Phi-, Psi+, Psi-. Their Pauli expansions are

| Frame | 4 times its projector on p,c |
| :--- | :--- |
| I | I + XX - YY + ZZ |
| Z | I - XX + YY + ZZ |
| X | I + XX + YY - ZZ |
| Y | I - XX - YY - ZZ |

Assign weights (f,f,f2,f1) in this order, with f1+f2=2f. Then

\[
 f\Pi_{W_I}+f\Pi_{W_Z}+f_2\Pi_{W_X}+f_1\Pi_{W_Y}
 =fI-\frac{f_1-f_2}{4}(X_pX_c+Y_pY_c).
 \tag{16}
\]

Thus the suggested coefficient is correct if the order carrying weights (f1,f2,f,f) is (Psi-,Psi+,Phi+,Phi-). Assigning those weights to (Phi+,Phi-,Psi+,Psi-) would instead give XX-YY. The choice f1=1, f2=3, f=2 has all weights strictly positive and gives 2I+T_port=2I+T-L.

There is a gate-order qualification. Throughout, the rightmost gate acts first. With Pi_W=W^dagger Pi W, the words in (15) apply the Pauli to c before the SWAP. If W=P_c S_sc means “SWAP followed by a Pauli on c” in chronological circuit order, then

\[
 \Pi_{P_cS_{sc}}=S_{sc}^\dagger\Pi S_{sc}
\]

for every P_c, since P_c commutes with the original Pi. Those four words give the same projector. Applying a Pauli on s after SWAP, or writing the conjugated Bell frames directly as in (15), gives the intended construction.

## 4. Explicit positive recipe for general m and Q

The following indexing rules specify the complete finite list. There is no search or continuous parameter selection.

For every A in {X,Y} and every P in P_m commuting with F_A, except P=F_A, write

\[
 F_AP=\epsilon R,\qquad \epsilon\in\{+1,-1\},\quad
 R\in P_m\setminus\{I\}.
\]

Set

\[
 D=P_s\otimes A_c,\qquad \sigma=\epsilon\eta(P),\qquad
 t=\sigma\eta(R).
 \tag{17}
\]

Thus this target term is (2/q^2) sigma R_p tensor D_B. Obtain C from the deterministic gate list in Section 5, so C^dagger(R_s tensor I_c)C=D.

For this pair (A,P), include all the words

\[
 W_{A,P,U}=U_s C,
 \qquad U\in P_m\text{ such that }\chi_R(U)=t,
 \qquad w_{A,P,U}=\frac4{q^2}.
 \tag{18}
\]

There are q^2/2 such U because R is nonidentity. To verify the resulting sum, write it as an unrestricted character average:

\[
\begin{aligned}
 \sum_{\chi_R(U)=t}\frac4{q^2}\Pi_{U_sC}
 &=\frac2{q^2}\sum_U[1+t\chi_R(U)]\Pi_{U_sC}\\
 &=\frac2{q^2}I+\frac{2t\eta(R)}{q^2}R_p\otimes D_B\\
 &=\frac2{q^2}I+\frac{2\sigma}{q^2}R_p\otimes D_B.
\end{aligned}
\tag{19}
\]

Summing (19) over the N target terms proves (2), including the absence of any other terms. The identity trace also checks the normalization. The total weight is 2N, each projector has rank two, and both sides of (2) have trace 4N.

The list length is M=q^2(q^2-2)/2, namely 4, 112, and 1984 for m=1,2,3. The number of frames grows as q^4. Only the length of each gate word is linear in m.

## 5. Deterministic short Clifford gate lists

The elementary gate families are H, S, S^dagger, single-qubit Paulis, CNOT on the pairs (s_l,c) in either direction, and the allowed SWAP(s_1,c). Here S=diag(1,i). All are Clifford. These are fixed gate types independent of m, with the addressed qubits specified by the recipe.

For any unsigned nonidentity Pauli D on s,c, construct L(D) so that

\[
 L(D)D L(D)^\dagger=Z_c.
 \tag{20}
\]

The following lists are chronological.

1. On each X factor, apply H. On each Y factor, apply S^dagger and then H. On each Z or I factor, do nothing. The Pauli is now a product of Zs, with no minus sign. In particular, S^dagger Y S=X.
2. If c is absent from the support, choose the lowest-index supported s_j and swap it with c. Implement this swap as CNOT(s_j,c), CNOT(c,s_j), CNOT(s_j,c). This avoids requiring a native SWAP other than the one allowed in the question.
3. For every remaining supported s_j, apply CNOT(s_j,c), in increasing index order.

Step 3 removes each Z_sj. Under that CNOT, Z_c maps to Z_sj Z_c and Z_sj maps to itself, so their product maps to Z_c. This proves (20), including its sign.

Now define

\[
 C=L(R_s\otimes I_c)^\dagger L(D).
 \tag{21}
\]

Using (20),

\[
 C^\dagger(R_s\otimes I_c)C
 =L(D)^\dagger Z_cL(D)=D.
\]

Thus the chronological word for each frame (18) is

\[
 [\text{gates of }L(D)],\quad
 [\text{reversed inverse gates of }L(R_s\otimes I_c)],\quad
 [\text{nonidentity single-qubit factors of }U_s].
 \tag{22}
\]

This is the promised explicit gate list for arbitrary m and Q.

For the source R_s tensor I_c, L uses at most 2m local gates, three CNOTs for the swap, and m-1 eliminating CNOTs, hence at most 3m+2 gates. The target D=P_s tensor A_c has c in its support, so its L uses at most 2(m+1) local gates and m CNOTs, again at most 3m+2. Appending U costs at most m gates. Consequently,

\[
 \text{gates per frame}\le 7m+4,
 \qquad \text{gates for the full direct list}\le M(7m+4).
 \tag{23}
\]

The count treats Paulis and S^dagger as elementary Clifford gates. They can instead be compiled into H and S with a constant-factor overhead. No cancellation between neighboring frames is assumed in these bounds.

## 6. Why Pauli-frame uniformity is false

A partial-trace invariant gives an obstruction independent of the construction. Directly from the target definition,

\[
 \operatorname{tr}_sT_Q
 =\frac2q\sum_{A=X,Y}(A\otimes Q)_p\otimes A_c.
 \tag{24}
\]

Conjugating any operator by any unitary acting only on s leaves its partial trace over s unchanged. The free terms cannot change the nonidentity-p, nonidentity-c component of this trace. For F1, the trace is I_p tensor an operator on c. For F2, each normalized Bell projector has s-trace I_p/q, so the trace is a scalar multiple of I_pc.

For Q != I, the Pauli strings appearing in (24) are orthogonal to those for Q=I. It follows that no s-only conjugation of a realization of g T_I, with g != 0, can yield a realization of any nonzero multiple of T_Q modulo F1, F2, and the scalar identity. In particular,

\[
 (Q_s^T)^\dagger T_I Q_s^T\ne T_Q
 \quad\text{modulo F1 and F2}.
 \tag{25}
\]

This remains impossible for an arbitrary s-only unitary, not just a Pauli. A modified Bell check would not change this invariant of the proposed conjugation procedure.

The fact that Pi_(Q_s^T) is itself a Bell-diagonal F2 term is true. It does not imply (25) is an equality. In the Bell-state picture, T_I has links from the identity-labelled Bell state to the A-labelled states. A spectator Pauli moves both endpoints to Q and AQ, whereas T_Q requires endpoints I and AQ.

There is also no alternative rescue by conjugating the gate words themselves with this ordinary spectator Pauli. Under W_i -> Q_s^T W_i Q_s^T, the sum is conjugated jointly by Q_p and Q_s^T. This joint conjugation fixes both Omega and F_A=A_1 for the Q=I target, so T_I is unchanged. F1 remains F1 and F2 is merely permuted.

Thus the literal Pauli-frame uniformity requirement cannot be met. This is a proved obstruction, not an unverified part of the construction.

## 7. Corrected uniformity using one controlled-Pauli Clifford

A small change to the uniformity operation works without adding a new gate type. On m qubits define

\[
 V_Q=|0\rangle\langle0|_1\otimes I
       +|1\rangle\langle1|_1\otimes Q,
 \qquad K_Q=V_Q^*\text{ on }s.
 \tag{26}
\]

Since Q is a Hermitian Pauli, V_Q is a Hermitian unitary and

\[
 V_QX_1V_Q^\dagger=X_1Q,\qquad
 V_QY_1V_Q^\dagger=Y_1Q.
 \tag{27}
\]

For example, the two off-diagonal control blocks of X_1 become Q, and the corresponding blocks of Y_1 become -iQ and iQ. This proves both identities directly. V_Q is Clifford because it is a product of controlled single-qubit Paulis with the same control.

Let U_Q=(V_Q)_p tensor (K_Q)_s tensor I_c. The maximally entangled state is invariant under V tensor V* for every unitary V, so

\[
 U_Q\Pi U_Q^\dagger=\Pi,\qquad
 U_QT_IU_Q^\dagger=T_Q,\qquad
 U_QL_IU_Q^\dagger=L_Q.
 \tag{28}
\]

For any gate word W on s,c, define

\[
 W^{(Q)}=K_QWK_Q^\dagger.
 \tag{29}
\]

Then

\[
 \Pi_{W^{(Q)}}=U_Q\Pi_WU_Q^\dagger.
 \tag{30}
\]

To see the order in this identity, use the invariance of Pi to write V_p Pi V_p^dagger=K_Q^dagger Pi K_Q. Since W commutes with V_p,

\[
 U_Q W^\dagger\Pi WU_Q^\dagger
 =K_QW^\dagger K_Q^\dagger\Pi K_QWK_Q^\dagger
 =\Pi_{K_QWK_Q^\dagger}.
\]

Applying (30) to the fixed Q=I recipe proves (2) for every Q with unchanged weights. The V_Q on p is only a device in the proof. No implemented gate acts on p. Also, (29) is a conjugation of the gate words, not a conjugation of the synthesized operator by K_Q alone, so it does not contradict Section 6.

F1 is preserved by the joint conjugation U_Q. F2 is preserved as well. In fact,

\[
 U_Q(I\otimes P)|\Omega\rangle
 =(I\otimes K_Q P K_Q^\dagger)|\Omega\rangle,
\]

and Clifford conjugation sends a Pauli to a signed Pauli, hence permutes the Bell projectors. Pauli factors U_s used in the seed construction also remain Paulis under conjugation by K_Q. However, K_Q itself is generally not an F2 frame. This is a controlled-Clifford replacement of the requested ordinary Pauli frame.

An explicit chronological circuit for K_Q is as follows. If Q=I use the empty circuit. Otherwise first apply SWAP(s_1,c), apply the following gates for every nonidentity factor Q_j in increasing j, and apply SWAP(s_1,c) again.

| Q_j | Chronological gates between the two swaps |
| :--- | :--- |
| X | CNOT(c,s_j) |
| Y | S(s_j), CNOT(c,s_j), S^dagger(s_j) |
| Z | H(s_j), CNOT(c,s_j), H(s_j) |

The Y line implements controlled(-Y), as required by Q*, because S^dagger X S=-Y. The other two lines implement controlled X and controlled Z. The swaps move the control from s_1 to c and restore it afterward, so the completed K_Q acts trivially on c.

For n_X,n_Y,n_Z nonidentity factor counts, the K_Q word has

\[
 \ell(K_Q)=2+n_X+3n_Y+3n_Z\le2+3(m-1)
 \tag{31}
\]

gates when Q != I. The chronological list for W_i^(Q) is the inverse K_Q list, the seed W_i list, and the K_Q list. Thus the extra cost per frame is 2 ell(K_Q), with the bound

\[
 \text{gates per corrected-uniform frame}\le13m+2.
 \tag{32}
\]

The total list still has M frames. The full per-Q cost is at most M(13m+2), without optimizing transitions between frames. There is no non-Clifford rotation and no need to replace Pi by Pi_Omega tensor Y_c.

## 8. Alternative interpretation with only single-qubit Cliffords

Suppose the allowed gates are restricted to single-qubit Cliffords on s,c and SWAP(s_1,c), with no entangling Clifford. Every word is local gates followed or preceded by either the identity permutation or this swap. The spectator qubits never interact with s_1,c.

Without the swap, Pi_W has identity on c. With the swap, Pi_W has identity on s_1, because the identity originally on c has moved to s_1. Local conjugations preserve these identities. Therefore every obtainable linear combination, even allowing signed coefficients, has zero coefficient on every Pauli string that is nonidentity on both s_1 and c.

For m >= 2, choose a spectator Pauli P_0 commuting with Q and different from Q. Such a choice always exists. Use P_0=I if Q != I, and any nonidentity P_0 if Q=I. In (13), take A=X and P=X_1 tensor P_0. The resulting target term is

\[
 \frac2{q^2}
 [I_1\otimes QP_0]_p\otimes
 [X_1\otimes P_0^T]_s\otimes X_c.
 \tag{33}
\]

It is nonzero on p, s_1, and c. Its coefficient cannot be cancelled by another target term, since the Pauli expansion is unique. F1 cannot supply it because its p factor is nonidentity. F2 cannot supply it because its c factor is nonidentity. This proves impossibility under the strictly local gate interpretation, for every Q at m >= 2.

Adding the single entangling gate type CNOT, addressable on the star pairs (s_l,c), suffices by Sections 4 and 5. At least an entangling gate type is necessary under this interpretation. No claim is made that the star uses the minimum number of physical couplers or that these circuits minimize gate count.

## 9. Numerical verification

The scripts form the target directly from the normalized rank-one Bell projector. Independently, they construct each W from its gate list and build Pi_W as a rank-two projector. For example its two columns are (I_p tensor W^dagger)(Omega tensor |a>_c), a=0,1. The residual tested is

\[
 E_Q=\sum_iw_i\Pi_{W_i}
 -\left(2-\frac4{q^2}\right)I+L_Q-T_Q.
 \tag{34}
\]

F2 is zero. The following are maxima over every Q at the indicated m, with double-precision arithmetic.

| Direct recipe | Q values tested | Frames per Q | Largest gate count observed | max entrywise residual | max Frobenius residual | max relative Frobenius residual |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| m=1 | 1 | 4 | 8 | 4.441e-16 | 1.410e-15 | 4.986e-16 |
| m=2 | 4 | 112 | 15 | 4.441e-16 | 2.960e-15 | 1.047e-15 |
| m=3 | 16 | 1984 | 22 | 2.221e-16 | 2.604e-15 | 9.205e-16 |

The corrected uniform recipe was checked separately by replaying every primitive gate, including the K_Q sandwiches, for all 21 targets.

| Corrected uniform recipe | Q values tested | Largest gate count observed | max entrywise residual | max Frobenius residual | max relative Frobenius residual |
| :--- | ---: | ---: | ---: | ---: | ---: |
| m=1 | 1 | 8 | 4.441e-16 | 1.410e-15 | 4.986e-16 |
| m=2 | 4 | 25 | 8.882e-16 | 5.543e-15 | 1.960e-15 |
| m=3 | 16 | 38 | 6.662e-16 | 9.116e-15 | 3.223e-15 |

Additional checks and their outcomes follow.

- The Pauli target expansion agreed with the direct target minus L_Q to at most 2.221e-16 entrywise. All targets had Frobenius norm sqrt(8), as follows from the orthogonality of their two Bell-state links and c Paulis.
- Every normalizer sign and every mapping C^dagger(R^T tensor I)C=eta(R)D was checked. The largest errors were 6.662e-16 and 8.882e-16, respectively.
- The implemented K_Q agreed with controlled-Q* tensor I_c to at most 3.331e-16.
- The trace sum rule had absolute error at most 2.843e-14, at total trace 248 for m=3.
- The four-frame m=1 example with weights (2,2,3,1) gave 2I+T_port to 4.441e-16. The incorrectly ordered words P_c SWAP produced identical projectors, with zero observed variation.
- The partial-trace obstruction was checked for every nonidentity Q. The Frobenius norm of tr_s(Q_s^T T_I Q_s^T-T_Q) was sqrt(32/q), namely sqrt(8) at m=2 and 2 at m=3. These nonzero values agree with (24), not with numerical roundoff.

The direct verification took about 7.9 seconds and the full corrected-uniform verification about 8.6 seconds. Both are local, small dense-matrix calculations. The initial invocation using `python` failed because that executable was unavailable. The reported runs used `python3` and completed successfully.

The reproducible artifacts are

- `verify.py`, the direct construction and verification;
- `recipes.py`, a generator of chronological gate lists for either construction, and the independent full-gate replay of the corrected uniform recipe;
- `explicit_recipes_m1_m2.json`, every gate and positive weight for the corrected uniform recipe at m=1 and m=2, for every Q;
- `verification.json` and `uniformity_verification.json`, per-Q numerical results;
- `PROGRESS.md`, the work record.

In the JSON gate lists, sites 0 through m-1 mean s_1 through s_m, and site m means c. `Sd` means S^dagger. Gates are in chronological order. For arbitrary m, `make_recipe(m, qlabels, uniform=True)` in `recipes.py` returns the complete list; labels 0,1,2,3 mean I,X,Y,Z. Setting `uniform=False` uses the direct Q-dependent recipe (18).

Reproduce the numerical checks with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 recipes.py
```

## 10. What is proved, assumed, and not claimed

The orbit span, the positive-weight identity (2), the gate decomposition, the failure of ordinary Pauli-frame uniformity, the controlled-Clifford repair, and the strictly local gate obstruction are proved above. Numerical agreement is a separate consistency check.

The construction assumes the entangling-Clifford interpretation of the allowed gate set. The alternative interpretation is handled by the obstruction in Section 8. As in the question, the mathematical synthesis assumes arbitrary signed F1 and F2 corrections and an exact weighted sum of conjugated Bell-check terms. Physical implementations with finite-duration control pulses or a finite-step Hamiltonian averaging protocol would require a separate error analysis. None is claimed here.

There is no unproved mathematical step in the stated construction. Optimal frame count, optimal gate count, and the largest possible normalized coupling are not determined. The literal Pauli-frame requirement remains impossible; the controlled-Q sandwich is the explicit replacement.
