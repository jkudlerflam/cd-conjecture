# Task: adversarial referee of a claimed theorem (uniform decay of the ring energy defect), then an attempt to improve its rate

You are an independent referee. Below are (I) the definitions, (II) a claimed proof produced by another model, reproduced verbatim, and (III) your tasks. Assume the proof may contain an error. Your first duty is to find it. Do not trust any step because it is tagged [rigorous].

## (I) Definitions

Fix d >= 2 and a Hermitian bond h on C^d ⊗ C^d. Ring Hamiltonian on (C^d)^{⊗L}: H_L = sum_{j in Z_L} h_{j,j+1}, lowest eigenvalue E_L. e_0 = lim E_L / L = infimum of tr(h rho_2) over nearest-neighbour two-site marginals rho_2 of translation-invariant states on Z. M(h) = min over Hermitian X on C^d of (lambda_max - lambda_min)/2 of h + X ⊗ 1 - 1 ⊗ X. Omega_L and Omega_inf are the sets of nearest-neighbour two-site marginals of translation-invariant states on the L-ring and on Z. g_d(L) = sup over h with M(h) = 1 of (L e_0(h) - E_L(h)).
(F1) Trivial bound: L e_0 - E_L <= 2M for all h and L (cut the ring state at one bond and take tensor powers, averaged over translations).
The claim under review: g_d(L) -> 0 as L -> infinity at fixed d, with the explicit rate (1) below.

## (II) The claimed proof, verbatim

[rigorous given the universal recovery theorem stated below] T0 and T1 hold for every fixed \(d\). The proof gives an explicit rate, although the lengths at which it improves on \(g_d(L)\le 2\) are very large. It applies to arbitrary translation-invariant ring states. For permutation-invariant states, a simpler construction gives two-site marginal distance \(O_d(L^{-2})\), including for mixed states and Dicke states. I do not obtain T2 for general bonds.

The construction joins long intervals using a shared classical measurement record. At each join, a channel inserts a shorter interval and reproduces its entire interior marginal exactly. Only the two bonds connecting that interval to the retained sites incur an error. Exact reproduction of the interior matters because the inserted interval can be much longer than the buffers on which the channel acts.

## 1. Statement and rate

[rigorous] Fix an integer \(n\ge 1\). Define \(b_1=1\) and, for \(i=1,\ldots,n\), set
\[
 t_i=b_i+1,\qquad J_i=4n d^{5t_i},\qquad
 b_{i+1}=3t_i\,2^{J_i}.
\]
Then
\[
 \boxed{\quad
 L\ge 2b_{n+1}
 \quad\Longrightarrow\quad
 g_d(L)\le
 \frac{3}{n}+2\sqrt{\frac{2}{n}}+
 8\sqrt{\frac{\ln d}{n}} .
 \quad}                                                     \tag{1}
\]
One can take the minimum with \(2\) from F1. For every fixed \(n\), the length threshold is finite and independent of the bond. Sending \(n\) to infinity proves T1, and homogeneity gives T0, including for bonds that have vacuum-incompatible lengths.

[rigorous] There is a normalization correction in the question. If \(D_L\) denotes the one-sided Hausdorff distance from \(\Omega_L\) to \(\Omega_\infty\), then
\[
 g_d(L)=L D_L .
\]
Indeed, the support-function difference in direction \(-h\) is \(e_0(h)-E_L(h)/L\). Taking its supremum over the dual unit ball gives the directed distance by convex separation. On differences of stationary marginals, the \(M\)-dual norm is the trace norm. Such differences annihilate scalars and telescoping terms, so centering and optimizing the bond gives one inequality. The reverse follows because every bond of operator norm at most one has \(M\le1\). Thus (1) gives the kinematic inclusion in I3 with distance \(o_d(1/L)\).

## 2. Two scales selected from the ring state

[rigorous] Fix any translation-invariant state \(\rho\) on the \(L\)-ring. Write \(S_j\) for the entropy in nats of a \(j\)-site interval and \(u_j=S_{j+1}-S_j\). Strong subadditivity on overlapping intervals gives \(u_{j+1}\le u_j\). Subadditivity and the Araki–Lieb inequality give
\[
 -\ln d\le u_j\le \ln d .
\]
Consequently there is an \(i\in\{1,\ldots,n\}\) with
\[
 u_{b_i}-u_{b_{i+1}}\le \frac{2\ln d}{n}.                 \tag{2}
\]
Otherwise the \(n\) differences would sum to more than \(2\ln d\). Fix this \(i\), and write \(b=b_i,\ b'=b_{i+1},\ t=b+1,\ q=d^t,\ J=J_i\).

[rigorous] We next select a spatial averaging scale for \(t\)-site windows. The following finite informationally complete measurement makes the constants explicit. Take a Hilbert–Schmidt orthonormal basis \(G_\alpha\) of Hermitian \(q\times q\) matrices. Every \(G_\alpha\) has operator norm at most one. Use the effects and associated reconstruction matrices
\[
 F_{\alpha,\pm}=\frac{I\pm G_\alpha}{2q^2},
 \qquad X_{\alpha,\pm}=\pm q^2G_\alpha .
\]
The effects sum to \(I\). For any state \(\gamma\),
\[
 \sum_{\alpha,\pm}\operatorname{tr}(\gamma F_{\alpha,\pm})
 X_{\alpha,\pm}
 =\sum_\alpha\operatorname{tr}(\gamma G_\alpha)G_\alpha
 =\gamma,\qquad \|X_{\alpha,\pm}\|_2^2=q^4.               \tag{3}
\]
The tensor product of two such measurements reconstructs the joint state on two disjoint windows. All expectations below are finite sums with the measurement weights.

[rigorous] For \(R\) consecutive, disjoint \(t\)-site windows, let \(\bar X\) be the arithmetic mean of their reconstruction matrices. Put \(a_R=\mathbb E\|\bar X\|_2^2\). Translation invariance makes \(a_R\) independent of the starting site. If \(U,V\) are two adjacent groups of \(R\) windows each, the parallelogram identity gives
\[
 \mathbb E\|\bar X_U-\bar X_V\|_2^2
 =4(a_R-a_{2R}).                                         \tag{4}
\]
Since \(a_1=q^4\) and \(a_{2^J}\ge0\), there exists \(R\in\{1,2,\ldots,2^{J-1}\}\), a power of two, with
\[
 a_R-a_{2R}\le q^4/J.                                    \tag{5}
\]
Take three consecutive groups \(U,V,W\) of this size. They fit in the ring since \(3Rt\le b'\le L/2\). Equation (4) applies to both adjacent pairs.

[rigorous] Number the windows of \(U\) and \(W\) from \(0\) through \(R-1\), in spatial order. Compare two states on a pair of \(t\)-site windows:
\[
 \tau_{\rm ends}=\frac1{R^2}\sum_{i,j}\rho_{U_iW_j},
 \qquad
 \eta_{\rm ends}=\sum_v p_v\,\gamma_U(v)\otimes\gamma_W(v). \tag{6}
\]
Here \(v\) records the measurement outcomes in \(V\) only, and \(\gamma_U(v)\) is the conditional state averaged over the \(R\) individual windows of \(U\). Define \(\gamma_W(v)\) in the same way. By (3),
\[
 \tau_{\rm ends}=\mathbb E[\bar X_U\otimes\bar X_W],
 \qquad \gamma_U(v)=\mathbb E[\bar X_U\mid v],
 \qquad \gamma_W(v)=\mathbb E[\bar X_W\mid v].
\]
Their difference is therefore
\[
 \mathbb E\big[
   (\bar X_U-\gamma_U(v))\otimes(\bar X_W-\gamma_W(v))
 \big].
\]
Conditional expectation minimizes mean squared Hilbert–Schmidt error among functions of \(v\). In particular, its error for either end is bounded by the error of the estimator \(\bar X_V\), which is a function of \(v\). Applying \(\|Y\|_1\le\sqrt q\,\|Y\|_2\) to each factor, followed by Cauchy–Schwarz and (4)–(5), gives
\[
 \|\tau_{\rm ends}-\eta_{\rm ends}\|_1
 \le 4q(a_R-a_{2R})\le \frac{4q^5}{J}=\frac1n .          \tag{7}
\]
This is where long-range order enters the argument. The two ends can retain order-one correlations through the shared record \(v\).

## 3. The interval to be inserted

[rigorous] For each \(i,j\), take the full interval from the beginning of \(U_i\) to the end of \(W_j\), in the direction passing through \(V\). Its length is
\[
 m_{ij}=(2R+j-i+1)t\le 3Rt\le b'.
\]
Write its first \(t\) sites as \(A_\ell B_\ell\) and its last \(t\) sites as \(B_r A_r\), with \(A_\ell,A_r\) single sites and \(B_\ell,B_r\) of length \(b\). Write \(A=A_\ell A_r\) and \(B=B_\ell B_r\). Form the uniform mixture of the complete interval states, retaining \(i,j\) in a classical register. Call this state \(\tau_{ABC}\), where \(C\) contains the register and all sites between the two \(B\) buffers. Variable lengths are represented by a finite direct sum. In particular, \(\tau_{AB}=\tau_{\rm ends}\) after reordering tensor factors.

[rigorous] Conditional entropy subadditivity and strong subadditivity give
\[
 S(A|B)_\tau
 \le S(A_\ell|B_\ell)_\tau+S(A_r|B_r)_\tau=2u_b .
\]
When \(BC\) is given, the interval label is known. Tracing out the two exterior sites leaves a contiguous interval of length \(m_{ij}-2\). Hence
\[
 \begin{aligned}
 S(A|BC)_\tau
 &=\frac1{R^2}\sum_{i,j}(S_{m_{ij}}-S_{m_{ij}-2})\\
 &=\frac1{R^2}\sum_{i,j}(u_{m_{ij}-1}+u_{m_{ij}-2})
 \ \ge\ 2u_{b'} .
 \end{aligned}
\]
Together with (2), this gives
\[
 I(A:C|B)_\tau\le \frac{4\ln d}{n}.                      \tag{8}
\]

[rigorous given the universal recovery theorem of Junge, Renner, Sutter, Wilde and Winter] For any finite-dimensional tripartite state \(\tau\), there is a channel \(\mathcal R_{B\to BC}\) satisfying
\[
 \mathcal R(\tau_B)=\tau_{BC},\qquad
 -2\ln F\big(\tau_{ABC},(\mathrm{id}_A\otimes\mathcal R)(\tau_{AB})\big)
 \le I(A:C|B)_\tau ,
\]
where \(F(\alpha,\beta)=\|\sqrt\alpha\sqrt\beta\|_1\). This is the relative-entropy recovery theorem applied to the reference \(\tau_A\otimes\tau_{BC}\) and the channel tracing out \(C\). By the fidelity to trace-norm inequality and \(1-e^{-x}\le x\), (8) implies
\[
 \|(\mathrm{id}_A\otimes\mathcal R)(\tau_{AB})-\tau_{ABC}\|_1
 \le 4\sqrt{\ln d/n}.                                     \tag{9}
\]
One may dephase the output length register, since this fixes the target and contracts the error. [Universal recovery maps and approximate sufficiency of quantum relative entropy](https://arxiv.org/abs/1509.07127)

### Correcting the input marginal

[rigorous] The actual input at a join will be \(\eta_{AB}=\eta_{\rm ends}\). Feeding it directly to \(\mathcal R\) would leave a small error on every interior bond. We need a channel on \(B\) that first changes its marginal to \(\tau_B\) exactly.

For states \(r,s\) on \(B\), set \(z=\frac12\|r-s\|_1\). There exists a channel \(\mathcal D\) with \(\mathcal D(r)=s\) such that every extension \(\xi_{AB}\) of \(r\) satisfies
\[
 \|(\mathrm{id}_A\otimes\mathcal D)(\xi)-\xi\|_1
 \le 2\sqrt z+z .                                         \tag{10}
\]
To prove this, let \(T=s+(r-s)_+\), so \(T\ge r,s\) and \(\operatorname{tr}(T-s)=z\). Put \(K=\sqrt s\,T^{-1/2}\) on the support of \(T\), extending by zero. Then \(K^\dagger K\le I\), \(KTK^\dagger=s\), and \(KrK^\dagger\le s\). Complete the operation \(X\mapsto KXK^\dagger\) by preparing the normalized residual \(s-KrK^\dagger\) on the other outcome. A zero residual permits any fallback state. This channel sends \(r\) to \(s\).

The weight of the other outcome is at most
\[
 \operatorname{tr}[(I-K^\dagger K)T]=z.
\]
Also,
\[
 \begin{aligned}
 \|(I-K)\sqrt r\|_2^2
 &\le \|(I-K)\sqrt T\|_2^2\\
 &=\|\sqrt T-\sqrt s\|_2^2
 \le \operatorname{tr}(T-s)=z .
 \end{aligned}
\]
The last inequality follows from \(\sqrt T\ge\sqrt s\). On any extension, the squared Hilbert–Schmidt norm of \((I_A\otimes(I-K))\sqrt\xi\) is the same quantity. Expanding \(K\xi K^\dagger-\xi\) and applying Hilbert–Schmidt Cauchy–Schwarz bounds the trace norm of this difference by \(2\sqrt z\). Adding the other outcome proves (10).

[rigorous] Apply (10) with \(r=\eta_B,\ s=\tau_B\). By (7), \(z\le 1/(2n)\). The join channel \(\mathcal J=\mathcal R\circ\mathcal D\) therefore has both properties
\[
 \mathcal J(\eta_B)=\tau_{BC},
 \qquad
 \|(\mathrm{id}_A\otimes\mathcal J)(\eta_{AB})-\tau_{ABC}\|_1
 \le \zeta_n,                                              \tag{11}
\]
where
\[
 \zeta_n=\frac{3}{2n}+\sqrt{\frac{2}{n}}+
                 4\sqrt{\frac{\ln d}{n}} .
\]
For the inequality, use (10), then (7), then (9), with trace-norm contractivity. The equality in (11) controls the energy of every bond wholly inside the inserted region. The inequality is needed only at its two ends.

## 4. Joining the long intervals and counting the energy

[rigorous] Measure \(V\) on the original ring. Given its outcome \(v\), retain the long interval starting at the beginning of \(W_i\) and proceeding around the other side of the ring to the end of \(U_i\). Average uniformly over \(i\), identifying these intervals in their spatial order. Their common length is
\[
 N=L-(2R-1)t .
\]
Write \(\xi_v\) for the resulting conditional \(N\)-site state. Its rightmost \(t\)-site marginal is \(\gamma_U(v)\), and its leftmost one is \(\gamma_W(v)\). Every internal bond of \(\sum_v p_v\xi_v\) has marginal \(\rho_2\), since the nonselective measurement was outside the retained interval.

[rigorous] For an integer \(k\ge2\), put \(k\) long intervals in cyclic order in the state
\[
 \sum_v p_v\,\xi_v^{\otimes k}.
\]
At every join, the last \(t\) sites on the left and first \(t\) sites on the right have joint state \(\eta_{\rm ends}\) from (6). Apply \(\mathcal J\) to the two \(b\)-site buffers at each join. The channels have disjoint inputs. In particular, summing over the outcomes of the other channels leaves the marginal calculation (11) valid at each join.

[rigorous] In branch \(i,j\) of the output of one join, the number of new sites is
\[
 m_{ij}-2t=(2R+j-i-1)t .
\]
The first equality in (11) gives the prescribed uniform label distribution at each join. Thus the mean number of new sites per join is \((2R-1)t\), and the mean total ring length is \(kL\). Independence of output lengths is unnecessary.

Choose a gauge and scalar shift for which \(\|h\|=M(h)\), and put \(e=\operatorname{tr}(h\rho_2)\). The expected energy can now be counted bond by bond.

- [rigorous] Each retained core has \(N-2b-1\) internal bonds. Their unconditional marginals are still \(\rho_2\), since the channels act outside the core.
- [rigorous] At each join the interior, including both buffers, has unconditional labeled state \(\tau_{BC}\). In each label sector every internal bond has marginal \(\rho_2\). Its expected internal energy is \(e\) times its expected number of internal bonds.
- [rigorous] Two bonds per join connect this interior to the retained core sites \(A_\ell,A_r\). Each has energy at most \(e+M\zeta_n\), by the second part of (11).

These account for every bond of every output ring. The expected total energy is therefore at most \(kLe+2kM\zeta_n\). On the other hand, in every output length sector \(\ell\), F1 bounds the energy below by \(\ell e_0-2M\). Averaging gives
\[
 kLe_0-2M\le kLe+2kM\zeta_n .
\]
Divide by \(k\) and send \(k\) to infinity. We have proved for the original arbitrary ring state that
\[
 L\big(e_0-\operatorname{tr}(h\rho_2)\big)\le 2M\zeta_n .
\]
Taking a translation-invariant ring ground state proves (1).


## (III) Tasks

1. Referee. Check every step: the entropy inequalities and the pigeonhole (2); the informationally complete measurement (3); the parallelogram identity (4) and the selection (5); the covariance identity and bound (7), including the claim that the conditional expectation given the record v on V has mean-square error at most that of the estimator built from V, and that the joint measurement on U, V, W is consistent with the definitions of tau_ends and eta_ends; the conditional-entropy bounds leading to (8), including the variable-length direct sum with a classical register; the precise statement of the Junge-Renner-Sutter-Wilde-Winter universal recovery theorem being used (does the map depend only on tau_{BC}, does it satisfy R(tau_B) = tau_{BC} exactly, is the fidelity convention right, does it hold for states that are not full rank and for direct sums); the marginal-correcting channel (10); the composition (11); and above all the energy accounting in Section 4. In Section 4 check: that the k long intervals with a shared record have the claimed end marginals; that the channels at different joins act on disjoint sites; that every bond of every output ring is counted once; that the two connecting bonds per join have the claimed marginals, taking into account that the site A_l is correlated with its own core and with the record; that the output is a legitimate mixture of states on rings of varying length and that the lower bound l e_0 - 2M may be applied in each length sector with the same gauge-fixed bond; and the limit k -> infinity. Check the length thresholds (3 R t <= b' <= L/2, the retained length N, whether N - 2b >= 2 is needed) and the constants in (1). State clearly: is the theorem correct as stated, correct after repairs (give the repairs), or wrong (give the counterexample or the broken step).
2. Sanity tests. Test the logic on extreme examples: a GHZ-type ring state with long-range order; a ring state that is an equal mixture of the L translates of a product state with a period that does not divide into the window structure; the ground state of a classical frustrated bond on an incompatible length (where E_L - L e_0 is positive and order M); a ring state whose two-site marginal is known to be at distance of order 1/L^2 from Omega_inf. Does the argument behave consistently in each?
3. Rate. The thresholds b_{i+1} = 3 t_i 2^{J_i} with J_i = 4 n d^{5 t_i} make L grow like a tower of exponentials of height about 2n in n. Find the best rate you can prove. In particular: (a) can the informationally complete step be replaced by an entropic argument (for example a chain rule for the mutual information between the end windows and successive measured groups, in the style of Brandao and Harrow) to reduce J; (b) can the dependence on q = d^{b+1} be removed or made polynomial in b, given that only two-site energies are needed at the joins; (c) is there any architecture that avoids tower growth, for example giving g_d(L) <= C(d) / (log log L)^{1/2} or a power of 1/log L. Prove what you can, and state precisely the obstacle to a power law L^{-p}.
4. Corollaries. State and prove the cleanest corollaries: (i) for every fixed h, limsup_L (L e_0 - E_L) <= 0; (ii) the kinematic inclusion Omega_L ⊂ Omega_inf + o(1/L) in trace norm, uniformly at fixed d; (iii) whether the same construction gives sup_h (E_{2L} - 2 E_L)/M -> 0, or a version of it, and if not why not; (iv) a correct statement of what this implies, at physics level, for a family of chains with N infrared sectors of velocity at least eta M and crossover length at most xi: under the usual finite-size form E_L = L e_0 - (pi/6L) sum_{a: xi_a <= L} c_a v_a + ..., give the implied constraint and note that a tower with deficits M N/(L + N^2) is not excluded.

Output: a referee report with a verdict, a list of every gap with its repair, the improved rate if any with a complete proof, and the corollaries. Tag steps [rigorous], [rigorous given cited theorem X] with X stated precisely, or [physics-level]. Plain declarative prose, no em-dashes, no probability estimates.
