# C56. Reference comparison and the counterterm map

[proved, status] This work does not prove C50 (7.1) for the matched action, either for general n or for n=1. The missing estimate is a normalized, marked renormalized-tree bound, not a Gaussian matrix identity. The results below strengthen the endpoint norm, construct the first-order counterterm, and prove continuity of the full counterterm conditional on the collaborator's uniform RG hypotheses. They do not assert a microscopic Luttinger theorem.

[proved, conventions] A proved statement is an implication from its explicitly stated hypotheses. Imported mathematical results are identified. Assumed statements are the remaining application inputs. All limits in the first-order calculation below are on the infinite spatial line. The all-scale conditional theorem uses thermodynamic beta functions and does not assert a new estimate uniform in arbitrary finite-ring aspect ratios.

## 1. Inputs and the reference theory

[assumed, imported model] Use C39/C42's finite spectral construction, C44/C45's rooted logical representation and exact charges, and C48's enlarged-model matching theorem. Put m=n+1. There are 2m chiral Fermi branches, not m chiral branches. Work in retained-bandwidth units W and write their signed velocities as v_a. At each fixed n,

    v_* = min_a |v_a| > 0.

The constants below may depend on v_* inverse, the largest speed, band-isolation data, and the nonresonance margins. No uniformity in n or in a limit v_* -> 0 is asserted.

[proved, inherited geometry] C42 certifies payload resonances. For the wrapper one also needs the mixed margin and the length-number symmetry in C45 section 6 and C46 section 5. Payload U(1) alone would not justify dropping all payload/length conversion terms. With the supplied two exact charges, the surviving local quartics are products of branch densities. Nonzero higher-degree resonances need no Diophantine exclusion for this local classification. Their RG summability is a separate assertion.

Choose a product P of closed intervals for trial roots p'. Shrink it once so that the roots stay separated, each bare dispersion e_a obeys

    |e'_a| >= v_*/2,

and all relevant quadratic and quartic momentum margins persist. If each root moves by at most r, a quartic mismatch changes by at most 4r. Thus r<=delta/8 preserves at least delta/2. Shrinking further gives one common matching energy E_d>0 and common smooth patches throughout P. E_d can be extremely small. It is fixed before the infrared limit.

[read, primary-source input] Mastropietro and Porta, [Theorem 8.5 and section 10](https://arxiv.org/html/2108.13708v3), supply the multichannel reference cancellation and describe its transfer using renormalized trees. The theorem bounds constant-coupling reference beta functions; it does not state the requested replacement-kernel theorem. Porta, Scola and Singh, [section 4.5, equations (4.111)-(4.119), Remarks 4.12-4.13](https://arxiv.org/html/2509.08665v1), retain field-normalization ratios in the tree bounds and describe short memory and covariance continuity. Their running-flow proposition is for their microscopic interaction class. These are imported source results, not independently reproved here.

[proved, use of the requested review] The collaborator's [MULTIVELOCITY_REVIEW.md](../../contrib/decisive_attack_20261028/exact/MULTIVELOCITY_REVIEW.md) already constructs the local Gaussian modes. In its notation D=(4 pi |V|)^(-1), and

    ||D^(1/2) Lambda D^(1/2)|| < 1

ensures a positive quadratic Hamiltonian and a real invertible current matrix of the original signature. For example, ||Lambda||<4 pi v_* is sufficient in the operator norm. The limiting velocities are eigenvalues of KH, not necessarily entries of V. This algebra is used without rederivation. It does not control a microscopic beta function or eliminate the anomalous wave-function factors Z_h. There is no inverse difference |v_a-v_b| in this stability condition.

## 2. Diameter bounds and a stronger C48 norm

### 2.1 A limitation of diameter bounds alone

[proved, norm counterexample] Fix positive a,b,q and field radii H,h. There is no constant C for which every all-degree coefficient measure satisfies

    ||K||_(h,polynomial terminal tree,q)
       <= C ||K||_(H,exp(a D_x+b D_t)).                    (2.1)

Here both norms use the same factorial convention. The left weight on a set S of external positions is

    w_q(S)=min_(terminal spanning trees T)
                   product_(ij in T)(1+|x_i-x_j|+|t_i-t_j|)^q.

Proof. For each integer s>=2 use the s^2 points

    (x,t)=(s i,s j),    0<=i,j<s,

and put one density psi^+ psi^- at each point, using just one species. The field degree is p=2s^2. Translate all points together to make a translation-invariant kernel. Its relative coefficient measure is a finite sum of point masses. Every distinct pair of positions has distance at least s. A terminal tree, after contracting the coincident creation/annihilation positions, needs at least s^2-1 nonzero edges. A nearest-neighbor grid tree attains the lower bound, so

    w_q=(1+s)^[q(s^2-1)],
    D_x=D_t=s(s-1).

Give this monomial coefficient

    c_s=epsilon 2^(-s) H^(-2s^2)
                         exp[-(a+b)s(s-1)].               (2.2)

Antisymmetrizing the kernel leaves this explicit coefficient-norm upper bound unchanged when the accompanying 1/p! is included. The right norm sums to at most epsilon. The s-th contribution to the left norm is

    epsilon 2^(-s) (h/H)^(2s^2)
       exp[-(a+b)s(s-1)] (1+s)^[q(s^2-1)].                (2.3)

Its logarithm divided by s^2 grows as q log s plus bounded terms. The contributions are unbounded. Time-reflect and adjoint-symmetrize each term if a reality relation is desired; this only changes fixed factors. Every term conserves the species charge. It therefore also respects both wrapper charges if it uses only the payload species. This proves (2.1) fails even at n=1.

[proved, scope of the counterexample] This is a counterexample to a norm implication, not to (7.1), to Gaussian stability, or to the actual C48 action. The example need not have the microscopic order-versus-degree structure of that action. The connected-support information in C48 is stronger than the displayed diameter norm and can be retained, as follows.

### 2.2 Retaining C48's connected network

[proved conditional on C48's polymer construction] Denote C48's support exponent L by L_0 to avoid confusion with the ring size. Keep its alpha>0, tau=T/Delta, and B_*. The same polymer majorant controls an exponential terminal-tree weight in the metric

    d_T(z,z')=a_T |x-x'|+k_T |t-t'|,
    a_T=alpha/8,       k_T=L_0/(128 tau).                  (2.4)

One can retain the same conservative B_* and the same small-input condition. The leakage-linear difference estimate survives with this stronger weight. This is an application of the collaborator's connected-network method in [SINGLE_CHANNEL_RG.md section 5](../../contrib/decisive_attack_20261028/demonstrator/SINGLE_CHANNEL_RG.md) to C48's occupied blocks, not a new massive integration theorem.

Proof. In C48 each node is one cell in one time block of width between tau and 2tau. An active interval occupies every cell between its endpoints. Its external legs lie in these occupied blocks. Passive auxiliary propagation fills intervening occupied nodes. Each polymer is connected by nearest-neighbor grid adjacency.

In a connected logarithmic cluster, incompatible polymers overlap or are adjacent. The union of their node supports is therefore connected. Let N be the sum of their node counts, with multiplicities. The union has at most N distinct nodes. For every node put a vertical segment through its entire block. All external times at that cell and block lie on this segment, regardless of their number. These segments have total temporal length at most 2tau N. Choose a grid spanning tree of the union. A spatial adjacency can be joined by one horizontal segment at a time in the common block; temporal adjacencies already touch. The total horizontal length is at most N-1. The resulting network contains every external point.

Double this connected network, take an Euler tour, and shortcut to the external points using the metric triangle inequality. Removing unnecessary terminal edges gives a terminal spanning tree. Its weighted length is at most

    2 a_T N + 4 k_T tau N.                               (2.5)

Circle metrics obey the same argument. If an interaction is point split at block boundaries, a small enlargement of the occupied intervals and at most a fixed number of neighboring block layers gives the safe bound

    T_terminal <= (8 a_T+16 k_T tau)N
               = (alpha+L_0/8)N.                         (2.6)

Alternatively take the ordered zero-splitting limit, where (2.5) suffices. The number of external legs never multiplies the right side: one vertical segment covers all times in a cell's occupied block.

C48's equation (15) already budgets exp[(alpha+L_0/8)N] and retains the root and boundary reserves. Thus its positive cluster majorant also bounds exp(T_terminal). The coefficient estimates and connected-series convergence are unchanged. The analytic difference argument is in the same stronger norm: the K-only zero-temperature functional vanishes, and a Cauchy estimate in the residual insertion gives a bound linear in that residual. The layer argument for finite-temperature winding terms is unchanged. This proves the stated strengthening within the imported polymer representation.

[proved, normalization and crossover] Under t'=Wt the anchored norm changes by W^(-1), and k_T becomes k_T/W. For W<=Delta,

    k_T/W >= L_0/(128T)>0.

There is therefore a positive lower bound on the normalized time rate at fixed n. A fixed loss in field radius converts this exponential terminal-tree weight to any fixed polynomial-tree weight. Indeed, if eta=min(a_T,k_T/W)>0 and

    C_(q,eta)=sup_(r>=0)(1+r)^q exp(-eta r),

then a terminal tree with p legs gives

    w_q <= C_(q,eta)^(p-1) exp(T_terminal).                (2.7)

Choose the input field radius larger by C_(q,eta). C50 section 2A explains how to obtain that radius from C48 by increasing the microscopic support-weight budget before choosing lambda. This does not reweight an already fixed radius-4 bound for free. The collaborator's UV and compact-patch finite integrations can then be used in their stated norms. Their polynomial-tree conversion still assumes the appropriate localization-routing estimate when it is used in an infinite RG expansion.

## 3. What a proof of (7.1) must establish

[proved, sufficient comparison identity] Let beta_h^lat be the normalized marginal beta vector, with the exact carried marginal term removed. It contains velocities and forward couplings, but not Z_h itself. Write its dependence on a marginal history G=(g_(h+1),...,g_0), a relevant history a, and the fixed initial kernels K. Suppose a renormalized-tree construction proves

    |beta_h^lat(G,a;K)-beta_h^ref(G)|
                                      <= A_1 epsilon 2^(theta h),  (3.1)
    |beta_h^ref(u,...,u)|              <= A_2 epsilon 2^(theta h),  (3.2)
    ||partial_(g_k) beta_h^ref||       <= b epsilon 2^[-kappa(k-h)], (3.3)

uniformly on a convex history domain, with 0<theta<kappa. The relevant history in (3.1) must already belong to its proved decaying critical ball, or its contribution must be included explicitly and controlled simultaneously. Equations (3.1) and (3.3) are assumptions here.

Freeze the history at u=g_(h+1). Integrating the derivative of beta_h^ref along the segment from the constant history to G gives

    |beta_h^ref(G)-beta_h^ref(u,...,u)|
       <= b epsilon sum_(k=h+1)^0 2^[-kappa(k-h)] |g_k-g_(h+1)|.

With B_j=|g_j-g_(j+1)|,

    |g_k-g_(h+1)| <= sum_(j=h+1)^(k-1) B_j
                  <= sum_(j=h+1)^k B_j.

Adding (3.1)-(3.2) yields C50 (7.1) with A=A_1+A_2. This algebra identifies the analytic estimates needed; it does not prove them by renaming them.

[assumed, missing tree statement] To establish (3.1)-(3.3), the normalized tree sum must be defined uniformly for the enlarged initial-kernel class. Trees with curvature, retarded remainders or an initial irrelevant endpoint must retain a positive root-to-insertion scale gain after all local subtractions. A history mark must preserve that gain and leave the extra small interaction factor. The dependence of products of Z ratios must be included. The nonconstant-history derivative is not bounded by merely setting Z_h=1.

[proved, relation to prior work] The collaborator's [GN_MARK_CALCULUS.md](../../contrib/direct_attack_20261028/vacancy/GN_MARK_CALCULUS.md) already proves the transfer from suitably marked normalized trees to the absolute and difference beta bounds G1-G2. It explicitly leaves the normalized covariance/tree hypotheses open. Its [DEGREE_REPARAMETERIZATION.md](../../contrib/direct_attack_20261028/vacancy/DEGREE_REPARAMETERIZATION.md) supplies the endpoint hierarchy under a physical perturbative-order/support bound. These results are not repeated here. Section 2 repairs one norm issue but does not prove the missing routing and marked-tree theorem.

[proved, why a finite-scale estimate is insufficient] A finite integration theorem spends a fixed field radius in each use. Applying that theorem at infinitely many scales would require an unbounded initial radius unless a renormalized dimensional estimate replaces that loss. Similarly, |beta_h|<=C epsilon^2 has divergent sum over h<=0. Neither estimate proves (3.1). Convergence of the exact Gaussian current covariance is unrelated to this summation problem.

## 4. A first-order all-degree counterterm

[proved, proposition] Consider the compact-patch action on the infinite line after a finite initialization, with trial roots p' in P. Keep the bare mismatch e(p')-mu separate from the interaction correction V. Suppose the kernels of V have a uniformly small field-weighted norm epsilon with two relative-coordinate moments, at a sufficiently large fixed field radius. Suppose their p' dependence comes from uniformly bounded smooth spectral projections, centered dispersions and Fermi phases of a fixed physical kernel family, or has the corresponding dominated differentiability property with derivative norm bounded by a fixed constant times epsilon after taking the stated moments. Then the coefficient nu^(1)(p') of the first-order self-energy at the trial roots is a real, small, continuous map. It is locally Lipschitz with constant C epsilon, where C depends on the fixed patches and v_* inverse. This statement concerns the coefficient linear in an artificial multiplier z of V, not an expansion whose higher orders have been shown to converge uniformly in the infrared.

Proof. Use the patch covariance

    C_(a,p')(q)=chi_a(q)/[i q_0+xi_(a,p')(q_1)],
    xi_(a,p')(k)=e_a(p'_a+k)-e_a(p'_a).                   (4.1)

Choose chi with compact frequency and momentum support in the fixed patches, equal to one near zero. It can be chosen independent of p' in centered coordinates. Transversality gives

    |i q_0+xi_(a,p')(q_1)| >= c (|q_0|+v_* |q_1|).

The absolute Fourier integral of C is finite, uniformly on P, since 1/|q| is locally integrable in two dimensions. The direct sum over the finite branch labels has a Gram bound

    |det C(X_i,Y_j)| <= D^r                               (4.2)

for r contracted pairs, with D finite and independent of an infrared cutoff on the infinite line.

For a root shift,

    partial_(p'_a) xi_(a,p')(k)
                  =e'_a(p'_a+k)-e'_a(p'_a)=O(|k|).

Consequently

    |partial_(p'_a) C_(a,p')(q)| <= C/|q|                 (4.3)

near zero, with an integrable majorant. Smooth projector derivatives and cutoff derivatives add integrable terms. Thus C and its root derivative are uniformly bounded as coordinate kernels by their Fourier L1 norms. This uses the centered subtraction in (4.1). Differentiating an uncentered denominator would not give the same estimate.

The derivative at z=0 of the effective interaction is E_C V(psi+zeta). Retain two external fields at a degree-p endpoint and contract all the others. The contractions form a determinant of order (p-2)/2. Let k_p denote the anchored total-variation coefficient norm, including its factorial convention. Choosing the external fields costs at most p^2. With R=max(1,sqrt(D)), the absolute quadratic coefficient is bounded by

    sum_(p even>=2) p^2 R^p k_p <= C_H epsilon,           (4.4)

when H>R and sum H^p k_p<=epsilon. Evaluating its Fourier kernel at zero external deviations cannot increase this bound. This is an explicit formula for nu^(1): sum the two-external-leg projections of E_C V and include the first-order finite-initialization quadratic coefficient, using the convention G^(-1)=G_0^(-1)+Sigma. All even degrees contribute; a quartic-only Hartree-Fock formula would omit terms.

A derivative of an r-by-r determinant is a sum of r^2 cofactors times differentiated entries. Equation (4.2) bounds the cofactors and (4.3) bounds the differentiated entries. This costs at most a fixed power of p times a geometric factor in p. Differentiating the Fermi phases costs at most p times the relative diameter. The moment bound and a fixed further field-radius reserve absorb both costs.

There is no uncontrolled absolute anchor in this one-endpoint calculation. A diagonal quadratic root is branch neutral. Since the covariance preserves each branch, a single endpoint contributing to that root must itself have zero net charge in every branch. Its Fermi phase therefore depends only on relative coordinates. This observation is specific to this coefficient; it does not assert anchor-supremum differentiability of arbitrary branch-changing initial kernels.

Dominated differentiation now proves continuity and the stated Lipschitz estimate. If the action and covariance obey the physical adjoint/time-reflection relation with a consistent equal-time prescription, the quadratic coefficient obeys Sigma(iq_0,k)^*=Sigma(-iq_0,k). Its diagonal zero-frequency value is real. This proves the proposition.

[proved, first-order root equation only] For small epsilon, C50's inverse-dispersion argument applies to nu^(1), giving a solution of

    e_a(p'_a)-mu+nu_a^(1)(p')=0,
    |p'_a-p_a| <= 2 C_H epsilon/v_* .                    (4.5)

If its Lipschitz constant is below v_*/2 this truncated equation has a unique solution. No uniform O(epsilon^2) bound on the error between this solution and the true interacting roots is established here. Such a statement would require the higher-order infrared construction. At epsilon=0, (4.5) returns the bare roots.

## 5. Continuity of the full counterterm from a unique bounded RG trajectory

[assumed, exact hypotheses] The following theorem assumes the collaborator's [COUNTERTERM_INDUCTION.md](../../contrib/direct_attack_20261028/vacancy/COUNTERTERM_INDUCTION.md), G1-G2, uniformly over p' in compact P. Assume also fixed-scale continuity of the causal, reduced normalized beta maps and of the initial marginal data. This is the continuity condition H2 in [FERMI_SHIFT_EXISTENCE.md](../../contrib/direct_attack_20261028/vacancy/FERMI_SHIFT_EXISTENCE.md). Its branch-neutral-root domination criterion is applicable; differentiability in p' is unnecessary. These beta hypotheses remain unproved for the full action.

Use a_h for normalized relevant coefficients, to distinguish them from the final physical self-energy convention nu. The imported recursion theorem gives a unique trajectory in a uniform ball, with marginal increments b_h,

    a_(h-1)=2a_h+r_h,       g_(h-1)=g_h+b_h,
    |a_h|<=A epsilon 2^(theta h),
    |b_h|<=B epsilon 2^(theta h),
    |r_h|<=R epsilon 2^(theta h).                         (5.1)

It also gives the finite-cutoff solutions with a_(-N-1)=0 and the same bounds. Only uniqueness within this specified small ball is used. We do not assume uniqueness of physical Fermi points.

[proved, continuity theorem] Under these hypotheses the infinite-scale map p' -> a_0(p') is continuous. The continuous finite-cutoff maps a_0^(N)(p') converge to it uniformly on P. More generally, for every 0<=theta'<theta the finite-cutoff relevant and increment trajectories, extended by zero below their last scale, converge uniformly in p' in the theta'-weighted sequence norm. A proof in the stronger theta norm is neither needed nor claimed.

Proof, existence and bounded tail. Existence, uniqueness and (5.1) are the imported recursion theorem. Iterating the relevant recursion and using its decaying boundary condition gives

    a_h(p')=-sum_(j=-infinity)^h 2^(j-h-1) r_j(p'),
    a_0(p')=-sum_(j<=0) 2^(j-1) r_j(p').                  (5.2)

The j<-M contribution to the second sum is bounded by

    (R epsilon/2) 2^[-(1+theta)(M+1)]
                          /[1-2^(-(1+theta))].           (5.3)

This is an explicit physical-energy-scale tail, not an assumption of a differentiable critical self-energy.

Proof, continuity. Let p'_l -> p' in P and consider the infinite trajectories at p'_l. At each fixed h the vectors a_h,b_h,g_h belong to compact finite-dimensional sets by (5.1) and summability of b. A diagonal subsequence converges at every fixed integer scale. The causal equations at scale h use only finitely many higher-scale coordinates. Fixed-scale continuity therefore passes every recursion to the limit, with initial datum g_0(p'). Bounds (5.1) pass to the limit as well. Backward iteration and the decaying bound recover (5.2), so this limit is in the uniqueness class of the imported theorem. It must be the trajectory at p'. Every convergent diagonal subsequence has that same limit. If a fixed coordinate failed to converge, a subsequence staying a fixed positive distance away would have a diagonal subsequence with that limit, a contradiction. In particular a_0(p'_l)->a_0(p'). This proves continuity.

Proof, uniform cutoff removal. Suppose a fixed-coordinate convergence were not uniform in p'. There would be a positive d, cutoffs N_l -> infinity and p'_l in P with a coordinate error at least d. Pass to p'_l -> p'. Use (5.1) to extract a diagonal limit of the finite-cutoff trajectories. For every fixed h the terminal boundary has moved below h, so the limit obeys the infinite causal equations and the same bounds. It is the unique trajectory at p'. The preceding continuity proof also sends the infinite trajectories at p'_l to that trajectory. Their coordinate error tends to zero, contradicting d. This proves uniform convergence at every fixed scale, including scale zero. The same argument applies to b_h.

Finally, below h=-M the difference of extended finite and infinite relevant trajectories obeys

    sup_(h<-M) 2^(-theta' h)|a_h^(N)-a_h|
       <= 2 A epsilon 2^[-(theta-theta')(M+1)].           (5.4)

The corresponding b estimate uses B. Above -M there are finitely many coordinates and their convergence is uniform. First choose M, then N. This proves the weaker weighted-norm convergence. Reality follows because the finite RG maps preserve the real subspace and their unique fixed points are real; the limits stay real. QED.

[proved, improvement and attribution] FERMI_SHIFT_EXISTENCE already obtains physical roots by finite-cutoff Brouwer and diagonal compactness. The theorem here additionally produces a continuous infinite-scale counterterm map, using the uniqueness in COUNTERTERM_INDUCTION. It does not redo the conditional recursion theorem, and it does not remove its G1-G2 assumptions.

## 6. The physical map and Brouwer, with the matching convention retained

[assumed, exact physical initialization] Keep the original Hamiltonian, chemical potentials and matched kernels fixed. Recenter at p'. Before imposing an artificial relevant coefficient, suppose the finite initialization gives the exact physical relevant vector

    a_0^phys(p')=A(p') d(p')+u(p'),
    d_a(p')=e_a(p'_a)-mu,                                (6.1)

with A,A^(-1),u continuous, ||A^(-1)||<=M_A and ||u||<=C_u epsilon. All other initial data are held at their physical p'-dependent values while only the local relevant coefficient is varied in the auxiliary RG problem. This is H3 in the collaborator's Fermi-shift note. It must be checked in an actual initialization; its validity is not implied by a bound on a_0 alone.

[proved conditional construction] Define

    nu(p')=A(p')^(-1)[u(p')-a_0(p')].                     (6.2)

By section 5 this is real and continuous, and

    ||nu(p')|| <= M_A(C_u+A) epsilon = C_nu epsilon.      (6.3)

Equality a_0^phys=a_0 is equivalent to d+nu=0, with the sign in C50 section 8. The original physical couplings have not been changed. In the simplest A=I,u=0 convention, (5.2) gives

    nu(p')=sum_(j<=0) 2^(j-1) r_j(p').                    (6.4)

An artificial relevant parameter at each branch is only a way to solve these equations. It is not an assumption of independently tunable microscopic chemical potentials.

For completeness, let I_a=[p_a-r,p_a+r] be the intervals defining P. Monotonicity and |e'_a|>=v_*/2 imply e_a(I_a) contains [mu-v_*r/2,mu+v_*r/2]. If

    2 C_nu epsilon/v_* < r,

then

    F_a(p')=e_a^(-1)(mu-nu_a(p'))                         (6.5)

is a continuous self-map of P. Brouwer gives a fixed point. At any fixed point,

    |p'_a-p_a| <= 2 C_nu epsilon/v_* .                    (6.6)

Taking this bound below delta/8 retains all the reserved quartic margins. No assertion of uniqueness follows from continuity. This completes the counterterm map and Brouwer step conditional on the analytic RG and initialization inputs.

## 7. The n=1 case and what changes with n

[proved, operator count] At n=1 there are four chiral fields, payload R,L and length R,L, and four real relevant bookkeeping coefficients. There are at most six independent off-diagonal entries in the symmetric density-coupling matrix. No local same-branch quartic survives Grassmann antisymmetry. The two within-fluid forward terms and the four cross-fluid density terms all conserve the two exact charges. The mixed margin excludes the charge-preserving exchange assignments that would otherwise survive. A pure half-filled length umklapp, if that filling is retained, has repeated same-chirality fields, so its contact quartic vanishes and its two Pauli derivatives must be counted. One must not call that case an unrestricted positive pair-sum margin. One can instead choose a generic permitted length filling.

[proved, unchanged obstruction] Fixing n=1 makes every branch sum finite and reduces the counterterm map to four variables, but it does not bound the time-separated external degree of the matched action. Even a single species admits arbitrarily many legs at distinct times. Thus the norm and repeated-localization problem remains, as section 2.1 demonstrates. The marked-tree theorem is not proved here for this case either.

[proved, general-n changes] Replace four branches by 2(n+1), four relevant coordinates by 2(n+1), and six possible forward entries by binomial(2(n+1),2). The velocity, patch, kernel and field-radius constants change with the fixed reference data. The compactness/continuity proof and the Brouwer construction are otherwise unchanged. C42's payload certificate must always be combined with the wrapper's mixed margin and exact length symmetry. There is no need for a uniform-in-n constant for the stated fixed-n objective.

[proved, final boundary] The requested estimate (a) remains unresolved. The new results are the explicit stronger C48 endpoint weight, the failure of a diameter-only norm shortcut, the first-order all-degree counterterm, and a continuous full counterterm map conditional on the already stated RG hypotheses. The outstanding proof is (3.1)-(3.3), or the collaborator's sufficient normalized marked-tree bounds yielding G1-G2, for the actual retarded endpoints on one common history domain. No thermal coefficient, full source scaling limit, compactification lattice or physical finite-size spectrum is established here.

[numerical, checks] The script checks.py verifies the connected-network bound on 120 finite connected supports with multiple external times per node, the grid spanning-tree formula in section 2.1, a solvable relevant recursion against its closed form, and a four-branch root equation against an independent linear-algebra solution in sine coordinates. All passed. The recurrence formula discrepancy is 8.7e-19 and the root discrepancy is 6.8e-21. These controls do not verify G1-G2 or (3.1)-(3.3). No independent referee has checked the new proofs.
