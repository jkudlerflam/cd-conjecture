# C58. Transfer to the original physical Form W ring

## 1. Statement, conventions and provenance

[proved, scope] This note closes estimate (c) in its thermodynamic version, conditional on the imported finite construction and C48 massive matching theorem. The physical bond is the original C44 wrapped C39 bond, with dimension d=4464. The positive completion in the collaborator's notes is used as an intermediate comparison Hamiltonian. It is not silently substituted for the original bond. The estimates concern free energy, with L tending to infinity before T tends to zero. They do not prove o(1/L) spectral matching.

[proved, notation] Fix n and all the earlier C42/C44/C45 parameters. Write b=n+3, N=n+1, d_P=2^N, and ell_bar=b+1/2. Let H^o_L be the original periodic physical Hamiltonian and H^c_L its positive pin completion. The fixed coexistence energy e is negative. Superscripts o and c mean original and completed. All rooted Hamiltonians below are centered by subtracting e times total physical length. The weak marker chemical adjustment already chosen in C45 is included in the Hamiltonian. An additional marker source is denoted mu. Positive s penalizes physical length. The original core-scale parameter called s in C45 will be denoted s_core here.

The exact rooted pressures per cell are

    p_x(T,s)=lim_(M->infinity) (T/M)
       log Tr exp[-(H^x_root,M-e L_M+s L_M)/T],
    x=o,c.

The physical marker source changes this to p_x(T,s)+mu. Define p_m(T,s) by the exact matched thermal logical action in Section 6, including its scalar normalization. It is an enlarged-Hamiltonian pressure, so it is real and convex for real s. It is not the pressure obtained by dropping the massive thermal endpoint.

[proved, main transfer] At fixed sufficiently small positive couplings there are epsilon>0, T_0>0 and fixed C,delta>0 such that

    sup_(|s|<=epsilon) |p_o(T,s)-p_m(T,s)|
                         <= C T exp(-delta/T),  0<T<T_0.       (1)

Both physical and matched pressure zeros lie in this interval. The physical free energy with marker number summed is

    f_o(T)=e-s_o(T),             p_o(T,s_o(T))=0,

and its matched counterpart is e-s_m(T), with p_m(T,s_m(T))=0. Thus

    |f_o(T)-(e-s_m(T))|<=C T exp(-delta/T).                    (2)

For a fixed marker-density sector M/L tending to rho>0 whose length saddle stays inside the interval, put

    f_m(T,rho)=e+sup_s {-rho p_m(T,s)-s}.

Then the physical fixed-marker free energy exists and

    |f_o(T,rho)-f_m(T,rho)|<=rho C T exp(-delta/T).             (3)

The supremum in (3) can be taken on any compact subinterval containing the saddles. Section 8 supplies the canonical argument. No fixed-L coefficient is extracted from an inequality between fugacity sums. A marker number literally held constant while L grows is not the selected critical regime. The fixed-marker statement means a canonical sequence M_L with M_L/L tending to an interior selected density.

[imported, finite results] The finite core, complete gaps, two selected lengths, positive readout coefficient, compensated ports, exact charges and marker transport are C39/C42 and C44 PROOFS Sections 1-8. C45 PROOF Section 3 supplies the finite local unitary. C48 PROOF Sections 2-4 supplies operator-valued massive matching. C49 PROOF Section 4 supplies the physical-support lift for its generators. These results are used as proved local inputs, not independently rederived here.

[imported, collaborator results] The all-string completion and token-count split are in `../../contrib/thermal_ensemble_20261028/ensemble/ILLEGAL_SECTORS.md`, independently audited in `exact/COMPLETION_AUDIT.md`. The pressure-zero and tilted reference inequalities are in `exact/LENGTH_FUGACITY.md`. The embedding, separate fermionization and application are in `TILTED_EMBEDDING.md`, `exact/SEPARATE_JW.md` and `MOBILE_MATCHING_APPLICATION.md`. Complex source derivatives are in `../../contrib/direct_attack_20261028/infrared/TILTED_MATCHING.md`. Their proof hypotheses are checked below. The new ingredients here are the removal of the completion, fixed-marker pinning, canonical extraction and their assembly.

## 2. Fixed-marker pinning without a positivity assumption

[proved, lemma] Let H=sum_i h_(i,i+1) be a translation-invariant ring, ||h||<=J, with conserved M=sum_i p_i, where p_i is the one-site marker projection. Restrict to the sector with integer 1<=M<=L. Write Z for its partition function and Z_pin for the trace of p_1 H p_1 on its range in this sector. Then

    (L/M) Z_pin <= Z <= L exp(4 beta J L/M) Z_pin.             (4)

There is no assumption that p_1 commutes with H. In particular (4) is a boundary and polynomial comparison when M/L stays positive.

Proof. Scalar spectral convexity in an eigenbasis of p_i H p_i gives

    Tr_(p_i) exp(-beta p_i H p_i)
                  <= Tr p_i exp(-beta H).

Translation symmetry and sum_i p_i=M give the first inequality in (4). For the other direction define the isometry

    V psi=M^(-1/2) direct_sum_i p_i psi,
    K=direct_sum_i (p_i H p_i)|_(ran p_i).

Its compression is

    V* K V=(1/M) sum_i p_i H p_i
       =H-(1/(2M)) sum_i [p_i,[p_i,H]].                       (5)

Only two bonds meet p_i, and ||[p_i,[p_i,h]]||<=4||h||. Consequently ||V*KV-H||<=4JL/M. Scalar convexity again, now in an eigenbasis of V*KV, proves

    Tr exp(-beta V*KV)<=Tr exp(-beta K)=L Z_pin.

Min-max comparison with H proves (4). This uses trace convexity, not operator monotonicity of exp. All spaces in this lemma are finite dimensional. The argument works for the original and completed bonds, and with a marker chemical term, which commutes with every p_i.

[proved, rooted coordinates and seam] Here p=|#><#| tensor I_1116. The pinned fixed-M space is exactly

    K_(M,L) subset K^tensor M,
    K_ell=C^1116 tensor (C^3348)^tensor(ell-1),
    sum_j ell_j=L.

This is a basis permutation and includes all internal marker states and all illegal words. Its dimension is binom(L-1,M-1)1116^M 3348^(L-M). The collaborator's PINNED_MARKER Section 3 proves that restoring motion of the root changes this operator by norm at most B independent of M,L and cell lengths. Each physical boundary matrix unit becomes a bounded partial isometry between two cells. One may use B=2d^4J, or the much smaller sum of actual root-crossing terms. Thus if Z_root,M,L denotes the periodic rooted trace,

    (L/M)e^(-beta B) Z_root,M,L <= Z_phys,L,M
       <= L e^[beta(B+4JL/M)] Z_root,M,L.                    (6)

The error per physical site is bounded by

    (B+4J/rho_L)/L + T log L/L,       rho_L=M/L.              (7)

It vanishes for every fixed T when rho_L stays positive. This handles the actual cyclic identification of roots, including noncommon length-dependent bonds. It does not replace the physical Hilbert space by a momentum-zero projection or impose independent periodic boundary conditions on all fermion species.

[imported and checked, unrestricted marker number] If marker number is summed, PINNED_MARKER (7) gives the stronger constant-factor comparison with sum_M Z_root,M,L, with d/q=4. Its hypotheses hold because every literal term preserves the number of # letters and the physical bond is bounded. The marker-free sector is included in that comparison. Separately, its energy is nonnegative, since ports and marker motion vanish and the readout is positive, while the selected phase retains negative physical energy density for the C45 small-coupling schedule.

## 3. Removing the positive completion for the original bond

[proved, issue] The original endpoint term is conditioned on X# or R#. The incoming A-track pin of a cell therefore depends on the final selector letter of its predecessor on malformed strings. The original all-string unperturbed Hamiltonian is not the product onsite Hamiltonian assumed in the collaborator's completed application. Positivity of the added pins alone gives only one inequality between partition functions. A second estimate is needed.

[proved, invariant cell flags] Let z_j count the selector mismatches inside ordered cell j and on its outgoing bond into the next marker. Let Q_bad,j=1_(z_j>=1). These commute with the full original and completed rooted Hamiltonians. The selector arrows are #:(1,0), L:(0,0), X:(0,1), R:(1,1). Head hopping LX<->XR changes neither outer arrow and has a legal internal bond. It therefore preserves z_j. For R#->#L, the old possible mismatch aR becomes a#, with the same incoming arrow 1; the old #c mismatch becomes Lc, with the same outgoing arrow 0. The new #L bond is legal. The left and right ordered cells retain their respective mismatch counts. The reverse move is identical. Ports, readout, filters, helix operations and wrapper rewrites do not change selector letters. Ordered markers do not cross. The same argument applies across the rooted seam.

Let S be the invariant set of bad cells and D=S union (S+1), with cyclic indices. Then |D|<=2|S|. Every cell outside D has both original endpoint pins, independently of its neighbors. Dropping both endpoint helix pins on D produces a lower operator bound for the original unperturbed Hamiltonian. No large pin norm is charged to the small perturbation.

[proved from imported local lower bounds] On a bad cell with both pins dropped, C44's fragment argument gives, before subtracting e ell,

    A_raw,bad^- >= mu_sel+K_sel z-(2t_sel+4U)(z+1)
                    >= U,        z>=1,

since K_sel=4t_sel+9U and mu_sel>0. Therefore

    A_bad,ell^- >= U+|e|ell.                               (8)

On a selector-legal cell, dropping endpoint helix pins cannot lower the minimum of the other tensor factors. The helices are positive and possess a zero vector even without pins. All internal bonds of a legal cell have the same selector enablement, so this statement is valid as an operator lower bound, including superpositions of selector-head positions. Thus

    A_good,ell^- >=0                         at ell=b,b+1,
    A_good,ell^- >=gamma(1+|ell-b|)           otherwise,      (9)

with gamma=t_sel eta_b-U>0. There is no selector-legal length-one cell. The selected zero space in (9) may have additional helix degeneracy. We count its whole finite dimension rather than assert a gap there. Equations (8)-(9) remain valid with malformed core or wrapper states because the imported open-wrapper bound includes them.

[proved, separating the unbounded readout] Use the collaborator's conserved token count N_tok,j. The code has one token. Ports change its label, and marker motion transfers only a hole, so N_tok,j is invariant also in malformed words. The positive readout splits into a bounded part on N_tok<=1 and a positive onsite part C_j on N_tok>=2. Put lambda^2 C_j into the onsite operator. All remaining perturbations have an interaction incidence norm

    w=sup_j sum_(X contains j)||W_X||=O(lambda)+O(lambda^2).

This is before SW, and W includes the full port, marker, good readout and weak marker chemical adjustment. Its constants include all compensators and overlap inverses. In the comparison below C_j is retained outside D and dropped on D, which only lowers the Hamiltonian.

[proved, heat bounds for damaged cells] Choose epsilon satisfying

    epsilon<=min(gamma/2, |e|/2, U/(8 ell_bar), U/4),
    w<=U/32.                                                (10)

Use the centered tilt s(ell-ell_bar); the common scalar s ell_bar M cancels from all trace ratios. For |s|<=epsilon, (8) gives

    A_bad^-+s(ell-ell_bar)>=7U/8+(|e|/2)ell.

Since dim K_ell<=d^ell, for beta>=2log(2d)/|e| its full trace z_bad satisfies

    z_bad<=exp(-7beta U/8)<=exp(-beta U/2).                  (11)

For legal damaged cells, selected lengths contribute at most
(d^b+d^(b+1))exp(beta epsilon/2). For all other lengths,
|ell-ell_bar|<=1+|ell-b|, so (9) leaves at least
(gamma/2)(1+|ell-b|). The elementary geometric tail used in C46 and LENGTH_FUGACITY then gives, for beta>=2log(2d)/gamma,

    z_good^- <= C_g exp(beta epsilon/2),
    C_g=d^b+d^(b+1)+b d^(b-1)+2d^(b+2).                     (12)

These traces are independent of any maximum cell length. Positive bad-token readouts can be discarded for these upper bounds.

[proved, comparison to the full selector-legal trace] Let Z_g be the full interacting rooted trace with every selector cell legal, including all legal excited states and all lengths. It is common to H^o and H^c. Fix a bad set S, let k=|S| and d_D=|D|, and delete every weak interaction touching D. Their norm is at most w d_D. Outside D the remaining Hamiltonian is the same legal Hamiltonian in the numerator and denominator. Call its trace Z_rest.

The lower operator bound just proved and min-max give

    Z_o,S <= exp(beta w d_D) z_bad^k
                            (z_good^-)^(d_D-k) Z_rest.       (13)

For a lower bound on Z_g, compress the cells in D to the complete code and take their normalized partial trace. The onsite operator has zero code energy, the centered code length has mean zero, and the averaged weak touched interactions have norm at most w d_D. Scalar convexity in an eigenbasis of the averaged Hamiltonian gives the Peierls/partial-trace inequality

    Z_g >= d_P^d_D exp(-beta w d_D) Z_rest.                  (14)

This does not require the code to be invariant under the interacting Hamiltonian. One can first compress finite-length spaces, apply the finite trace inequality, and pass to the trace-class limit. Combining (11)-(14), and using k<=d_D<=2k, gives

    Z_o,S/Z_g <= [A_def exp(-beta delta_def)]^k,
    A_def=d_P^(-1) max(1,C_g/d_P),
    delta_def=U/2-epsilon/2-4w >= U/4.                       (15)

The only entropy factor in A_def is fixed at fixed n. Sum over the binomial family of bad sets, including the empty set. The original and completed traces obey

    1<=Z_o/Z_g<=[1+A_def exp(-beta delta_def)]^M,
    Z_g<=Z_c<=Z_o.                                          (16)

Consequently, uniformly in the tilt window,

    0<=p_o-p_c<=T log[1+A_def exp(-delta_def/T)].             (17)

This is the missing direction of the completion comparison. It holds at fixed nonzero weak couplings, not just for the unperturbed selector. It does not use an extensive ground-energy error to estimate a dilute defect. The cost is paid separately for every bad cell and its at-most-one damaged successor.

## 4. The completed full-space hypotheses and physical dressing

[imported and checked] For H^c the unconditional endpoint pins produce a genuine onsite operator on the complete all-string K. The collaborator proves a complete code P at lengths b,b+1 and, after the positive readout absorption, a positive operator A_lambda with the same code. Choose

    Delta<=min(complete selected gaps, gamma, U, |e|),
    kappa=Delta/(b+2).

Then the imported legal and illegal estimates imply

    A_lambda>=Delta Q,
    A_lambda|K_ell>=kappa(ell+1)Q_ell,
    Tr_Q exp(-t A_lambda)
       <=d exp(-2t kappa)/(1-d exp(-t kappa)).                (18)

The minimum includes the charge tunneling gap and the core/filter combination, not U alone. Delta is deliberately reduced if necessary. This makes the ratio Delta/kappa=b+2 explicit and avoids importing the old legal-only gamma into an all-string heat trace. At t_0=log(2d)/kappa the last bound is 1/(2d).

The bad-token readout commutes with every code-connected SW generator by two-sided annihilation. Its positivity improves (18); it does not have to commute with the rest of A_lambda. The commuting operator required for the source argument is total cell length, which it does preserve. C45 therefore gives, at any chosen support weight a,

    U(H^c_root-e L_M)U*=A_lambda+B_lambda+R_lambda,
    [B_X,P_X]=0,
    ||B_lambda||_a<=K_B lambda,
    ||R_lambda||_a<=K_R lambda^3,
    ||P_all(A_lambda+B_lambda)P_all||<=M K_P lambda^2.        (19)

The constants include scalar port shifts as well as the narrow nonscalar band. They are fixed before lambda. This is a finite-order exact unitary identity with a remainder, not an exact diagonalization.

[imported and checked, physical support] Each generator term on r cells has P_X on one side and preserves total length on X. Both its input and output therefore have length at most (b+1)r. It vanishes outside these finite length sectors and has physical support at most (b+1)r+1. C49 (15) converts its exponential cell norm to a physical exponential norm at any rate mu_phys with mu_phys(b+1)<a. Combining all marker patterns at each fixed physical length is a direct sum, so no exponential pattern count is introduced. The finite-time unitary is physically quasi-local and preserves total marker number, physical length, payload charge and cyclic gluing. It can be restricted to every finite physical sector before taking its trace. Arbitrary auxiliary operators on r cells do not share this physical support property. The proof uses it for the generators only.

## 5. Length source, two charges and the C48 norm interface

[imported and checked, source embedding] Embed K into F_N graded-tensor E, with E=C Omega direct_sum QK. Code states map to their logical states tensor Omega, and excited states map to logical vacuum tensor their original excited vector. The original length becomes

    ell_ext=b+N_len+q_aux,
    q_aux|QK=ell-b,          q_aux Omega=0.                   (20)

There are two independent integer charges, payload and length minus b per cell. Both are retained in the auxiliary grading. They commute with A_lambda and with the total physical dynamics. In particular |q_aux|<=bQ+A_lambda/kappa. For

    |Re s|<=s_a=min(kappa/4,Delta/(4b)),
    a_s=a_lambda+s q_aux,

one has Re a_s>=A_lambda/2>=Delta Q/2 and

    Tr_Q |exp(-t a_s)|<=Tr_Q exp(-t A_lambda/2).              (21)

Take epsilon smaller than s_a/2 as well as (10). The unbounded source stays in a_s. The scalar sb and bounded logical sN_len stay in the logical Hamiltonian. Since U commutes with total length, neither B nor R acquires a new unbounded perturbation under this tilt.

[imported and checked, derivatives] The complex-tilt proof in TILTED_MATCHING applies because A_lambda and q_aux strongly commute and (18) bounds |q_aux| relative to A_lambda. All polymer estimates are uniform on the strip (21). On an inner strip, Cauchy's formula therefore bounds every fixed derivative of the massive influence kernels, including their scalar density, by k! times a fixed inverse strip-width power. The leakage-dependent part remains O(K_R lambda^3), and derivatives of its thermal periodization error remain exponentially small in 1/T. These are derivatives of the matched endpoint, not an assertion of differentiability of the fully integrated gapless pressure. The latter belongs to (a).

[imported and checked, fermionization] SEPARATE_JW Sections 2-6 uses separate payload and length strings formed from total logical-plus-auxiliary parities, with a Klein factor. Neutral interval operators retain their interval support and operator norm, even when they exchange logical charge with an auxiliary. The retained payload and hard-core length hoppings are both quadratic in this representation. The reference extension E_JW(H_P) includes auxiliary parity factors when written in the original tensor variables. Using a naive tensor-identity extension would give the wrong quadratic reference. With the audited extension,

    K_D,X=J_X B_X J_X*-E_JW(P_X B_X P_X),
    ||K_D,X||<=2||B_X||,

and every coefficient of K_D annihilates the auxiliary vacuum on both sides. This checks C48's crucial cancellation hypothesis. On periodic rings the two Jordan-Wigner maps have their charge-parity seam twists. We first cut the bounded seam and compute the thermodynamic pressure. We never assert equality of finite-ring spectra with an independently periodic fermion trace.

[proved from C48, strengthened norm] For any fixed desired field weight H and output spatial rate alpha, choose the input exponent before lambda, for example

    L_*=8(alpha+log 32+10),
    a>=L_*+N log(8H^2)+1.

The matrix-unit conversion costs at most (8H^2)^(Nr) on r cells, independent of dim E. Thus the massive input norms satisfy k<=2K_B lambda and r<=K_R lambda^3. Equation (18), with t_0=log(2d)/kappa, supplies a uniform coarse-time heat trace for (21) at time 2t_0. C48 (8) can be applied with auxiliary gap Delta/2 and

    C_0=max(1,exp(t_0 Delta)/(2d)),
    tau=T_*/(Delta/2),
    T_*>=max(2L_*,Delta t_0,log C_0+2).

Indeed Tr_Q exp(-t A_lambda/2)<=C_0 exp(-t Delta/2) for t>=2t_0. This uses the full excited space, not a noninvariant length cutoff. C48's constants A_match,B_match are finite and its temporal decay is positive. With the above reduced Delta, Delta t_0=(b+2)log(2d), so the dimensionless heat constants can be chosen at fixed n independently of the small charge gap.

Choose lambda last so that

    2K_B lambda+K_R lambda^3 <= v_*/4,
    ||B_lambda+R_lambda||_0<Delta/8,
    A_match K_R lambda^3 <= eta W,
    W=c_n lambda^2,    c_n=J0/s_core>0,                       (22)

as well as C45's SW inequalities and (10). All right sides before the factors of lambda are fixed. The third condition is lambda<=eta c_n/(A_match K_R). No parameter depends on L or T. The C42/C44 static finite-g and unequal-length errors must already fit eta before lambda is chosen. Their scalar density interactions cannot be removed by a marker chemical adjustment or by reducing lambda.

## 6. Exact thermal logical action and activated embedding error

[imported and checked, C48 conclusion] With the normalization of (19)-(22), the enlarged rooted partition function has an exact logical Grassmann representation S_(beta,s). Keep its free logical determinant, the scalar sbM, all scalar port terms, and the complete auxiliary scalar logarithm. Its zero-temperature massive endpoint differs from the full retained second-order Hamiltonian by at most A_match K_R lambda^3 in the chosen all-degree exponential kernel norm. The finite-temperature massive endpoint differs from its antiperiodized zero-temperature endpoint by at most

    C exp(-delta_w/T),           delta_w>0,                  (23)

including scalar densities and every fixed source derivative on a smaller strip. Time moments follow by reserving exponential weight. This is the massive winding conclusion of C48, not an application of a gapless single-shell winding lemma to the whole RG.

Define Z_m,M(T,s) to be this exact normalized logical partition function. By construction it equals the enlarged trace. Hence its pressure p_m is a physical real convex function of s despite its retarded Grassmann presentation. All logical integrations, including vacuum terms, are still to be performed. Retaining S_(beta,s) avoids assuming that a small endpoint error automatically gives a small fully integrated pressure error.

[imported and checked, invalid-state comparison] The embedding's invalid states are not illegal physical cells. For w_c=||B_lambda+R_lambda||_0, TILTED_EMBEDDING (9)-(11) gives

    1<=Z_m,M/Z_c,M
       <=[1+exp(beta w_c) zhat_I(beta,s)/d_P]^M,
    zhat_I(beta,s)<=(d_P-1)Tr_Q exp(-beta A_lambda/4).

For beta>=4t_0, (18) implies

    0<=p_m-p_c<=T log[1+C_I exp(-delta_I/T)],
    C_I=((d_P-1)/d_P)exp(t_0 Delta)/(2d),
    delta_I=Delta/4-w_c>0.                                  (24)

This uses the tilted partial-trace proof, which remains valid when the centered length source is present. An untilted embedding estimate alone would not suffice. Combining (17) and (24) proves (1), with any convenient larger C and delta=min(delta_def,delta_I). The original and enlarged pressures can lie on either side of one another. Only their absolute difference is asserted.

[proved, source limitation] The uniform real-tilt pressure bounds suffice for pressure zeros and canonical duality. They are not derivative bounds on arbitrary normalized observables. Bounded neutral microscopic sources can be dressed and embedded with C48's source theorem, and the unbounded uniform length source has just been handled. Charged fields acquire their Jordan-Wigner strings. A theorem for arbitrary charged physical source functionals is neither used nor claimed.

## 7. Interior pressure zero and physical free energy

[imported and checked] All rooted interactions have finite cut strength. For the original bond this includes the finite, possibly large conditional pin terms; smallness of those pins is not needed for cutting. Onsite terms and the unbounded length source do not cross cuts. The heat bounds above give a common trace-class interval. LENGTH_FUGACITY Sections 1-2 therefore construct p_o,p_c and identify the physical pressure with their unique interior zero after adding mu. Restoring e gives

    p_phys,o(T,mu)=s_o(T,mu)-e,
    p_o(T,s_o)+mu=0.                                        (25)

Their strict monotonicity needs no compressibility estimate. Every original cell has length at least one, so

    p_x(T,s+t)<=p_x(T,s)-t,        t>0.                      (26)

The enlarged length is also at least one, including invalid states, so (26) holds for p_m.

[proved, placement] The code has b<=ell_P<=b+1. The finite scalar bound in (19) gives

    p_P(T,+epsilon)<=K_P lambda^2+T log d_P-b epsilon,
    p_P(T,-epsilon)>=-K_P lambda^2+b epsilon.

The collaborator's tilted reference comparison gives an exponentially small difference p_D-p_P uniformly on this interval. The bounded R changes the pressure by at most K_R lambda^3. Choose lambda still smaller so that

    K_P lambda^2+K_R lambda^3<=b epsilon/4,

and then T small enough that T log d_P and all activated errors total at most b epsilon/4. This brackets the roots of p_c,p_o,p_m strictly inside (-epsilon,epsilon), with a fixed margin. The same holds for mu in a sufficiently small fixed neighborhood of zero. The norm error K_R lambda^3 is used only for placement. It is never used as the final thermal error.

If |p_o-p_m|<=d_T on an interval containing their zeros, (26) implies |s_o-s_m|<=d_T by evaluating one function at the zero of the other and shifting s by d_T. Equations (1),(25) prove (2). Every seam error has already vanished when L tends to infinity at fixed T. There is no exchange of limits.

[proved, matched renewal interpretation] The matched zero is also the thermodynamic pressure of the actual positive coefficients of the enlarged trace, summed over marker counts at fixed total length. To see existence rather than just a convergence abscissa, set S_L=sum_M Z_m,M,L in an open-cell convention. For fixed total cell count K, the subspaces in which a prefix of M cells has physical length L_1 and its complement has length L_2 are orthogonal for different M: positive cell lengths make a specified prefix sum occur at most once. Peierls compression onto these subspaces and deletion of the crossing interactions give

    S_(L_1+L_2)>=exp(-beta C_cut) S_L_1 S_L_2.

Each fixed-length enlarged cell space is finite dimensional and its dimension grows at most exponentially in length, since its states have original length ell or ell-1 and a finite logical factor. The onsite operator is positive before the tilt, and the bounded perturbation has norm at most a constant times M. Thus S_L has an exponential upper bound. Every length is available. The elementary superadditive blocking argument proves the limit of T log S_L/L. The positive generating identity sum_L exp(-beta sL)S_L=sum_M Z_m,M(s) then identifies this limit with s_m. This is the proof in LENGTH_FUGACITY Section 2.1 extended to the present all-length space; it requires only positive lengths, finite fixed-length dimensions and local length conservation, all checked here. Restoring e gives the matched free energy in (2).

## 8. Canonical length and fixed marker density

[proved, ensemble lemma] Consider a translation-covariant rooted Hamiltonian on a countable direct sum of finite-dimensional length spaces, with local conservation of total length, finite cut strength C_cut, and a trace-class open interval of real tilts. Suppose the negative one-sided slopes at two interior tilts strictly bracket ell. Then the constrained free energy per cell at L_M/M tending to ell exists and equals

    f_cell(T,ell)=sup_t {-p(T,t)-t ell}.                     (27)

For the application it suffices to work in a compact interval of selected densities with this strict bracket. Differentiability of the finite-temperature pressure is unnecessary. Here is a blocking proof that does not assume a local central limit theorem.

Cut a chain of qk+r cells into q blocks of k cells and a remainder, with 0<=r<k. The sum of the deleted interactions has norm at most C_cut(q+1), including a fixed periodic seam correction if needed. This is also a bound after restriction to fixed total length, since all pieces conserve it. Thus constrained logarithmic partition functions differ by at most beta C_cut(q+1).

For the independent k-cell blocks, write w_l=Tr_(L_k=l)exp(-beta H_k). Their tilted sum is finite in an open interval. The coefficient of total length is a positive convolution. Chernoff's elementary bound, obtained by multiplying every term by exp(-beta t L), gives the upper coefficient bound at every t. For a lower bound choose a finitely supported rational frequency distribution r_l of block lengths with mean ell k. The multinomial contribution has, by factorial asymptotics, logarithm per block

    sum_l r_l log w_l-sum_l r_l log r_l.

Finite distributions approximate the tilted distribution w_l exp(-beta t l)/Z_k(t) in this expression and in mean, because an open exponential-moment interval controls both tails. Adjust a small part of a finite distribution on two lengths bracketing the desired mean to make that mean exact. Integer multiplicities approximate the resulting frequencies. Bounded changes in finitely many block lengths make the total length exact along any sequence L/M tending to ell, without changing the rate. In the all-string spaces used here every sufficiently large integer length is available, so no sublattice obstruction occurs. The finite remainder block costs o(M). This proves the product-block canonical rate as the Legendre transform of (T/k)log Z_k(t).

Deleting cuts also gives |p(T,t)-(T/k)log Z_k(t)|<=C_cut/k on compact tilt intervals. Convex finite-difference bounds and the strict slope bracket confine the optimizing tilts of both pressures to a common compact interval for large k. The finite-block pressure is smooth there and its negative slope ranges across ell, so its tilted distribution used above exists. The two Legendre transforms differ by at most C_cut/k. The preceding constrained norm bound gives the same error for the true chain. Let k tend to infinity. This proves (27), existence, and continuity at the interior densities under discussion, including possible finite-temperature nondifferentiability. The proof also covers the enlarged Hamiltonian: its length spectrum is positive integer, every block length at least k is represented, and (21) supplies the exponential moments. No finite dimension for the whole cell was used.

[proved, application to physical fixed M] Apply (27) to the original rooted Hamiltonian. The original conditional pins are bounded intercell terms, and (8)-(9) imply trace-class tails even before completing them, so the lemma applies. Equation (6) identifies its constrained thermodynamic free energy with the actual periodic physical one. Since M/L tends to rho and ell=1/rho,

    f_o(T,rho)=e+sup_s {-rho p_o(T,s)-s}.                    (28)

Apply the same lemma to the enlarged representation to define f_m. If the optimizing tilts of both pressures are in the uniform comparison interval, taking suprema in (1) proves (3). The following finite-difference check verifies this confinement for a nonempty interval of the scheduled densities, without an infrared susceptibility theorem.

[proved from the finite schedule, canonical saddle bracket] Let p_G,0(s) be the ground pressure of the common quadratic payload plus free length reference, with all common scalars and the chosen marker chemical potential retained. In an interior length-filling interval its length contribution has

    -p_G,0'(s)=b+nu_G(s),
    p_G,0''(s)=1/[2 pi W j_len sin(pi nu_G(s))]>0,

where j_len>0 is the fixed normalized hopping. This is C44's already computed strict length convexity, expressed in the dual variable. C42/C44's finite error bounds give ||H_P-H_G||_loc<=epsilon_stat W after their ordered choices. The reference thermal comparison, the R norm and (17),(24) imply, for x=o,m,

    |p_x(T,s)-p_G,0(s)|<=E_T,
    E_T=epsilon_stat W+K_R lambda^3+T log d_P+O(T exp(-c/T)),

uniformly on a selected source interval. This norm error is used only for a saddle bracket, not for a thermal coefficient.

Choose a<s_0<d inside that interval and h>0 with a+h<d-h. For example all distances can be fixed positive multiples of W. Choose ell in a compact interval strictly between

    -[p_G,0(d)-p_G,0(d-h)]/h
       and -[p_G,0(a+h)-p_G,0(a)]/h.

The ordering is strict by the displayed positive curvature. Reserve a positive margin m to both endpoints. Convexity and the uniform pressure error give

    -p_x,+'(T,a)>=-[p_G,0(a+h)-p_G,0(a)]/h-2E_T/h,
    -p_x,-'(T,d)<=-[p_G,0(d)-p_G,0(d-h)]/h+2E_T/h.

Make epsilon_stat small enough in the earlier finite C42/C44 schedule, then lambda small enough, then T small enough, that 2E_T/h<m. This is possible since h is a fixed multiple of W, K_R lambda^3/W=O(lambda), and T tends to zero after lambda is fixed. These inequalities supply the strict slope bracket for both exact pressures. All optimizing tilts lie in (a,d), so (3) follows for this explicit nonempty density interval. At a fixed reference target one chooses a,d,h around its conjugate tilt and allocates the same finite margin. This also shows why no extra physical-ensemble hypothesis is being hidden in (a). The differentiable dressed response in (a) is needed later to identify the coefficient and the density-dependent speeds, not for the activated canonical comparison itself.

[proved, relation to the marker source] The same formula follows from physical canonical/grand-canonical duality: p_o(T,s)+mu=0 sets mu=-p_o(T,s), and sup_mu {mu rho-p_phys(T,mu)} becomes (28). The block proof supplies the existence and equivalence required for this manipulation. A fixed marker number removes one global zero-mode sum. It does not remove the propagating length oscillator or change a thermal coefficient per site in this density regime.

## 9. Final conditional theorem with only (a) and (b) open

[assumed, estimate (a), stated in the form actually needed] For the C48 matched family at the fixed microscopic parameters, the all-scale constructive RG and counterterm map exist on a small compact real interval I of physical length tilts around the selected root. They preserve the C42/C44 nonresonance margins and yield N=n+1 stable modes with positive dressed chiral speeds v_(a,R)(s),v_(a,L)(s), continuous in s. The source part gives a continuously differentiable ground pressure p_0(s), with mean length -p_0'(s) in (b,b+1). For fixed-density statements use a smaller interval on which this mean is strictly varying, with the selected density in the interior of its range. This is the length compressibility/source component of the missing RG theorem. C50's source-free beta-function estimate alone is not being renamed as this stronger conclusion. If (a) is interpreted as source-free only, this source component remains an additional required input and the fixed-density coefficient cannot yet be asserted.

[assumed, estimate (b)] With the same temperature-independent microscopic parameters and the same source interval, the subtracted scalar RG estimate of C51 Section 7 holds uniformly in s for the actual thermal matched action S_(beta,s). All scalar determinants, field-renormalization constants and the exponentially small massive thermal endpoint are included. In particular it yields

    p_m(T,s)=p_0(s)+A(s)T^2+r_T(s),
    sup_(s in I)|r_T(s)|/T^2 ->0,
    A(s)=(pi/12) sum_(a=1)^N [1/v_(a,R)(s)+1/v_(a,L)(s)].   (29)

One sufficient scalar estimate is |delta_h(T,s)|<=C E_h^2 eta_h min(1,(T/E_h)^p), with p>2, eta_h bounded and tending to zero uniformly, together with its UV and actual thermal-history bounds. The collaborator's WINDING_ESTIMATE proves the single-scale fixed-history determinant comparison, including arbitrary field degree. It does not prove the infinite-scale history or common analytic-radius/moment control. Those remain in (a),(b), as also recorded by C55-C57. A claim that winding by itself has proved (29) would be incorrect.

[proved conditional theorem] Assume precisely these two infrared inputs, and the imported finite construction and C48 theorem already used above. For the original d=4464 periodic Form W bond, with all physical sectors retained and marker number summed, let s_* solve p_0(s_*)=0 and ell_*=-p_0'(s_*)>0. Then

    f_o(T)=e-s_*-[A(s_*)/ell_*]T^2+o(T^2).                 (30)

Proof. By (1), p_o has the same uniform expansion (29). Equation (26) gives s_o(T)-s_*=O(T^2). At s_*, differentiability gives p_0(s_o)=-ell_*(s_o-s_*)+o(|s_o-s_*|), while continuity gives A(s_o)=A(s_*)+o(1). Insert these into p_o(T,s_o)=0 and restore e. This proves (30) without differentiating the thermal remainder.

For a fixed density rho in the source interval, let s_rho be the unique solution -p_0'(s_rho)=1/rho. Then

    f_o(T,rho)=e+sup_s {-rho p_0(s)-s}
                           -rho A(s_rho)T^2+o(T^2).         (31)

Indeed strict exposure confines all finite-T optimizers to s_rho+o(1). Evaluate the supremum once at s_rho for a lower bound, and once at a maximizing point for the upper bound. The zero-temperature objective there cannot exceed its value at s_rho, while A differs by o(1). Uniformity of r_T proves (31). No assumed O(T^2) motion of the canonical saddle and no Hessian of the remainder are needed.

In both cases the physical speeds are u_(a,R/L)=ell v_(a,R/L), with ell=ell_* or 1/rho as appropriate. Therefore the coefficient is

    f_o(T)=f_o(0)-(pi/12) sum_(a=1)^N
                  [1/u_(a,R)+1/u_(a,L)]T^2+o(T^2),           (32)

and likewise in each fixed-density sector, with its own ground energy and speeds. At the equilibrium rho=1/ell_*, (30) and (31) agree. Distinct right and left speeds are allowed. No common Lorentz speed or bare-speed substitution is made.

[proved, meaning of the conditional boundary] Estimate (c) adds no new infrared assumption to (a),(b) in the source-uniform form just stated. Its non-RG errors are activated, its cyclic and seam errors vanish first in the thermodynamic limit, and its canonical transformation is proved. If (a) additionally establishes the complete compact Gaussian scaling theory, C51 identifies c_R=c_L=n+1. The free-energy transfer alone does not construct the stress tensor or prove exhaustion by those fields. No unconditional c-d counterexample is claimed.

## 10. Checks, limitations and hypothesis ledger

[proved, consistency checks] At lambda=0 the completion comparison still pays a positive cost per illegal cell, despite unpinned neighboring helix degeneracy. It never subtracts a large J_h pin norm from Delta. At s=0 the tilted embedding reduces to C48's valid/enlarged comparison. At the free common-bond reference, (31) reduces to C46's canonical length-fluid formula and the speed conversion ell=b+nu. The sign of the pressure root is consistent with restoring H=eL+H_centered: p_phys=s_*-e and f_phys=e-s_*. The dimensions of T exp(-delta/T), pressure, and the pinning boundary terms are energies per the stated length unit.

| Imported result | Hypothesis at risk | Check or repair |
|---|---|---|
| C44 selected gaps | Complete gap, not U | Delta in (18), including charge tunneling |
| Completed all-string cell theorem | Original pins depend on predecessor | Completion used internally; (8)-(17) transfer back |
| Full-space bounded perturbation | Multi-token readout unbounded | Positive invariant onsite split |
| C45 SW | Strong support norm and unbounded onsite domain | Exponent chosen before lambda; bounded commutator construction; readout commutes with generators |
| C49 physical lift | Arbitrarily long cells | P_X on one side plus conserved length bounds generator support |
| Separate JW | Two charges, including auxiliary exchange | Total parities and isometric CAR reference extension |
| C48 massive matching | Heat trace, coefficient norm, vacuum-annihilation, small input | (18)-(23), no auxiliary-mode expansion or cutoff |
| Tilted embedding | Unbounded source and invalid states | (20)-(24), uniform real/complex tilt window |
| Pinning and root quotient | Original theorem sums M | New fixed-M isometry (4)-(7) |
| Pressure-zero theorem | Interior sign change | Explicit scalar-inclusive bracket in Section 7 |
| Canonical extraction | A pressure error is not a coefficient bound | Blocking/Legendre proof and exposed-density confinement |
| Thermal winding | Fixed-history versus actual RG history | Massive winding is proved; gapless history remains in (b) |

[proved, spectral limitation] The seam bound is O(1), and the pinning entropy bound is O(log L). Neither is o(1/L) in energy. C49's finite prethermal spectral window also does not extend to every L at fixed coupling. No finite-ring tower multiplicities or Casimir coefficient follow from this note. The thermodynamic alternative requested in C58 is the conclusion proved here.

[numerical, local controls] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes 75 fixed-marker trace comparisons on number-conserving quantum rings, with compression-identity residual 3.29e-16; 6400 selector rewrites preserve each ordered cell's mismatch count; and nine completion comparisons on a noncommuting flag/pin model satisfy the activated bound. The largest matrix in these checks is 64 by 64. These support the algebraic proofs and are not microscopic Form W simulations. No large local computation, cluster job, network shell command, delegated agent or write outside work_C58 is used. The imported results and the new derivations remain without human review.
