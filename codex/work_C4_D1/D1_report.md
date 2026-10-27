D1. Closed-book attack on the c–d conjecture

No in-scope counterexample was established, for either the Virasoro charge or the effective charge. The two explicit non-compact candidates reconstructed below fail Hermiticity in the ordinary tensor-product inner product. This failure survives the periodic sum and, for the staggered chain, blocking. At the parameters tested they also have complex eigenvalues, so they cannot be made Hermitian by any positive-definite finite-size metric.

This is not a proof that no counterexample exists. In particular, the compact-family survey is incomplete as a classification of all anisotropies and regimes. I have not verified the literature identifications externally. The instruction to work closed book was enforced: no preexisting project files, papers, websites, or other agents were consulted. Everything read as task data was created in this run. Installed numerical libraries were used.

The status labels below distinguish four kinds of statements. P means derived here, including algebraic proofs. N means checked numerically, with finite sizes and residuals stated. G means supplied in the question. A means an additional assumption or a recalled identification that was not proved or externally verified. An A attached to a continuum identification applies to its central charge and unitarity as a lattice limit, not to the algebraic evaluation of a displayed WZW formula.

1. What the non-compact loophole would mean

For the proposed Euclidean black-hole spectrum, take k>2 and write the continuous representations in the convention j=1/2+is. The assumed conformal weights are

\[
h(s,m,w)=\frac{s^2+1/4}{k-2}+\frac{(m+kw)^2}{4k},\qquad
\bar h(s,m,w)=\frac{s^2+1/4}{k-2}+\frac{(m-kw)^2}{4k}.
\]

The identification of these representations with either lattice model is A. Conditional on that identification, the following calculation is P. At the bottom of the uncharged, unwound continuum,

\[
h_{\min}=\bar h_{\min}=\frac{1}{4(k-2)},\qquad
c=2+\frac{6}{k-2},\qquad
c_{\rm eff}=c-24h_{\min}=2.
\]

The minimum is an infimum of a continuous spectrum, not necessarily a normalizable eigenstate. The identity operator is not a normalizable vacuum state of this Hilbert space. The scaling dimensions relative to the continuum threshold are

\[
x-x_{\min}=\frac{2s^2}{k-2}+\frac{m^2}{2k}+\frac{k w^2}{2}.
\]

The cylinder energy explains the distinction without an entanglement argument. For one velocity,

\[
E=Le_\infty+\frac{2\pi v}{L}\left(h+\bar h-\frac c{12}\right)+o(L^{-1}),
\]

so the ground-state coefficient is

\[
E_0=Le_\infty-\frac{\pi v}{6L}\bigl(c-12(h_{\min}+\bar h_{\min})\bigr)+o(L^{-1}).
\]

Thus periodic finite-size scaling measures c_eff. For several velocities it measures the sum of v_a c_eff,a, not the sum of the charges alone. The statement that a particular non-compact ground-state entanglement prescription also measures c_eff requires assumptions about the state and regulator. No entanglement calculation was performed here.

For d=3 the hypothetical Virasoro charge exceeds d-1 for every finite k>2. For the blocked d=4 chain it exceeds 3 exactly when 2<k<8. At k=8 it equals 3. In both cases the hypothetical c_eff obeys the bound. These comparisons are P; obtaining such a continuum limit from an in-scope Hamiltonian is the missing requirement.

The conjecture as written uses Virasoro central charges. One should not silently replace them by c_eff. There are two separate possible conjectures: (i) a bound on the sum of Virasoro charges, and (ii) a bound on the sum of effective charges. A normalizable vacuum with h_min=0 makes them coincide. Without it, finite-size and standard ground-state entropy evidence for (ii) does not establish (i).

There is a further limitation to the proposed verification. A background-charge improvement can change a non-compact stress tensor and the continuum threshold together, leaving c_eff unchanged. Observing c_eff=2 and logarithmically dense levels does not by itself determine the Virasoro charge or prove the Euclidean black-hole operator algebra. The works named in the question are therefore leads for a subsequent open-book check, not independent confirmations made in this run.

2. Explicit spin-1 twisted R-matrix candidate

This construction provides all entries of a fixed 9 by 9 nearest-neighbor density without invoking an unspecified transfer matrix. Let q=exp(i theta), 0<theta<pi/6, and define

\[
[a]_q=\frac{\sin(a\theta)}{\sin\theta},\quad
Z=\mathop{\rm diag}(1,0,-1),\quad
E=\sqrt{[2]_q}\begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix},\quad F=E^T,\quad K=q^Z.
\]

On two sites set

\[
E_{12}=E\otimes K+K^{-1}\otimes E,\quad
F_{12}=F\otimes K+K^{-1}\otimes F,\quad
M=Z\otimes1+1\otimes Z,
\]

\[
C=F_{12}E_{12}+[M]_q[M+1]_q,\qquad
\lambda_J=[J]_q[J+1]_q,\quad
P_J=\prod_{J'\ne J}\frac{C-\lambda_{J'}}{\lambda_J-\lambda_{J'}}\quad(J=0,1,2).
\]

Direct symbolic multiplication gives (C-lambda_0)(C-lambda_1)(C-lambda_2)=0 and traces tr(P_0), tr(P_1), tr(P_2)=(1,3,5). The distinct roots in the stated interval imply P_J P_K=delta_(J,K)P_J and sum_J P_J=1. These are algebraic projectors, not orthogonal projectors in the standard positive inner product when q lies on the unit circle. Define

\[
\check R(x)=P_2+
\frac{q^4x-1}{q^4-x}P_1+
\frac{q^6x+1}{q^6+x}P_0.
\]

It is regular at x=1. Its logarithmic derivative in the chosen normalization is

\[
h_{\rm IK}=i\left.x\frac{d\check R}{dx}\right|_{x=1}
=\cot(2\theta)P_1-\tan(3\theta)P_0,
\qquad H_L=\sum_j h_{{\rm IK},j,j+1}.
\]

The identification of this baxterization, sign, and parameter interval with the conventional Izergin–Korepin regime III Hamiltonian is A. The commonly recalled dictionary is gamma=2 theta and k=2pi/gamma=pi/theta, but it was not derived here. Changing the overall sign changes the low-energy problem and does not cure Hermiticity. A different model or boundary prescription must be checked independently rather than being identified by name alone.

Consistency checks N: the projectors sum to the identity to 4.6e-16 and their largest idempotence residual is 2.5e-15 at theta=0.19. The multiplicative braid Yang–Baxter equation at x=1.17, y=0.83 has relative residual 3.6e-16. This is a numerical check of the proposed R-matrix, not a symbolic proof of the full Yang–Baxter identity.

The Hermiticity obstruction is P and does not depend on the model name or continuum identification. In the ordered basis (+,0,-), the two-site matrix is complex symmetric. An entry computed directly from the displayed projectors is

\[
\langle +,-|h_{\rm IK}|0,0\rangle
=\langle0,0|h_{\rm IK}|+,-\rangle
=\frac{iq}{q^6+1}
=\frac{\sin(2\theta)+i\cos(2\theta)}{2\cos(3\theta)}.
\]

Consequently

\[
\langle +,-|(h_{\rm IK}-h_{\rm IK}^\dagger)|0,0\rangle
=i\frac{\cos(2\theta)}{\cos(3\theta)}\ne0.
\]

Take two many-body basis states differing by this change on one bond and identical elsewhere. Both sites change. No other nearest-neighbor term connects the states, for L>=3. The same nonzero entry therefore occurs in H_L-H_L^dagger. This is not a telescoping imaginary boundary field.

At theta=0.19, the full L=3 spectrum contains

\[
E=3.643351653224842+1.059723971336765i
\]

with eigenvector residual 4.8e-15. The L=4 spectrum contains 3.181324822007003-1.156920550606275i with residual 1.1e-14. These are N. If a positive-definite metric eta obeyed H^dagger eta=eta H, then eta^(1/2) H eta^(-1/2) would be Hermitian and all eigenvalues would be real. Thus these particular finite-size Hamiltonians have no such metric. A real low-energy subsector would not establish quasi-Hermiticity of the full chain.

Verdict: d=3, NN, one-site translation invariant, explicitly not Hermitian. Conditional black-hole values are c=2+6/(k-2) and c_eff=2, with an intrinsically unitary coset for the stated k if that identification is correct. Neither lattice unitarity nor the continuum identification follows from the ED. This implementation is out of scope for both versions of the conjecture.

3. Explicit staggered six-vertex / Temperley–Lieb candidate

In the basis (00,01,10,11), put

\[
e=\begin{pmatrix}
0&0&0&0\\
0&e^{-i\gamma}&-1&0\\
0&-1&e^{i\gamma}&0\\
0&0&0&0
\end{pmatrix},\qquad 0<\gamma<\pi/2.
\]

One directly checks e_j^2=2cos(gamma)e_j and e_j e_(j+1)e_j=e_j. The three-site density and Hamiltonian used here are

\[
g_{123}=2\cos\gamma\,e_{12}-e_{12}e_{23}-e_{23}e_{12},\qquad H_L=\sum_jg_{j,j+1,j+2}.
\]

This representative is translation invariant already by one microscopic site. Its relation to the staggered six-vertex construction can be checked without a continuum claim. Set

\[
\check R_6(u)=\sin(\gamma-u)1+\sin u\,e,\qquad R_6(u)=P\check R_6(u),
\]

and form the auxiliary trace with alternating inhomogeneities 0, pi/2. At L=4, gamma=0.4, a direct numerical logarithmic derivative gives

\[
t(0)^{-1}t'(0)+t(\pi/2)^{-1}t'(\pi/2)
=\frac{H_L}{\sin\gamma\cos\gamma}-2L\cot(2\gamma)1.
\]

The relative residual is 4.7e-11 using a centered derivative with step 1e-6. This transfer-matrix check is N; the equality for arbitrary L has not been proved here. Identifying its scaling limit with a particular Ikhlef–Jacobsen–Saleur, Frahm–Seel, or Bazhanov–Kotousov–Koval–Lukyanov continuum prescription remains A. In particular, a common recalled dictionary k=pi/gamma is not independently established by this check.

The bulk Hermiticity question can be settled algebraically. Let X,Y,Z denote Pauli matrices, c=cos(gamma), s=sin(gamma), and T_j=X_j X_(j+1)+Y_j Y_(j+1). Then

\[
e_j=\frac12[-T_j+c(1-Z_jZ_{j+1})-is(Z_j-Z_{j+1})].
\]

Expanding the anticommutator and summing periodic translates cancels the one-site imaginary terms and the nearest-neighbor XX+YY terms. What remains is

\[
H_L=\sum_j\left[
\frac{\cos2\gamma}{2}+s^2Z_jZ_{j+1}
-\frac12(X_jX_{j+2}+Y_jY_{j+2}+Z_jZ_{j+2})
+\frac{is}{2}(T_jZ_{j+2}-Z_jT_{j+1})\right].
\]

Every operator inside the last parentheses is Hermitian, since the Z acts off the XX or YY bond. Its coefficient is imaginary. Summing does not cancel it: the imaginary hopping on bond (j,j+1) is

\[
\frac{is}{2}T_j(Z_{j+2}-Z_{j-1}).
\]

For a configuration with z_(j+2)=+1, z_(j-1)=-1 and an antiparallel pair on (j,j+1), the matrix element swapping that pair is 2is. The reverse matrix element is also 2is, not its conjugate. For L>=5 no next-nearest-neighbor real hopping connects the same two configurations. This proves non-Hermiticity for 0<gamma<pi/2.

Blocking sites (2j,2j+1) gives d_block=4 and the fixed NN density

\[
h_{\rm block}=g\otimes I_2+I_2\otimes g,
\]

where the four microscopic sites are ordered along the chain. Thus g acts on the first three sites in the first term and on the last three in the second. Blocking regroups the basis and cannot change H-H^dagger or its spectrum.

At gamma=0.4 and L=6, ED finds E=2.053095607496681+0.847585247891031i with residual 5.1e-15. Thus the three-block periodic d=4 chain cannot be made Hermitian by any positive metric. This excludes a local product metric as well as an arbitrary finite-size positive metric at this parameter. It does not classify every other staggering or boundary choice.

As a normalization check P, at gamma=0 the formula reduces to H=sum_j(1-P_(j,j+2)), two ferromagnetic permutation chains. The opposite sign is two antiferromagnetic spin-1/2 chains and has c=c_eff=2 (G), with blocked d=4. This is a useful compact control, not a black-hole counterexample.

Verdict: d=2 with three-site terms, or d=4 with NN terms after pair blocking. The explicit Hamiltonian is non-Hermitian. The conditional black-hole values are c=2+6/(k-2), c_eff=2. At 2<k<8 the Virasoro comparison would violate the blocked bound, but this implementation fails a hypothesis. Hermitianizing by (H+H^dagger)/2 deletes a bulk interaction and produces a different model; its continuum theory cannot be inherited from H.

4. What the finite-size calculations establish

The code uses charge sectors and sparse non-Hermitian eigensolvers, not eigsh on a non-Hermitian matrix. Translation eigenvalues are resolved within returned degenerate eigenspaces. Incomplete multiplets are explicitly flagged by translation residuals. Separate momentum-projected calculations remove this ambiguity for the reported p=0 and p=2pi/L levels.

The largest spin-1 charge-zero sector has 73,789 states at L=12. The largest spin-1/2 charge-zero sector has 184,756 states at L=20, corresponding to ten blocked sites. All nonnegative magnetization sectors were checked at these largest sizes. Spin reversal maps the negative-charge spectra to the conjugate positive-charge spectra for these densities, so the smallest real parts are covered. Within numerical accuracy the lowest-real-part state is real and is in the zero-charge sector. Calling it a ground state is a spectral convention for these non-Hermitian examples, not a variational principle.

The independent checks also reproduce the analytic IK anti-Hermitian entry to 1.3e-15, the periodic TL Pauli expansion to 4.3e-15, and the blocked Hamiltonian with zero numerical difference. The R-matrix finite-difference derivative matches the IK density to 2.8e-10. The graded-impurity checks have zero residual.

The momentum-projected eigenvector residuals used below are at most 2.3e-11. The broad charge-sector runs include a few residuals around 1.7e-9. Results are reported at much lower precision than either tolerance.

Spin-1 candidate, theta=0.19. Delta_0 is the first nonzero neutral p=0 gap. The velocity proxy is L[E_min(p=2pi/L)-E_0]/(2pi).

| L | E_0 | L Delta_0 | velocity proxy |
| --- | ---: | ---: | ---: |
| 4 | -1.541015008 | 1.153434 | 2.264122 |
| 6 | -1.978497565 | 0.890155 | 1.973996 |
| 8 | -2.483091580 | 0.776859 | 1.816249 |
| 10 | -3.014309909 | 0.711110 | 1.729463 |
| 12 | -3.558843449 | 0.666678 | 1.678031 |

At L=12 the first four nonzero neutral p=0 gaps, divided by the first, are 1, 3.9674, 9.6527, 18.4694. The first three are roughly consistent with a quadratic radial sequence. This is not an exact n^2 tower, and the fourth differs substantially from 16. Complex low-lying momentum-one levels also occur; they cannot be discarded when assessing microscopic unitarity.

Staggered/TL candidate, gamma=0.4. Sizes divisible by four are used in this table because the other even-size subsequence has a different ground-state degeneracy.

| microscopic L | E_0 | L Delta_0 | velocity proxy |
| --- | ---: | ---: | ---: |
| 4 | -1.213173163 | 2.426346 | 2.932644 |
| 8 | -1.645701810 | 0.850955 | 2.988505 |
| 12 | -2.267394249 | 0.583152 | 2.558486 |
| 16 | -2.927003630 | 0.477359 | 2.317459 |
| 20 | -3.602519996 | 0.419757 | 2.181484 |

The low neutral levels appear at both p=0 and p=pi. They coincide in momentum after blocking pairs. At L=20 the first p=0 gap ratios are 1, 4.7921, 12.9203, 26.6481. These are not enough to identify exact black-hole conformal towers or a scattering phase.

A logarithmic radial quantization would give Delta_0 approximately A/[L(log L+B)^2]. As a short-size diagnostic, fitting 1/sqrt(L Delta_0) linearly against log L gives the following N results, excluding the smallest size. For IK, slope 0.2381 and intercept 0.6358, with maximum ordinate residual 0.0036. For TL, slope 0.5031 and intercept 0.0464, with residual 0.0129. A power-law crossover can mimic this over the accessible interval. These fits do not prove a logarithmic asymptotic or a continuous spectrum.

For ground energies, a two-size estimate of v c_eff is

\[
(vc_{\rm eff})_{L_1,L_2}
=-\frac6\pi\frac{E_0(L_2)/L_2-E_0(L_1)/L_1}{L_2^{-2}-L_1^{-2}}.
\]

It gives 3.03816 for IK using L=10,12 and 3.81867 for TL using L=16,20. If the first momentum-one level is a dimension-one descendant, dividing by the velocity proxies gives apparent charges 1.81055 and 1.75049, respectively. That descendant assignment is A, and finite-size velocity corrections are large.

Two recalled or guessed velocity formulas, neither derived nor independently verified here, are v_IK=pi/(pi-6theta) and v_TL=2pi sin(2gamma)/(pi-2gamma). They give 1.56955 and 1.92488 at the parameters used. Inserting them instead gives apparent c_eff=1.93569 and 1.98385. These numbers are compatible with 2, but relying on these velocities to claim a confirmation would be circular without the missing derivation. All these diagnostics are saved with their assumption labels.

What was not established: the k dictionary, an asymptotic central charge, the full conformal towers, the density-of-states scattering phase, continuum unitarity, or equality of the particular continuum prescriptions in the named papers. No finite-size energy fit here measures the larger Virasoro c directly.

5. Compact integrable families

An affine-algebra label alone is not a complete Hamiltonian specification. One must also specify the representation, deformation parameter, real form, sign, boundary prescription, and regime. The following table separates useful compact benchmarks from facts proved about a specified density. Every compact unitary entry has c_eff=c if its stated CFT identification holds. All ranges in this table are NN with no blocking.

| Candidate / branch | d | Hermiticity in the ordinary inner product | CFT unitarity | c and c_eff | Verdict and status |
| --- | ---: | --- | --- | --- | --- |
| SU(d) fundamental ULS, h=P | d | yes | yes | d-1 | saturation, G |
| B_n^(1), rational SO(2n+1) vector antiferromagnet | 2n+1 | yes for h=P-2K/(2n-1) | yes if SO(2n+1)_1 | n+1/2 | below 2n; WZW identification A |
| D_n^(1), rational SO(2n) vector antiferromagnet | 2n | yes for h=P-2K/(2n-2) | yes if SO(2n)_1 | n | below 2n-1; WZW identification A |
| C_n^(1), rational Sp(2n) vector candidate h=P+K_sp/(n+1) | 2n | yes | compact identification not established here | if pure Sp(2n)_1, n(2n+1)/(n+2); actual limit not determined | conditional benchmark below bound, not a full resolution |
| B_n^(1) trigonometric branches | 2n+1 | depends on parameter and gauge; not checked generally | not classified | n for a pure n-boson branch; n+1/2 if an additional Majorana survives | these proposed compact values obey the bound; branch identifications A |
| C_n^(1) trigonometric branches | 2n | not classified | not classified | n for a proposed n-boson branch; other regimes undetermined | conditional only |
| D_n^(1) trigonometric branches | 2n | not classified | not classified | n for a pure n-boson branch | conditional only |
| D_n^(2), fundamental vector vertex model | 2n | not checked generally | not determined | undetermined | unresolved closed book; no claimed exclusion |
| A_(2n)^(2), compact regimes | 2n+1 | not checked generally | not determined | recalled candidate values n and n+1/2, not verified | both values obey 2n, but the family is not classified |
| A_2^(2) explicit non-compact candidate | 3 | no, proved above | hypothetical coset yes; lattice identification A | conditional c=2+6/(k-2), c_eff=2 | out of scope |
| G_2^(1), seven-dimensional fundamental | 7 | a real combination of orthogonal compact G_2 projectors is Hermitian; the integrable coefficients were not reconstructed | yes if G_2,1 | if pure G_2,1, 14/5 | below 6, conditional identification only |
| SU(N), symmetric p-box fused antiferromagnet | binomial(N+p-1,p) | Hermitian rational compact projector construction, A for the standard integrable realization | yes if SU(N)_p | p(N^2-1)/(p+N) | inequality proved below, CFT identification A |
| SU(N), antisymmetric p-box fundamental fused representation, 1<=p<N | binomial(N,p) | Hermitian rational compact projector construction, A for the standard realization | yes if SU(N)_1 | N-1 | inequality proved below, CFT identification A |
| SU(2), spin S Takhtajan–Babujian branch | 2S+1 | yes for the real polynomial density | yes if SU(2)_(2S) | 3S/(S+1) | below 2S, equality only S=1/2; identification A |
| qutrit XX and finite-density t–J examples in the question | 3 | yes, G | yes, G | 2 | saturation, G |
| two next-neighbor spin-1/2 Heisenberg chains | 4 after blocking | yes | yes | 2 | below 3, G |
| three-site spin-1/2 c=3/2 example in the question | 4 after blocking | yes, G | yes, G | 3/2 | below 3 |

Here K=|Omega><Omega| with |Omega>=sum_a|a,a>. Thus the SO(M) density in the table is manifestly Hermitian. K_sp is the positive rank-one operator |J><J| associated with the real antisymmetric symplectic form, of norm squared 2n. Its stated density is also manifestly Hermitian. The osp construction below supplies the same rational densities up to the appropriate grading and sign. These algebraic checks do not establish their continuum limits.

I have deliberately not equated Sp(2n) symmetry with a proved Sp(2n)_1 continuum theory. Competing recollections about several velocities or additional massless sectors are not a calculation. The same caution applies to G_2 and the twisted families. The unresolved rows prevent a claim that every family requested in the question has been disposed of.

Several useful inequalities are nevertheless P. For a compact Sp(2n)_1 benchmark,

\[
(2n-1)-\frac{n(2n+1)}{n+2}=\frac{2(n-1)}{n+2}\ge0.
\]

For SO(M)_1, c=M/2<=M-1 for M>=2. The familiar vector embedding ceiling SO(M)_2 has c=M-1, consistent with the supplied embedding argument. For G_2,1, c=14/(1+4)=14/5<6. These checks rule out the displayed fixed points, not arbitrary levels or additional commuting sectors.

For the SU(N) symmetric-p benchmark,

\[
c=\frac{p(N^2-1)}{p+N}
=p(N-1)\frac{N+1}{N+p}\le p(N-1).
\]

There are at least 1+p(N-1) degree-p monomials in N variables: x_1^p and x_1^(p-r)x_j^r for 2<=j<=N and 1<=r<=p. Hence d>=1+p(N-1) and c<=d-1. For antisymmetric p, binomial(N,p)>=N for 1<=p<N, so the assumed SU(N)_1 charge N-1 obeys the bound. This includes conjugate fundamentals. The singlet endpoints p=0,N are not critical chains.

For spin S, direct subtraction gives

\[
2S-\frac{3S}{S+1}=\frac{S(2S-1)}{S+1}\ge0.
\]

For S>=3/2, any pure finite-level unitary SU(2) WZW model has c<3<=2S, independently of the level, as in the supplied background. None of these arguments excludes an additional gapless sector not included in the assumed WZW theory.

Other compact benchmarks one might test are F_4,1 in dimension 26 with c=26/5, E_6,1 in dimension 27 with c=6, and E_7,1 in dimension 56 with c=7. The arithmetic follows from the WZW formula. No nearest-neighbor lattice-to-CFT identification for these extra examples was established here, so they are not additional verified cases.

6. Graded permutations and orthosymplectic chains

There is an exact classification relevant to the unrestricted bare graded-permutation Hamiltonian. Assign m even colors and n odd colors, d=m+n, and define on the ordinary positive-definite tensor product

\[
P_g|a,b\rangle=(-1)^{p(a)p(b)}|b,a\rangle,\qquad H_+=\sum_j P_{g,j,j+1}.
\]

P_g is a real symmetric involution. Therefore H_+ is Hermitian and H_+>=-L. If n>=1, the state with the same odd color f on every site attains -L. This is an exact ground state.

If m>=1, put one even color b into this background and call its position state |j>. Only the two adjacent bonds change. Their diagonal energy changes from -2 to 0, and they hop the impurity with amplitude +1. Thus

\[
(H_++L)|j\rangle=2|j\rangle+|j-1\rangle+|j+1\rangle,
\qquad\epsilon(p)=2+2\cos p.
\]

At p=pi+q, epsilon=q^2+O(q^4). If instead n>=2, use a second odd color as the impurity. The hopping is -1 and epsilon(p)=2-2cos p=q^2+O(q^4) near zero. These states give a gap at most O(L^-2), or exact extra ground states and O(L^-2) excitations next to them. They exclude an unrestricted low-energy description consisting only of ordinary relativistic CFTs with finite nonzero velocities.

This proof is for the ordinary tensor-product periodic convention entering the conjecture. A graded supertrace boundary convention can shift allowed momenta or introduce a parity-sector twist. Such a shift does not turn the quadratic minimum into a relativistic branch, but it must be specified if comparing detailed Bethe spectra.

| Bare graded-permutation case | d | range / Hermiticity | c, c_eff and unitarity verdict |
| --- | ---: | --- | --- |
| H_+, n=0, m>=2 | m | NN, Hermitian | ordinary SU(m) ULS, c=c_eff=m-1, unitary, G |
| H_+, n>=1, d>=2 | m+n | NN, Hermitian | quadratic low-energy branch proved; a purely conformal c is not applicable |
| H_+, d=1 | 1 | constant Hamiltonian | trivial, no positive conformal charge |
| H_-=-sum P_g, m=0, n>=2 | n | NN, Hermitian | ordinary SU(n) antiferromagnetic ULS, c=c_eff=n-1 |
| H_-, m>=1, d>=2 | m+n | NN, Hermitian | analogous polarized even-color ground state and quadratic branch |

Fixing a nonzero density or adding chemical potentials is a different ground-state problem. A fixed-density constrained Hilbert space is not literally the unrestricted tensor-product ground-state problem in the conjecture. Chemical potentials can be encoded in a NN h, so those models can be in scope and require a separate analysis. In particular, the finite-density t–J example supplied in the question is not contradicted by the polarized-ground-state result. I did not classify all chemical potentials for all m,n, and I have not assigned the formal superalgebra WZW central charge to a positive-metric chain merely from its name.

For a conventional rational osp(m|2n) density, write d=m+2n, delta=m-2n and define

\[
|b\rangle=\sum_{a=1}^{m}|a,a\rangle,\qquad
|f\rangle=\sum_{r=1}^{n}(|u_r,v_r\rangle-|v_r,u_r\rangle),
\]

\[
E=|b+f\rangle\langle b-f|,\qquad
E^2=\delta E,\qquad
h_{\rm osp}=P_g-\frac{2}{\delta-2}E.
\]

This specifies a matrix, and its Hermiticity analysis is P. Its identification with a particular graded transfer-matrix convention is A. The displayed normalization is singular at delta=2. The formal rescaled density (delta-2)P_g-2E has a finite value -2E there and is still non-Hermitian when both m and n are nonzero. Whether that rescaling is the appropriate regular integrable chain at delta=2 was not established.

For m,n>0,

\[
E-E^\dagger=2(|f\rangle\langle b|-|b\rangle\langle f|)\ne0.
\]

The anti-Hermitian entries change both colors on a bond, from an even repeated pair to an odd symplectic pair. As in the IK argument, the periodic sum cannot cancel them. Thus the displayed mixed osp density is not an in-scope Hermitian Hamiltonian for delta!=2.

At n=0 one obtains the Hermitian SO(m) density. At m=0, E=-|f><f| and P_g=-P, so reversing the overall sign gives the Hermitian symplectic candidate P+|f><f|/(n+1) used above. The mixed cases were checked numerically for osp(1|2), osp(2|2), and osp(3|2); their periodic L=3 anti-Hermitian ratios are approximately 0.6426, 0.8944, and 1.1429.

Recalled continuum statements such as a Goldstone or supersphere value c=m-2n-1 are A and are not a uniform classification of osp chains, boundary sectors, and effective charges. They need not describe a unitary CFT in a positive inner product; c_eff cannot be read off without h_min. The arithmetic c=m-2n-1<=m+2n-1 would obey the Virasoro bound if that formula applied, but this does not establish either unitarity or a c_eff bound. Other exceptional delta=2 constructions and any alternative explicitly Hermitian realization remain unresolved here.

7. Attempts to construct a violation

First, SU(d-1)_1 plus one Ising sector would have c=(d-2)+1/2=d-3/2, below the target. SU(d-1)_1 plus one boson would have c=d-1, saturation. These suggested additions do not violate the inequality even if they can be embedded. SU(d)_1 plus any nonzero sector would violate it, but it requires a new construction retaining the original SU(d)_1 sector while adding degrees of freedom.

A tensor product of two already-bounded chains does not work. If the factors have dimensions d_1,d_2 and charges c_1,c_2, the on-site dimension is D=d_1d_2 and

\[
c_1+c_2\le d_1+d_2-2\le d_1d_2-1,
\]

since the difference is (d_1-1)(d_2-1). The density h_1 tensor identity plus identity tensor h_2 is NN and Hermitian, so this is a useful construction check rather than a failure of locality. A chain satisfying its own bound cannot be converted into a violation by this tensor-product operation.

A direct sum of local state spaces also does not automatically add central charges. Introduce projectors Q_1,Q_2 onto the two spaces and a sufficiently large positive domain-wall penalty on Q_1 tensor Q_2 + Q_2 tensor Q_1. The low-energy uniform-label sectors are a direct sum of two chain Hilbert spaces. They are not two simultaneous propagating CFTs in each state. Choosing a vacuum gives one CFT; a finite mixture of vacua contributes at most a constant entropy. One cannot sum their central charges in the conjecture's sense. Tuning domain walls to be gapless creates a new interacting model whose continuum must be computed.

An on-site isometric compression of a tensor product cannot inject D independent states into d<D states. Projection or a gauge constraint removes states and may change the spectrum, the local observable algebra, or the interaction range. This is not a theorem against all constrained realizations. Spin-charge separation already shows that distinct continuum sectors need not correspond to independent on-site tensor factors. The supplied t–J saturation is an explicit warning against overusing local state counting. A multicritical t–J-like chain retaining all saturated sectors and adding an Ising sector remains a possible research route, not a construction achieved here.

For the d=2 next-nearest-neighbor Heisenberg example, blocking produces d=4 and c=2<3. On-site changes of basis cannot shorten a term's support. Reordering odd sites followed by even sites displays two separate periodic chains, but does not embed both cycle graphs into a single connected NN ring with one fixed nonzero bond density. A nonlocal duality is a separate possibility.

The standard Ising duality illustrates the obstruction without proving a general no-go. With

\[
X_j=\tau^z_{j-1/2}\tau^z_{j+1/2},\qquad
Z_jZ_{j+1}=\tau^x_{j+1/2},
\]

one has

\[
Z_jZ_{j+2}=\tau^x_{j+1/2}\tau^x_{j+3/2},
\]

which becomes NN. But the Heisenberg XX term becomes

\[
X_jX_{j+2}=\tau^z_{j-1/2}\tau^z_{j+1/2}\tau^z_{j+3/2}\tau^z_{j+5/2},
\]

a four-site interaction. The YY term does not remove this independent operator string. Thus this particular duality does not produce the desired NN qubit chain. Duality also carries global constraints and boundary sectors that must be retained. No general classification of nonlocal or finite-depth dualities was performed, and no no-go theorem for all such transformations is claimed.

The given three-site c=3/2 model is likewise safe after blocking, since d_block=4. I did not construct a range-shortening transformation for it. Non-on-site symmetries can constrain possible fixed points, but they neither provide the missing NN Hamiltonian nor supply an upper bound on total central charge by themselves.

8. Definite conclusions and unresolved work

For the two explicit non-compact densities in this report, the answer is definite: they are not counterexamples to either in-scope conjecture, because they are not Hermitian. Their non-Hermiticity is a bulk effect. The complex-spectrum checks additionally rule out positive-metric quasi-Hermiticity at the tested parameters. Conditional on the black-hole identification they illustrate why c and c_eff must not be conflated, but they do not settle the conjecture for Hermitian chains.

For the bare positive-metric graded permutations, the unrestricted mixed-color ground-state problem has an explicit quadratic branch rather than the hypothesized purely conformal low-energy limit. For the compact fixed points whose identifications were supplied or explicitly assumed, all charge comparisons obey the bound. The mixed rational osp construction displayed here fails Hermiticity.

No in-scope counterexample was found for c. No in-scope counterexample was found for c_eff. I cannot honestly replace either statement by “no counterexample exists.” The outstanding steps are (i) verify the precise lattice and continuum dictionaries of the named non-compact models, including alternative Hermitian realizations if any, (ii) classify the unresolved twisted, symplectic, exceptional, and finite-density superalgebra regimes, and (iii) find or rule out a NN realization of an additional sector beyond a saturated compact chain. Those are substantive gaps, not numerical error bars.

9. Reproduction and files

The standalone d1.py contains the matrix constructions, sector and momentum ED, diagnostics, direct six-vertex transfer check, and symbolic algebra checks. It uses NumPy and SciPy, with SymPy needed only for the symbolic command. Run from this directory with one BLAS thread, for example:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py algebra
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py validate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py symbolic
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py transfer
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py ed
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py momentum
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py audit
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 d1.py analyze
```

The analyze command reads d1_momenta.json generated by the momentum command. The computation does not need input data from the project. JSON files preserve eigenvalues, residuals, magnetizations, momenta, Hermiticity checks, and the assumption-dependent diagnostics. The companion scripts are earlier runnable versions of routines also incorporated into d1.py. No project STATE.md outside this working directory was accessed or updated.
