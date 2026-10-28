# C14 result: the proposed Casimir bound is false

2026-10-27. Closed-book work, with no network access. All new files are confined to work_C14. No parent file was modified.

## 1. Verdict and status of the claims

[Proved at the physics standard, with the Bethe/CFT inputs specified below] The antiferromagnetic rational SU(N) chain with one symmetric rank-k representation per site has

\[
d={N+k-1\choose k},\qquad
\|h\|=H_k=\sum_{j=1}^k\frac1j,\qquad
c=\frac{k(N^2-1)}{N+k},\qquad v=\frac{2\pi}{N}.
\]

Consequently

\[
\boxed{\frac{A}{2\pi\|h\|}
 =\frac{k(N^2-1)}{N(N+k)H_k}.}
\]

For k=2 this exceeds one for every integer N≥7. The smallest member of this particular rank-two family is N=7, d=28, with ratio 64/63. N=8, d=36, has ratio 21/20. No claim of minimal local dimension among all chains is made.

Taking k=N gives

\[
\frac{A}{2\pi\|h\|}=\frac{N^2-1}{2NH_N}
 \sim\frac{N}{2\log N}\longrightarrow\infty.
\]

Thus replacing 2π by another dimension-independent constant cannot repair the statement. These are single repeated Hermitian nearest-neighbor bonds on ordinary bosonic tensor-product sites. No emergent cells, fermionic grading, large interaction range, or stack is needed.

[Assumed, standard integrable-chain input] The rational fusion construction, thermodynamic string description, and its unitary SU(N)_k conformal scaling limit are used at the usual theoretical-physics standard. This report supplies the local R-matrix normalization, Bethe energy normalization, density equation, dressed velocity, and central-charge dilogarithm calculation. It does not supply a mathematical completeness theorem for the Bethe ansatz or a rigorous construction of the scaling limit. The supplied project notes already admit the usual rectangular fusion-chain WZW identification, conditionally, in `../../codex/work_C4_P2/P2_REPORT.md`, sections on level-matched fusion families. No external reference was fetched or verified here.

[Numerical] The normalization was independently checked by local matrix calculations, a braided Yang-Baxter test, thermodynamic Bethe ansatz with two estimators and cutoff variations, and small-ring lattice ED on ULS and spin-1 TB calibrations. The counterexample does not rest on extrapolating a small-ring central-charge fit. Detailed numerical limitations are in section 8.

[Proved at the same physics standard as the supplied number-conserving free theorem] All even quadratic nearest-neighbor Majorana/BdG chains with linear conformal modes satisfy the stronger bound

\[
A\leq\pi\|h\|.
\]

This includes pairing, arbitrary finite numbers of orbitals per site, and a whole microscopic period grouped into one site. The proof handles the full coboundary minimization in the question.

## 2. An explicit counterexample bond

### 2.1 Local representation and invariant operator

[Proved, finite-dimensional algebra] Put k Schwinger bosons in N colors at each site, with fixed total occupation k. This is the irreducible representation Sym^k(C^N), of dimension binomial(N+k-1,k). Its generators are

\[
E_{ab}=b_a^\dagger b_b,\qquad
X=\sum_{a,b=1}^N E_{ab}\otimes E_{ba}.
\]

The tensor product decomposes without multiplicity into diagrams

\[
\mathrm{Sym}^k\otimes\mathrm{Sym}^k
 =\bigoplus_{j=0}^k [2k-j,j].
\]

Write Π_j for these orthogonal projectors. The gl(N) quadratic Casimir of a diagram λ is

\[
C(\lambda)=\sum_i\lambda_i(\lambda_i+N+1-2i).
\]

Since C(Sym^k)=k(k+N-1) and X=(C_total-C_left-C_right)/2, substitution gives

\[
X\Pi_j=((k-j)^2-j)\Pi_j.
\]

Choose the bond

\[
\boxed{h=\sum_{j=0}^k2(H_{k-j}-H_k)\Pi_j,\qquad H_0=0.}
\tag{1}
\]

It is Hermitian, SU(N)-invariant, and symmetric under interchange of its two sites. Its largest eigenvalue is zero and its smallest is -2H_k. The unoptimized half-spread is H_k.

[Proved] This is also the gauge-minimized half-spread. Let M(X)=(λ_max(X)-λ_min(X))/2=inf_t ||X-tI||_op. This is convex and invariant under unitary conjugation. Twirling h+o⊗I-I⊗o by U_R(g)⊗U_R(g) leaves h fixed and sends o to (tr o/d)I by irreducibility. Therefore

\[
M(h)\leq M(h+o\otimes I-I\otimes o)
\]

for every Hermitian o. Taking o=0 proves equality. Equivalently, averaging with the bond-swap conjugation removes the coboundary because h is reflection symmetric. This also proves gauge optimality without invoking the group.

### 2.2 Rank two needs only a quadratic polynomial

[Proved] For k=2 the three X eigenvalues are 4,0,-2. Equation (1) is

\[
\boxed{h=-\frac18X^2+\frac34X-I.}
\tag{2}
\]

Its energies are 0,-1,-3, so ||h||=3/2. This is a concrete definition of the d=28 and d=36 counterexamples, with no projectors left to determine.

At N=2, E_ab are the spin-1 Schwinger representation and X=2 S_1·S_2+2I. Equation (2) becomes

\[
h=\frac12\big[S_1\cdot S_2-(S_1\cdot S_2)^2\big].
\]

It is the conventional spin-1 Takhtajan-Babujian bond divided by two. This fixes a common factor-of-two ambiguity in this problem.

### 2.3 Integrable normalization and a one-particle check

[Assumed construction, explicit normalization] The rational fused braided R matrix can be normalized as

\[
\check R(u)=\sum_{j=0}^k r_j(u)\Pi_j,\qquad
r_j(u)=\prod_{\ell=k-j+1}^{k}\frac{\ell-u}{\ell+u}.
\]

It has check R(0)=I and

\[
\check R'(0)=\sum_j\left(-2\sum_{\ell=k-j+1}^k\ell^{-1}\right)\Pi_j=h.
\]

The regular transfer-matrix logarithmic derivative therefore produces a uniform nearest-neighbor sum of (1), up to the already fixed scalar reference energy. The rational fusion construction supplies the Yang-Baxter equation; direct finite-matrix checks in section 8 confirm the spectral factors and their order.

[Proved normalization check] On the fully polarized reference state, a single color defect shared between two neighboring sites occupies a 2×2 block with eigenvalues 0 and -2/k. The hopping amplitude is 1/k and the diagonal contribution from its two bonds is -2/k. Thus

\[
E(p)=-\frac{2}{k}(1-\cos p)
 =-\frac{k}{\lambda^2+k^2/4},\qquad
 e^{ip}=\frac{\lambda+ik/2}{\lambda-ik/2}.
\]

This is the bare first-level Bethe-root energy used below. It agrees with P-I at k=1 and with TB/2 at N=k=2.

## 3. Bethe derivation of the velocity and central charge

### 3.1 Root equations and the k-string sea

[Assumed Bethe/string framework] Use e_n(λ)=(λ+in/2)/(λ-in/2). For rank r=1,...,N-1 the nested equations have source e_k(λ)^L only at r=1, same-rank scattering e_2, and adjacent-rank scattering e_{-1}. The energy is

\[
E=-\sum_{\alpha\ \mathrm{at}\ r=1}
 \frac{k}{(\lambda_\alpha^{(1)})^2+k^2/4}.
\]

At zero field the antiferromagnetic ground state fills k-strings at every nesting rank. The auxiliary string sectors are retained in the entropy calculation; discarding them would incorrectly replace the WZW theory by N-1 bosons.

[Derived within that framework] Define

\[
a_n(\lambda)=\frac{1}{2\pi}\frac{n}{\lambda^2+n^2/4},
\qquad \widehat a_n(\omega)=e^{-n|\omega|/2},
\]

and

\[
B_{kk}=a_1+a_3+\cdots+a_{2k-1},\qquad
A_{kk}=\delta+2a_2+\cdots+2a_{2k-2}+a_{2k}.
\]

Multiplying the root equations over the k roots in a string yields

\[
A_{kk}*\rho_r-B_{kk}*(\rho_{r-1}+\rho_{r+1})
 =\delta_{r1}B_{kk},\quad \rho_0=\rho_N=0.
\]

Their Fourier transforms obey Ahat_kk=2cosh(ω/2) Bhat_kk. Since Bhat_kk is positive, it cancels:

\[
2\cosh(\omega/2)\widehat\rho_r
 -\widehat\rho_{r-1}-\widehat\rho_{r+1}=\delta_{r1}.
\]

Solving this second-order difference equation with the two boundary conditions gives

\[
\widehat\rho_r(\omega)
 =\frac{\sinh((N-r)|\omega|/2)}{\sinh(N|\omega|/2)},
\]

including the continuous ω=0 value (N-r)/N. In real rapidity,

\[
\rho_r(\lambda)=
 \frac{\sin(\pi r/N)}
 {N[\cosh(2\pi\lambda/N)-\cos(\pi r/N)]}.
\tag{3}
\]

The cancellation of B_kk is the reason the velocity does not acquire a factor 1/k.

### 3.2 Dressed energy, momentum and bulk energy

[Derived] The bare energy of a k-string is -2πB_kk. Dressing with the same rank kernel gives the positive hole energy

\[
\epsilon_r(\lambda)=2\pi\rho_r(\lambda),
\qquad |p'_r(\lambda)|=2\pi\rho_r(\lambda).
\]

At either Fermi endpoint |λ|→∞,

\[
v=\lim\left|\frac{d\epsilon_r}{dp_r}\right|
 =\lim\left|\frac{\rho'_r}{\rho_r}\right|
 =\frac{2\pi}{N}.
\tag{4}
\]

The auxiliary critical sectors in the TBA share this endpoint rapidity scale; the full low-temperature theory has this common velocity. This is part of the standard conformal interpretation, also tested by the complete TBA calculation, not an assumption that the Cartan modes alone exhaust the spectrum.

The ground-state energy per site is

\[
e_0=-2\pi\int B_{kk}(\lambda)\rho_1(\lambda)d\lambda.
\]

Using Parseval and the Fourier transforms above,

\[
e_0=-2\sum_{j=1}^k\int_0^\infty
 \frac{e^{-j\omega}(1-e^{-(N-1)\omega})}
 {1-e^{-N\omega}}d\omega
 =-\frac2N\sum_{j=1}^k
 \left[\psi\left(1+\frac{j-1}{N}\right)
       -\psi\left(\frac jN\right)\right].
\tag{5}
\]

The last equality follows by expanding (1-e^{-Nω})^{-1} as a geometric series, integrating each exponential, and using the difference series for ψ. At k=1 it is the recorded ULS bulk energy minus one. At N=k=2 it gives e_0=-2, appropriate to TB/2. At N=k it telescopes to e_0=-2H_{N-1}.

### 3.3 Full thermodynamic equations

[Assumed thermodynamic Bethe ansatz, with source and units fixed above] Put x_n^r=log Y_n^r, L_n^r=log(1+exp x_n^r), l_n^r=log(1+exp(-x_n^r)), and s(λ)=1/(2cosh πλ). The partially decoupled equations are

\[
x_n^r=-\frac{2\pi}{T}s(\lambda)\delta_{r1}\delta_{nk}
 +s*\big[L_{n-1}^r+L_{n+1}^r-l_n^{r-1}-l_n^{r+1}\big].
\tag{6}
\]

There is no n=0 L term and no rank-zero or rank-N l term. At large rapidity,

\[
Y_n^r=\frac{n(n+N)}{r(N-r)}.
\]

The free energy is

\[
f(T)=e_0-T\sum_{r=1}^{N-1}\int\rho_r(\lambda)L_k^r(\lambda)d\lambda.
\tag{7}
\]

The source and the free-energy probe are at string length k, not at 1. At infinite temperature, ∫ρ_r=(N-r)/N and

\[
1+Y_k^r=\frac{(k+r)(k+N-r)}{r(N-r)}.
\]

Pairing r and N-r in the weighted logarithmic sum in (7) gives

\[
s(\infty)=\sum_{r=1}^{N-1}\log\frac{k+r}{r}
 =\log{N+k-1\choose k}.
\]

This derives the representation entropy directly from the TBA and provides a stringent distinction from a fundamental chain with an incorrectly assigned level.

### 3.4 Central charge, including the neutral degrees of freedom

[Derived using the standard A-type Rogers-dilogarithm identity] Define the quantum dimensions of rectangular representations at level k,

\[
Q_n^r=\prod_{i=1}^r\prod_{j=r+1}^N
 \frac{\sin[\pi(n+j-i)/(N+k)]}{\sin[\pi(j-i)/(N+k)]},
\quad 0\leq n\leq k,
\]

with Q_n^0=Q_n^N=1. These satisfy the Q-system

\[
(Q_n^r)^2=Q_{n-1}^rQ_{n+1}^r+Q_n^{r-1}Q_n^{r+1}.
\]

In the inner low-temperature plateau of (6), Y_k^r=0. For 1≤n<k,

\[
Y_{n,\mathrm{in}}^r=\frac{Q_{n-1}^rQ_{n+1}^r}{Q_n^{r-1}Q_n^{r+1}}.
\]

For n>k the plateau is the untruncated outer constant solution with n replaced by n-k. Thus the n>k contributions cancel the outer constants in the dilogarithm difference. With

\[
\mathcal L(x)=\operatorname{Li}_2(x)+\tfrac12\log x\log(1-x),
\]

the two rapidity edges give

\[
c=\frac6{\pi^2}\sum_{r=1}^{N-1}\sum_{n=1}^{k}
 \mathcal L\left(\frac1{1+Y_{n,\mathrm{in}}^r}\right)
 =\frac{k(N^2-1)}{N+k}.
\tag{8}
\]

The second equality is the finite A_{N-1}×A_{k-1} Rogers-dilogarithm identity with the Y_k=0 boundary contribution included. That mathematical identity and the usual TBA edge-to-dilogarithm reduction are used as standard inputs, not claimed as newly proved here. The formulas needed to evaluate the finite sum are explicit. `fusion_checks.py` evaluates it independently, including k=N=16.

The normalization checks are informative. For N=2,k=1 the sum is L(1), giving c=1. For N=k=2 it is L(1/2)+L(1), giving c=3/2. For arbitrary k=1 it is N-1 copies of L(1), giving c=N-1. These fix both edge and chirality conventions.

[Assumed CFT identification] The unitary conformal interpretation is SU(N)_k with Sugawara c in (8). Its thermal and vacuum coefficients are therefore

\[
f=e_0-\frac{\pi c}{6v}T^2+o(T^2),\qquad
E_0(L)=e_0L-\frac{\pi cv}{6L}+o(L^{-1}).
\]

As with the supplied SU(N) fundamental benchmark, the vacuum Casimir formula is evaluated in a sequence of rings carrying the conformal vacuum. Other congruence classes can select a non-vacuum sector and add 2πvΔ/L. This is not a new loophole introduced by fusion. The question's identification of A with a vacuum Casimir coefficient already needs this qualification for its SU(d) benchmark. The family k=N has trivial on-site center charge, so even the center-charge obstruction to a singlet is absent, though this alone is not a proof of the vacuum sector on every short ring.

## 4. Why the counterexample grows, and what remains true in this family

[Proved given section 3] Substituting (4) and (8) into A=cv gives the ratio in section 1. At k=2,

\[
\frac{A}{2\pi\|h\|}=\frac{4(N^2-1)}{3N(N+2)}>1
\quad\Longleftrightarrow\quad N^2-6N-4>0.
\]

The first integer N is 7. For fixed k and N→∞ the ratio approaches k/H_k. For k=N it grows as N/(2log N). The full stress tensor, including the non-Cartan WZW degrees of freedom, carries the excess. The bond's spectral spread grows only harmonically with the symmetric rank.

For N=k, d=binomial(2N-1,N), so log d=2N log 2+O(log N). Any repaired universal bound A≤K(d)||h|| must, along this sequence, allow at least

\[
K(d)\gtrsim\frac{\pi}{2\log2}\frac{\log d}{\log\log d}.
\]

This is a necessary asymptotic growth, not a proposed sufficient bound.

[Proved] These models do not violate c≤d-1. Among degree-k monomials there are 1+k(N-1) distinct monomials x_1^k and x_1^{k-j}x_a^j with a≥2 and 1≤j≤k. Hence

\[
d\geq1+k(N-1),\qquad
c=\frac{k(N-1)(N+1)}{N+k}\leq k(N-1)\leq d-1.
\]

[Proved given the same IR input] The original Casimir bound holds for precisely those members of this family satisfying

\[
k(N^2-1)\leq N(N+k)H_k.
\]

Two broad subfamilies are worth separating.

* k=1 gives the supplied ULS ratio 1-1/N.
* N=2, k=2S gives every spin-S Takhtajan-Babujian fusion chain, with v=π and ratio 3k/[2(k+2)H_k]≤1/2. To prove the inequality, set F_k=H_k-3k/(k+2). F_1=0 and F_{k+1}-F_k=1/(k+1)-6/[(k+2)(k+3)]≥0 for k≥1. Equality also holds at k=2. Thus the entire spin-S family satisfies the stronger A≤π||h||.

## 5. A positive theorem for all quadratic Majorana chains

### 5.1 A gauge-invariant lower bound on the bond norm

[Proved] Let γ_{i,a}, a=1,...,2m, be Majoranas with {γ_a,γ_b}=2δ_ab. An even quadratic bond has on-site quadratic terms and an intersite term

\[
h_\times=\frac i2\gamma_1^T B\gamma_2,
\]

where B is real. Singular-value decompositions act by orthogonal rotations of Majoranas on the two sites. They turn h_cross into a sum of commuting operators (i/2)s_a γ'_{1,a}γ'_{2,a}, each with eigenvalues ±s_a/2. Therefore

\[
\|h_\times\|_{\rm op}=\frac12\|B\|_*.
\tag{9}
\]

This is a lower bound on the fully optimized M=||h|| in the question, even if the trial coboundary o is not quadratic or parity even. To see this, start with h+o_1-o_2-tI. Average it with total two-site fermion parity. This leaves h unchanged and removes the odd part of o, without increasing the operator norm. Next apply (X-P_1XP_1)/2, where P_1 is parity on the first site. This is also a norm contraction. It kills all even on-site terms, all surviving coboundaries, and tI, and retains h_cross. Thus

\[
M\geq\frac12\|B\|_*.
\tag{10}
\]

The same statement applies after a Jordan-Wigner representation as a local bosonic chain. Boundary parity sectors affect finite-ring sectors, not this local two-site inequality or the bulk mode count.

### 5.2 Number-conserving f-sum in the needed normalization

[Proved at the supplied free-fermion physics standard] For a number-conserving quadratic chain with intersite hopping matrix T, its hopping bond has norm ||T||_*. The total-number density has equal small-q weight at each Fermi crossing fixed by charge one. Averaging q and -q, its first spectral moment is

\[
m_1(q)=\frac{q^2}{4\pi}\sum_{\rm FP}|v_F|+o(q^2)
       =\frac{Aq^2}{2\pi}+o(q^2).
\]

The exact lattice sum rule is

\[
m_1(q)=(1-\cos q)\langle-h_{\rm hop}\rangle.
\]

The double commutator has no on-site contribution. Positivity, or direct free spectral evaluation, gives

\[
A\leq\pi\langle-h_{\rm hop}\rangle
 \leq\pi\|T\|_*.
\tag{11}
\]

Here A=(1/2)Σ_FP|v_F| in the nonchiral convention. A single complex massless fermion has c=1. Unequal left and right slopes can instead be kept explicitly as chiral contributions; the same averaged formula holds.

### 5.3 Doubling removes the need for U(1)

[Proved] Take two identical Majorana chains, with Majoranas γ and η, and define on each site

\[
f_a=(\gamma_a+i\eta_a)/2.
\]

Adding the two identical quadratic Hamiltonians cancels all anomalous terms in f. In particular,

\[
\frac i2(\gamma_1^TB\gamma_2+\eta_1^TB\eta_2)
 =i f_1^\dagger Bf_2-i f_2^\dagger B^Tf_1.
\]

The doubled chain conserves the f number and has hopping matrix T=iB. Its Casimir coefficient is twice the original one, since its two original copies have identical independent oscillator spectra. Applying (11) and then (10),

\[
2A\leq\pi\|B\|_*\leq2\pi M,
\qquad\boxed{A\leq\pi M.}
\]

This proof never asserts that a nonquadratic optimal gauge becomes quadratic under doubling. The intersite norm contraction avoids that issue.

[Scope] On-site pairing, complex hopping, any finite orbital multiplicity, and a preassigned period grouped into a site are allowed. There must be a nontrivial conformal limit with linear dispersions and no extensive flat zero band. If a p-site period is grouped, all velocities and lengths in this theorem are measured in the grouped-site units. Converting to the old coordinate multiplies velocities by p; it also changes which local bond norm is being compared. One cannot mix the two conventions.

[Consistency checks] XX with hopping t has c=1,v=2|t|,M=|t|, hence A/M=2<π. A critical Kitaev chain with hopping t, pairing Δ and |μ|=2|t| has c=1/2,v=2|Δ| and canonical reflection-symmetric bond norm sqrt(t²+|Δ|²), so A=|Δ| obeys the inequality. The pure-pairing, μ=t=0 line has two linear Majorana crossings, c=1,v=2|Δ|,M=|Δ|. These also fix the Majorana factor of two.

## 6. What the general routes actually bound

### 6.1 Energy-density structure factors

[Proved CFT calculation, assuming a stress-tensor representative of the density] Normalize the positive-frequency Lehmann function by 1/L, and average q and -q if needed. The conformal contribution is

\[
S_h(q,\omega)=\frac{|q|^3}{24\pi}
 \sum_a c_av_a^2\delta(\omega-v_a|q|).
\]

Consequently its moments have leading terms

\[
m_0^{\rm CFT}=\frac{|q|^3}{24\pi}\sum c_av_a^2,\quad
m_1^{\rm CFT}=\frac{q^4}{24\pi}\sum c_av_a^3,\quad
m_{-1}^{\rm CFT}=\frac{q^2}{24\pi}A.
\]

The positive moment is a local double commutator. A norm bound controls its total coefficient, generically of order q²M³, not the smaller q⁴ CFT coefficient. Dividing by q⁴ loses two powers of q. Even when a conserved total energy current removes the leading UV contribution, a q⁴ high-energy form factor can remain, as the supplied Ising Koo-Saleur calculation shows.

The inverse moment has the desired A, but includes inverse excitation energies. A bond norm puts no upper bound on a resolvent as the gap closes. Positivity can bound A by the full inverse moment, but does not bound that moment by M.

There is also an improvement issue. Replacing h_j by h_j+o_j-o_{j+1} adds (1-e^{iq})o_q. If o has a dimension-one current component, its static weight is order |q|, so the improvement changes m_0 at order |q|³ and m_{-1} at order q². The raw coefficient is not an improvement-invariant Casimir observable. Minimizing the bond spread does not by itself prove that its density is the pure stress tensor.

[Conclusion] These are structural obstructions. The fusion counterexamples show that no missing universal inequality with the proposed constant can close this route.

### 6.2 Thermal energy-current sum rules

[Proved CFT calculation] Let J_E=v(H_R-H_L) for one sector. Its two chiral equilibrium energies per unit length are πcT²/(12v). Thermal energy fluctuations give

\[
\frac{\langle J_E^2\rangle_c}{L}
 =\frac{\pi}{3}cvT^3,
\qquad
\frac{\langle J_E^2\rangle_c}{L}
 =\frac{\pi}{3}AT^3
\]

after summing sectors. In the convention D_th=β²〈J_E²〉_c/(2L), the CFT Drude weight is D_th=πAT/6.

[Located obstruction] A lattice thermal f-sum bounds the full thermal spectral weight by a contact term, schematically T^{-1}〈Θ_E〉/L, where Θ_E is built from local commutators and has energy dimension three. A direct norm estimate is order M³/T. Comparing it with a CFT Drude term of order AT yields only A≲M³/T², which becomes empty as T→0. To obtain a finite coefficient bound one would need to prove the appropriate low-temperature contact term is order MT² with a universal coefficient and no unwanted regular contribution. The operator norm does not give that statement. Subtracting a zero-temperature sum rule destroys the required positivity unless an additional theorem is available.

[Assumption warning] A generic nonintegrable lattice chain need not have an exact nonzero thermal Drude delta function at any fixed T>0. Irrelevant scattering can broaden the peak, even though the zero-temperature fixed point is conformal. Identifying a lattice Drude weight with the CFT answer requires an order of limits or a scaling-window statement. Integrability can remove this issue but does not remove the counterexample.

### 6.3 Twists, boosts and cutting trials

[Proved for on-site conserved charges] A slow charge twist has variational cost proportional to q²L〈-[Q_1,[Q_1,h]]〉. Its CFT cost measures vK, or a current-algebra level times v. It counts the charge stiffness, not every unit of stress-tensor central charge. The charge norm and representation Casimir cannot be dropped.

[Located obstruction] A lattice boost generated by the energy polarization has overlapping local generators. Its second-order energy cost contains nested commutators of h, of scale M³. Identifying its low-energy action with a conformal diffeomorphism requires projection and velocity normalization. The same stress-density improvements and high-energy matrix elements enter. No dimension-independent estimate for the projected anomaly coefficient follows from the unprojected variational norm.

[Proved scale mismatch] Cutting a periodic ring or joining block ground states changes O(1) bonds, costing O(M). The desired Casimir energy is O(A/L). Comparing these only gives A≲ML at a size already in the conformal regime. Choosing that size requires the conformal onset length. It does not yield A≤constant M independent of the window.

### 6.4 Lieb-Robinson and local current bounds

[Proved implication of a velocity bound] A Lieb-Robinson estimate v_a≤C_LR M gives A≤C_LR cM. It controls the fastest mode, not the number of modes. No particular optimized value of C_LR is needed for this conclusion.

[Proved local algebra] The bond energy current j=i[h_{12},h_{23}] has ||j||≤2M² for a centered representative. A chiral CFT thermal state carries πcT²/12, independent of its velocity, so a current norm bounds cT², not A. Identifying such a state on the lattice at a finite T is an additional conformal-window assumption.

The fusion family has small v=2π/N and growing c. It is fully consistent with these bounds. Its large A is not caused by a superluminal mode.

## 7. Other classes and stress tests

### 7.1 Multicomponent Luttinger liquids

[Proved conditional statement] Write a complete r-component Gaussian low-energy Hamiltonian as

\[
H=\frac1{2\pi}\int dx\big[(\partial\theta)^TV_J(\partial\theta)
 +(\partial\phi)^TV_N(\partial\phi)\big],
\]

with positive matrices V_J,V_N and microscopic on-site charges normalized to the densities -∂φ/π. Then D=V_J/π, κ=V_N^{-1}/π, and the velocities are the positive square roots of the eigenvalues of V_JV_N. Thus

\[
A=\operatorname{tr}\sqrt{V_J^{1/2}V_NV_J^{1/2}}.
\]

Let F_ab=〈-[Q_1^a,[Q_1^b,h]]〉, with its symmetric form understood. The matrix f-sum gives V_J≼πF. It does not upper-bound V_N=κ^{-1}/π. Even a scalar vK bound admits arbitrarily large v if K is allowed to decrease. A susceptibility lower bound, not a susceptibility upper bound, would be needed to close this inference.

For a useful sufficient condition, define the positive Luttinger matrix K by K V_N K=V_J. Then V_J^{1/2}K^{-1}V_J^{1/2} is positive, and its square is V_J^{1/2}V_NV_J^{1/2}. Therefore

\[
A=\operatorname{tr}(V_JK^{-1}).
\]

If K≽κ_0 I, it follows that

\[
A\leq\frac\pi{\kappa_0}\operatorname{tr}F.
\tag{12}
\]

For commuting on-site charges of spectral diameters Δ_a,

\[
\operatorname{tr}F\leq M\sum_a\Delta_a^2.
\]

Indeed, center each charge at the midpoint of its spectrum and apply the double-commutator norm bound 4||Q||²M. A symmetry average removes nonconserving coboundary gauges, and a conserving on-site coboundary commutes with the charges, so this holds with the optimized M. Hence a complete Gaussian theory with

\[
K\succeq\kappa_0 I,\qquad \sum_a\Delta_a^2\leq2\kappa_0
\]

obeys A≤2πM. A single unit-diameter charge and K≥1/2 is an immediate example. These are explicit sufficient hypotheses, not a universal inequality coupling all Luttinger modes. Charge basis and normalization must be kept fixed throughout.

[Proved independent dimension bound] If all r Gaussian modes are resolved by independent commuting on-site densities, r≤d-1 and c=r≤d-1, as recorded in the project. This mode-counting result does not determine their velocity sum.

### 7.2 On-site WZW currents and cosets

[Proved within the supplied current normalization] For one irreducible on-site representation R of a compact simple algebra,

\[
D=-\sum_a[t_1^a,[t_1^a,h]],\qquad
kv\dim\mathfrak g\leq2\pi\langle D\rangle.
\]

With Σ_a(t^a)²=C_2(R)I, write Φ(X)=Σ_at^aXt^a. For centered -MI≤h≤MI, positivity gives Φ(h)≤C_2(R)MI, and

\[
D=2\Phi(h)-2C_2(R)h\leq4C_2(R)MI.
\]

For a Sugawara-complete sector,

\[
A\leq\frac{8\pi C_2(R)}{k+h^\vee}M.
\tag{13}
\]

This proves the proposed 2π bound when 4C_2(R)≤k+h∨, or when the sharper expectation estimate 〈D〉≤(k+h∨)M is independently known. For several complete factors one must sum the actual expectation inequalities; assigning each factor the entire norm budget repeatedly is not an additive proof.

[Proved obstruction] In the counterexample R=Sym^k fundamental of SU(N),

\[
C_2(R)=\frac{k(N-1)(N+k)}{2N}.
\]

Equation (13) becomes A≤4πk(N-1)M/N. It grows with k and is compatible with the violation. Replacing this representation factor by a universal constant would be false, even though this is a Sugawara-complete theory.

For WZW times a neutral coset, the current f-sum sees only the charged Sugawara sector. It supplies no bound on A_coset. A separate microscopic estimate for the neutral sector is necessary. This was an obstruction in the supplied c-versus-d work and remains one here.

### 7.3 Graded swaps and the supersymmetric t-J distinction

[Proved] The mixed graded permutation is an ordinary Hermitian involution on C^{m+n}⊗C^{m+n}, with eigenvalues ±1 and norm one. For H=+ΣP_g and n≥1, a uniform fermionic-color product attains the bond minimum -1. Change one site to a bosonic color. The one-defect sector has diagonal energy shift 2 and hopping 1, hence excitation energy 2+2cos q, quadratic at q=π. For H=-ΣP_g and m≥1, use the uniform bosonic product instead; the defect dispersion is 2-2cos q, quadratic at q=0. Thus the mixed zero-field graded-swap chain has z=2 branches and is outside the stated pure CFT class. The cases with all colors of one grade reduce to ordinary swaps with a possible sign.

[Scope not resolved here] A finite-density t-J chain with chemical potentials is a different Hamiltonian. Its spin-charge CFT cannot be excluded by the preceding argument. Chemical potentials must enter both the thermodynamic ground state and the two-site norm, as half of the on-site term on each bond. A fixed-density spectrum with the chemical-potential term omitted from the norm is not the class in the question. No complete finite-density graded-swap survey was performed after the ungraded fusion counterexample settled the universal claim.

### 7.4 Orthogonal and symplectic rational bonds

[Proved local algebra] For the SO(N) vector rational bond, N≥3,

\[
h=P-\frac{2}{N-2}K,\qquad K=N|\Omega\rangle\langle\Omega|,
\]

the symmetric-traceless, antisymmetric and singlet eigenvalues are 1,-1,-(N+2)/(N-2). Thus M=N/(N-2), also gauge optimal by reflection symmetry. With the usual SO(N)_1 Bethe identification c=N/2 and v=2π/(N-2), its ratio is 1/2. This last IR input is assumed, not independently derived for every N here. At N=3 the bond is TB plus a scalar, so the velocity matches the explicit calibration in section 2. At N=4 the two SU(2) factors provide the corresponding normalization check.

For the Sp(2n) fundamental rational bond h=P+K/(n+1), with K=2n times the symplectic singlet projector, the symmetric, antisymmetric-traceless and singlet eigenvalues are 1,-1,(n-1)/(n+1). Hence M=1 for n≥2. Its full IR content and velocities were not established here. The project already warns against assigning a pure Sp(4)_k theory on the basis of on-site symmetry alone. No unsupported pure-WZW substitution is used to claim a Casimir bound for this branch.

### 7.5 Stacks and emergent cells

[Proved for reflection-symmetric bonds] For a decoupled stack h=Σ_a h_a acting on independent tensor factors, the canonical two-site spectra add, so their spectral widths add. If every h_a is reflection symmetric, so is the stack. Averaging any coboundary with full bond reflection removes it, proving

\[
M_{\rm stack}=\sum_aM_a,\qquad A_{\rm stack}=\sum_aA_a.
\]

The ratio is a weighted mean of the factor ratios. Stacking the free, ULS or TB examples cannot manufacture a violation. Stacking counterexamples preserves one. For general reflection-asymmetric bonds the upper inequality M_stack≤ΣM_a is immediate, but equality after allowing entangling on-site coboundaries is not proved here. It should not be asserted without that qualification.

[Heuristic, using the supplied construction's hierarchy] For frozen cells of length b with effective SU(N) exchange J, physical velocities multiply the logical velocity by b, giving A≈2πbJ(1-1/N). The recorded selector condition J≪t/b^4 therefore makes A/M small, of order at most b^{-3} if M is of order t. The first-order Majorana route is covered at the logical level by section 5 and is additionally reduced by the port visibility. Controlled weak cell couplings do not approach the new fusion violation.

The previous argument that emergent cells cannot violate the conjecture used, in part, a logical-chain Casimir bound as an input. That is valid for its specified free or fundamental-SU(N) targets. It is not a universal induction principle now that other logical targets violate that bound. No new assertion about fixed-d unbounded central charge follows from the fusion family, whose d grows with N and k.

## 8. Numerical survey and reproduction

All computations use Python 3, numpy and scipy, with standard-library file/JSON utilities. No network, external data, GPU or cluster was used. Each invocation below is a short single-machine calculation. Run inside this folder. Thread limits are recommended for predictable runtimes.

### 8.1 Local algebra and dilogarithms

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 fusion_checks.py
```

[Numerical] Output is `fusion_checks.json`. It independently constructs Schwinger generators, diagonalizes X and h for (N,k)=(2,1),(2,2),(2,3),(3,2),(7,2),(8,2), and tests the polynomial bond at k=2. The largest two-site matrix is 1296×1296. X spectral errors are below 8×10^-15. Braided Yang-Baxter residuals at (2,1),(2,2),(3,2) are below 8×10^-16. A centered numerical derivative of R agrees with h to 2×10^-10. The finite Rogers sum agrees with k(N²-1)/(N+k) to 2×10^-13 through N=k=16.

| N | k | d | c | M | A/(2πM) |
|---|---|---|---|---|---|
| 2 | 1 | 2 | 1 | 1 | 0.5 |
| 2 | 2 | 3 | 1.5 | 1.5 | 0.5 |
| 3 | 2 | 6 | 3.2 | 1.5 | 0.711111 |
| 7 | 2 | 28 | 32/3 | 1.5 | 64/63 |
| 8 | 2 | 36 | 12.6 | 1.5 | 1.05 |
| 5 | 5 | 126 | 12 | 137/60 | 144/137 |
| 8 | 8 | 6435 | 31.5 | 2.717857 | 1.448752 |
| 16 | 16 | 300540195 | 127.5 | 3.380729 | 2.357110 |

These entries evaluate the analytic formulas. They are not eight independent ED determinations of CFT data.

### 8.2 TBA: two estimators and cutoff control

`fusion_tba.py` is a documented local adaptation of the supplied ULS solver. It moves the drive and the free-energy probe to string k and uses (5) for the bulk energy. It solves (6) by Anderson iteration and obtains entropy from an independently solved temperature tangent equation. Define

\[
c_f(T)=-\frac{6v[f(T)-e_0]}{\pi T^2},\qquad
c_s(T)=\frac{3vs(T)}{\pi T}.
\]

The first is a free-energy estimator; the second differentiates the solution through the tangent equation. They share the TBA equations and are not statistically or logically independent evidence for Bethe completeness. They do test different temperature combinations and expose crossover corrections. A, in these common-velocity chains, is inferred as c_f v or c_s v. No many-velocity single-v fit is made.

```
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 2 --k 2 --T 100 0.1 0.04 0.02 --output tba_calibration_tb.json
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 3 --k 1 --T 0.02 0.005 --output tba_calibration_uls.json
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 8 --k 2 --T 100 0.04 0.02 0.01 0.005 --output tba_violation_n8k2.json
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 8 --k 2 --M 48 --K 1024 --R 48 --T 0.01 0.005 0.002 --output tba_violation_convergence.json
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 8 --k 2 --M 64 --K 1024 --R 48 --T 0.002 --output tba_violation_cutoff64.json
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 5 --k 5 --M 40 --K 1024 --R 40 --T 100 0.01 0.003 0.001 --output tba_violation_n5k5.json
OPENBLAS_NUM_THREADS=1 python3 fusion_tba.py --N 8 --k 8 --M 48 --K 1024 --R 48 --T 0.003 0.001 --output tba_violation_n8k8.json
```

Files append results when rerun. M in the command line is the string cutoff, unrelated to the bond norm M in the mathematics. K is the number of rapidity grid points and R the half-length of the rapidity interval. The default grid is 512 points with string cutoff 24.

[Numerical] Selected results are

| N,k | T | string cutoff, grid, half-length | c_f | c_s | limiting c |
|---|---|---|---|---|---|
| 2,2 | .02 | 24,512,16 | 1.504004 | 1.504968 | 1.5 |
| 3,1 | .005 | 24,512,16 | 2.001551 | 2.001878 | 2 |
| 8,2 | .005 | 24,512,32 | 12.651737 | 12.668871 | 12.6 |
| 8,2 | .005 | 48,1024,48 | 12.632709 | 12.650896 | 12.6 |
| 8,2 | .002 | 48,1024,48 | 12.614875 | 12.619620 | 12.6 |
| 8,2 | .002 | 64,1024,48 | 12.613358 | 12.618197 | 12.6 |
| 5,5 | .001 | 40,1024,40 | 12.026884 | 12.030203 | 12 |
| 8,8 | .003 | 48,1024,48 | 31.696772 | 31.745567 | 31.5 |
| 8,8 | .001 | 48,1024,48 | 31.634722 | 31.650391 | 31.5 |

At N=8,k=2,T=.002 on the refined grid, the corresponding ratios are 1.051240 and 1.051635. Increasing only the string cutoff further from 48 to 64 gives ratios 1.051113 and 1.051516 at the same temperature. They approach the analytic 1.05 while remaining well above one. Equation residuals and tangent residuals for the displayed runs are below 10^-9. These residuals are solver errors, not bounds on truncation or low-temperature errors. Marginal corrections and the finite string cutoff account for the remaining drift; numerical extrapolation is not used to prove the exact rational excess.

For the larger excess at N=k=8, the two T=.001 ratios are 1.454948 and 1.455668, approaching the analytic 1.448752. This run is corroboration, not a string-cutoff extrapolation at k=8.

At T=100, the entropies are 1.098567 for (2,2), 3.583458 for (8,2), and 4.836214 for (5,5), approaching log 3, log 36 and log 126. This independently checks the representation dimension and the source location. `fusion_checks.json` also records the infinite-temperature local trace of h for a check of the constant term in f(T).

### 8.3 Independent lattice ED calibration

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ed_calibration.py
```

[Numerical] Output is `ed_calibration.json`. The script directly builds periodic microscopic Hamiltonians in conserved charge sectors, finds ground states, and records both

\[
c_{\rm Cas}(L)=-\frac{6L[E_0(L)-e_0L]}{\pi v}
\]

and the spatial-entropy fit S(l)=(c/3)log[(L/π)sin(πl/L)]+constant, using l=2,...,floor(L/2). It also fits E_0(L)/L=e_fit-πA_fit/(6L²)+b/L⁴ when three sizes are available.

Sizes are L=8,10,12,14,16 for SU(2) fundamental; L=6,8,10,12 for TB/2; and L=6,9 for SU(3) fundamental. The largest sector is 73789 states. Sparse ground-state residuals are below 1.2×10^-10. The exact local bond convention is generated by `fusion_checks.py`, independently of the TBA implementation.

The spatial-entropy estimates decrease from 1.06866 to 1.03534 for SU(2), from 1.68157 to 1.59127 for TB/2, and from 2.56480 to 2.32555 for SU(3). These are recognizable finite-size overshoots, not precise estimates of the limiting c. The known-bulk Casimir estimates decrease from 1.03023 to 1.01045, from 1.60020 to 1.53481, and from 2.11346 to 2.05608, respectively. Fits with a free bulk energy and a 1/L⁴ correction give c=1.00538 for SU(2) and 1.51761 for TB/2. All data are recorded in the JSON. This calibration is included to expose, rather than hide, the limited reliability of small-size entropy fits. There is no direct lattice ED determination of c at d=28 or d=36 in this work.

### 8.4 Majorana algebra check

```
OPENBLAS_NUM_THREADS=1 python3 majorana_checks.py
```

[Numerical] `majorana_checks.json` records the exact two-site cross-term norm versus half the nuclear norm for m=1,2,3, with discrepancies below 2×10^-14. Random unrestricted Hermitian coboundaries are removed by the two parity projections with zero displayed residual. A direct two-copy matrix calculation agrees with the number-conserving hopping formula and commutes with total f number, both with zero displayed residual. These finite-dimensional checks test the coefficients and gauge extraction; the analytic proof in section 5 supplies the general statement.

## 9. Consequences for the conformal window and local dimension

[Proved from the entropy-window assumption] Let B=Σ_a c_a/v_a. The lower side of the entropy window and s≤log d imply

\[
B\leq\frac{3\log d}{\pi(1-\varepsilon)T_*}.
\]

Cauchy-Schwarz gives

\[
c^2=\left(\sum_a\sqrt{c_av_a}\sqrt{c_a/v_a}\right)^2\leq AB.
\]

In any class with an established A≤C M,

\[
\boxed{c^2\leq\frac{3C}{\pi}\frac{R\log d}{1-\varepsilon},
\qquad R=M/T_*.}
\]

For the full quadratic Majorana class and spin-S TB family, C=π, so

\[
c^2\leq\frac{3R\log d}{1-\varepsilon}.
\]

For any of the conditional classes satisfying C=2π, the requested 6R log d/(1-ε) follows. In the Majorana class c≤d-1 is therefore guaranteed already when T_*/M≥3log d/[(1-ε)(d-1)²], though the separate algebraic mode bound below is stronger.

[Proved independent free mode count] For m complex orbitals per site, a nearest-neighbor BdG symbol has size 2m and Laurent degree one. Multiplying its zero-energy determinant by z^{2m} gives a polynomial of degree at most 4m. Counting linear roots with multiplicity and their particle-hole pairing bounds the total nonchiral central charge by m. Thus c≤m=log_2 d for a full on-site Fock space, or c≤log_2 d with inert additional local states. Extensive flat zero bands are excluded by the CFT hypothesis. This is also the free bound already recorded in the project.

[Scope] An A bound alone supplies no c≤f(d) without a lower bound on T_*/M or on velocities. In the Gaussian charge-resolved class, c≤d-1 comes instead from the rank of commuting on-site densities. In the symmetric fusion family, the elementary dimension comparison in section 4 gives c≤d-1. None of these is a general argument for c≤f(d) in arbitrary chains.

[Conclusion about the proposed consequence] The universal Casimir premise is false, so it cannot prove the universal coefficient-6 conformal-window bound. This report does not demonstrate a violation of that entropy-window bound itself; it could conceivably hold for another reason. The proven entropy inequality, representation-dependent current bound, and all correctly scoped free results remain intact.

## 10. Ranked next steps and work not done

1. Independently audit the symmetric-fusion normalization and SU(N)_k infrared input, preferably against primary references or a separately implemented finite-size Bethe calculation. The local spectrum, k-string kernel cancellation and TBA dimension check make this a concrete audit, not an open-ended counterexample search. No network was permitted in this run, so literature verification remains outstanding.
2. If a direct finite-size counterexample is desired, target the SU(7) or SU(8) rank-two chain using representation-adapted singlet spaces or finite-size Bethe equations. Generic full-Hilbert-space ED is unnecessary and was not attempted. Keep vacuum-sector congruences explicit and compare Casimir energy with another estimator.
3. Retain the quadratic-Majorana theorem and seek sharp constants there. The proven coefficient π is not saturated by the elementary examples, whose best ratio here is 2. A proof or counterexample to A≤2M for the complete quadratic class would be a meaningful narrower problem.
4. Replace the failed dimension-independent conjecture by a representation- or channel-dependent statement. Any universal dimension-dependent coefficient must accommodate growth at least of order log d/log log d along the symmetric-fusion family. A sharp upper bound of that form is not established.
5. Finite-density graded swaps and the full symplectic rational branch remain useful surveys for smaller-d or different mechanisms. They were not completed after the ungraded analytic family settled the requested universal inequality. The supplied emergent-cell construction was assessed using its recorded scale hierarchy, not recomputed microscopically.

No project-level STATE.md update was made because the task expressly prohibits writes outside this folder. `PROGRESS.md` is the resumable record for this contribution.
