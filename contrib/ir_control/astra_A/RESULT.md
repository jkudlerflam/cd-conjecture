# Relevant and marginal operators under lattice translation and parity

[proved] Statements about zeros concern the continuous Bloch symbol, hence the thermodynamic limit. A finite ring samples discrete momenta and need not have any exact zero. Coincident occurrences in the multiset defining \(\delta _2\) give distance zero. Write \(N=2n\), \(d(u)=\min_{z\in\mathbb Z}|u-2\pi z|\). Results below separate algebraic proofs, consequences of the stipulated continuum theory, and predictions about the interacting lattice model. Several requested assertions need qualifications.

## A. Quadratic stability

[proved] With \(\Psi_k=(c_k,c^\dagger_{-k})^T\), the general quadratic Hamiltonian, up to a constant, is
\[
 H=\tfrac12\sum_k\Psi_k^\dagger\mathcal H(k)\Psi_k,\qquad
 \mathcal H(k)=\begin{pmatrix}h(k)&D(k)\\D(k)^\dagger&-h(-k)^T\end{pmatrix},
\]
where \(h=h^\dagger\), \(D(k)=-D(-k)^T\). Thus
\(\tau_x\mathcal H(-k)^T\tau_x=-\mathcal H(k)\).
These conditions follow by exchanging the two creation operators in the pairing term; conversely they give every Hermitian parity-even quadratic form. Exponential locality gives analytic Fourier coefficients. The reference symbol is diagonal with entries \(e_a(k),-e_a(-k)\).

[proved] An explicit sufficient constant is \(c_0=10^{-3}\), for the norm
\[
 \eta=\max\{\sup_k\|B(k)\|_{\rm op},\sup_k\|B'(k)\|_{\rm op}\},
 \qquad B=\mathcal H-\mathcal H_*.
\]
Assume \(\delta=\delta_2>0\). Here is a proof uniform in the hopping amplitudes. For any reference diagonal entry \(f\), with root separation \(D_0\le\pi\), root slope magnitude \(v\), and distance \(d_*\) to its nearest root,
\[
 |f(k)|=4t|\sin((k-p_+)/2)\sin((k-p_-)/2)|\ge v d_*/\pi.
\]
Indeed the chord triangle inequality gives
\(\sin(D_0/2)\le2\sin(d_{\rm farther}/2)\), and
\(\sin(d_*/2)\ge d_*/\pi\).

[proved] Surround each of the \(4n\) roots by an interval of radius \(r=\delta/8\). Outside these intervals the spectral distance from zero is at least \(v_{\min}\delta/(8\pi)>\eta\), so no perturbed zero occurs. Within one interval, separate its vanishing entry \(f\) from the other entries \(F\). They obey
\[
 |f'|\ge(1-\pi/8)v,\quad
 \|F^{-1}\|\le L^{-1},\quad L=7v_{\min}\delta/(8\pi),\quad
 \|F^{-1}F'F^{-1}\|\le {128\pi\over49v_{\min}\delta^2}.
\]
The first bound uses \(\sin k_a\ge D_0/\pi\ge\delta/\pi\).
The last uses \(|f'/f|\le1/d_++1/d_-\), with both distances at least \(7\delta/8\).

[proved] Write the perturbation blocks as \(b_0,b,C\), and put \(R=(F+C)^{-1}\). The zero condition is the scalar Schur complement
\(s=f+b_0-b^\dagger Rb=0\). Its correction has magnitude at most
\(\eta+\eta^2/(L-\eta)\), smaller than \(|f|\) at either endpoint. Differentiating gives
\[
 |s'-f'|\le\eta+{2\eta^2\over L-\eta}
 +{\eta^2\over(1-\eta/L)^2}
 \left({128\pi\over49v_{\min}\delta^2}+{\eta\over L^2}\right)
 <0.002v_{\min}.
\]
Here \(\eta<10^{-3}v_{\min}\delta\) and packing gives \(\delta\le\pi/2\). Thus \(s\) is strictly monotone, changes sign, and has one simple zero. Its kernel is one-dimensional, and the eigenvalue slope at the zero is \(s'/(1+\|Rb\|^2)\). This proves exactly \(4n\) simple zeros, preserving their chiralities.

[derived, bosonization input] Particle-hole redundancy halves this count: there remain \(n\) right and \(n\) left complex modes, with total central charge \(c=n\).

[proved] The opposite-chirality normal form is
\[
 H(q)=\begin{pmatrix}v(q-\delta/2)&h\\\bar h&-v(q+\delta/2)\end{pmatrix},
 \qquad q^2=\delta^2/4-|h|^2/v^2.
\]
Both zeros disappear for \(|h|>v\delta/2\); at equality they merge non-simply. It describes hybridization or an electron-hole pairing block, with its particle-hole companion understood. For equal co-propagating slopes, replace the second diagonal entry by \(v(q+\delta/2)\): the roots are instead \(q=\pm\sqrt{\delta^2/4+|h|^2/v^2}\). They separate. More generally a positive-definite velocity matrix cannot annihilate two co-propagating crossings locally. Thus the scale is sharp as a worst-case criterion, not as a necessary condition for every perturbation.

[proved] An exact nearest-neighbor periodic realization has diagonals
\(w\sin(k-k_0-\delta/2)\) and \(-w\sin(k-k_0+\delta/2)\), with constant off-diagonal \(h\). Its determinant is
\(w^2[\sin^2(\delta/2)-\sin^2(k-k_0)]-|h|^2\).
Taking \(k_0=\pi/4\) and \(0<\delta<\pi/4\) makes the full BdG reference spacing \(\delta_2=\delta\). Both opposite-chirality pairs disappear when
\(|h|>w\sin(\delta/2)=v_{\min}\tan(\delta/2)\).
This supplies a lattice example with the same asymptotic threshold and no linearization error.

## B. Complete low-dimension list

[derived, bosonization input] Choose \(0<\eta_0<1/4\), and write
\(b=\sup_{\|G\|\le g}\|\mathrm{Mmat}(G)-I\|\).
It suffices that
\[
 b\le\eta_0/2,\qquad 3(1-b)>2+\eta_0.
\]
If \(b\le C_Mg\) for \(g\le g_*\), take \(g\le\min(g_*,\eta_0/(2C_M))\). The input specifies neither \(C_M\) nor \(g_*\), so a numerical bound on \(g\) cannot be extracted.

[proved] For parity-even integer \(Q\), \(\sum Q_i^2\equiv\sum Q_i=0\pmod2\). The only squared norms below six are zero, two, and four. Their partitions give exactly
\[
 0;\qquad \sigma e_i+\tau e_j\ (i\ne j);\qquad
 \pm2e_i;\qquad \sum_{i\in I}\sigma_i e_i\ (|I|=4),
\]
with all signs independent. There are no additional vectors from repeated-index cancellations: cancel first, then classify the resulting charge.

[derived, bosonization input] Descendant levels are nonnegative integers. The allowed levels are \(0,1,2\) for \(Q=0\), \(0,1\) for squared norm two, and \(0\) for squared norm four. These entries stay below \(2+\eta_0\) throughout the specified ball. Every other charge has dimension at least \(3(1-b)>2+\eta_0\); the next oscillator levels are also excluded. The neutral entries are the identity, currents, \(J_iJ_j\), and \(\partial J_i\), including stress tensors. Spatial total derivatives vanish in the integrated Hamiltonian.

[derived, bosonization input] Pauli exclusion does not delete \(2e_i\). Point splitting gives
\(\psi_i(x)\psi_i(x+a)=a\psi_i\partial_x\psi_i+O(a^2)\), a nonzero dimension-two field. The lattice operator \(c_{j,a}c_{j+1,a}\) contains it with coefficient proportional to \(e^{ip_i}\). Distinct-charge bilinears and four-distinct-charge products are likewise nonzero. If isolation of branches is needed, a finite lattice stencil whose Fourier polynomial interpolates prescribed values at the distinct Fermi points supplies it.

## C. Translation selection

[proved] The required margins are precisely
\[
 \delta_2=\min\{d(2p_i),\ d(p_i+\sigma p_j):i<j,\ \sigma=\pm1\},
\quad
 \delta_4=\min_{|I|=4,\,\sigma_i=\pm1}d\left(\sum_{i\in I}\sigma_ip_i\right).
\]
For \(N<4\), set \(\delta_4=+\infty\). These reproduce the multiset definition of \(\delta_2\). Bilinears and their descendants require \(p_i\pm p_j\ne0\); doubled charges require \(2p_i\ne0\); four-distinct charges require the corresponding signed sum nonzero, all modulo \(2\pi\). Hermitian conjugation pairs \(Q\) and \(-Q\) without relaxing momentum conservation.

[derived, bosonization input] At Gaussian points whose Fermi points satisfy these conditions, the only uniform operators below the cutoff are neutral. Apart from the identity, the relevant ones are currents, which shift Fermi points. Marginal neutral terms change velocities and the current quadratic form; derivatives are redundant. Completing the square in the bosonic Hamiltonian removes linear currents, and a positive quadratic Hamiltonian remains Gaussian, with \(c=n\).

[heuristic] This classifies the continuum Hamiltonian through dimension two. A theorem establishing an interacting lattice phase would additionally have to control matching, irrelevant terms, the moving Fermi points, and the RG for arbitrary quasi-local perturbations. The supplied dictionary is not such a theorem.

## D. Oscillations, moving resonances, and measure

[heuristic] Use lattice spacing one, UV energy \(E_0\), and mismatch cutoff \(E_\delta\sim v\delta\). An operator of dimension \(2-x\) has dimensionless coupling
\[
 \lambda(E_\delta)\simeq\lambda_0(E_0/E_\delta)^x,
 \qquad \ell_\delta=\log(E_0/E_\delta).
\]
Rescaling grows the coupling by \(e^{x\ell}\) until the oscillating phase varies on the running cutoff length; subsequent local contributions average away. For dimension one, \(\lambda_0=\epsilon/E_0\), giving the sufficient criterion \(\epsilon\ll v\delta\). Part A disproves the proposed universal “iff.” For near-marginal operators the criterion is
\(\lambda_0e^{x_+\ell_\delta}\ll1\); sufficient stronger conditions are \(\lambda_0\ll1\) and \(x_+\ell_\delta\ll1\). The logarithm must contain a dimensionless ratio, not literally \(v/\delta\) when \(\delta\) is momentum.

[heuristic] Dressed points move by \(O(\epsilon/v_{\min})\), or \(O(\epsilon)\) at fixed velocities. Consequently a bare \(\delta_4\ll\epsilon\) does not protect the dressed margin. Checking only bare resonances is inadequate.

[proved] A quantitative measure statement follows under the stated regularity assumption. Let dressed points be a periodic lift \(F(p)=p+u(p)\), with \(\|Du\|\le b_*<1\). This is a bi-Lipschitz torus diffeomorphism, and its inverse volume distortion is at most \((1-b_*)^{-N}\). For each four-distinct charge, the torus character \(Q\cdot y\) is uniformly distributed: integrate any coordinate with coefficient \(\pm1\). Its strip \(d(Q\cdot y)\le r\) therefore has normalized measure \(r/\pi\). The union bound gives
\[
 \operatorname{Prob}\{\text{some dressed four-index mismatch}\le r\}
 \le{8\binom N4\over\pi}(1-b_*)^{-N}r.
\]
The factor eight removes overall-sign duplication. On bounded parameter charts the same argument has the corresponding Jacobian and boundary constants. The map \((\theta_a,k_a)\mapsto(p_{a,+},p_{a,-})\) has absolute determinant two per flavour; using \(\mu_a\) requires staying away from band edges.

[heuristic] If also \(g=O(\epsilon)\), then excluding \(r=e^{-c/\epsilon}\) leaves near-marginal couplings of size \(O(\epsilon e^{O(c)})\), small for fixed \(c\). This conclusion does not follow if \(g\) is fixed independently of \(\epsilon\). For fixed sufficiently small perturbation, strict dressed-margin inequalities define an open, large-measure sufficient stability region, within the region already satisfying the quadratic criterion. These properties are compatible. They neither establish that the entire physical gapless phase is open nor classify the resonant surfaces. The measure estimate is proved; lattice stability remains a Wilsonian prediction.

## E. Optimal orders and numerical examples

[proved] Packing the \(4n\) points of \(X\) gives \(\delta_2\le\pi/(2n)\). For \(n\ge2\), sort the \(2n\) points and take the \(n\) disjoint successive gaps \(h_j=p_{2j}-p_{2j-1}\). Their sum is at most \(2\pi\). If their minimum pairwise difference is \(d\), their sorted sum is at least \(dn(n-1)/2\). Differences of two such gaps use four distinct indices, hence
\[
 \delta_4\le {4\pi\over n(n-1)}\le {8\pi\over n^2}.
\]
There is no finite four-index upper bound for \(n=1\), since that class is empty.

[proved] The following construction attains both optimal orders. Choose an odd prime \(N\le q\le16N\), and define
\[
 A_i=2qi+(i^2\bmod q),\quad 0\le i<N,\qquad
 d_A=\min_i(A_{i+1}-A_i),\quad B=\lceil d_A/2\rceil,
\]
\[
 T=2A_{N-1}+4B+2,\qquad p_i={2\pi\over T}(A_i+B+1/6).
\]
For \(N\ge4\), \(d_A\ge q+1\), \(T\le6qN\le96N^2\). The latter follows from
\(A_{N-1}<2q(N-1)+q\) and
\(d_A\le A_{N-1}/(N-1)<2q+q/(N-1)\).
The points lie in \((0,\pi)\). Their internal gaps, and both gaps joining their positive and negative copies, are at least \(2\pi d_A/T\). Therefore
\[
 \delta_2\ge{\pi\over3N}={\pi\over6n}.
\]

[proved] The integers \(A_i\) are Sidon. Equality of two pair sums forces equal index sums, since the residue-sum difference has magnitude less than \(2q\). Reducing modulo \(q\) then gives equal sums of squares, hence equal products because two is invertible. The roots of the resulting quadratic give the same unordered pair. Also \(2A_{N-1}<T\), so a two-plus/two-minus sum cannot wrap by \(T\). For other sign patterns, the shift \(1/6\) contributes \(\pm1/3\) or \(\pm2/3\) modulo integers. Thus
\[
 \delta_4\ge {2\pi\over3T}\ge {\pi\over144N^2}={\pi\over576n^2}.
\]
This proves the optimal \(n^{-2}\) order jointly with \(n^{-1}\) spacing; the constants are not optimized.

[proved] For completeness, the prime choice needs no unproved prime-density assumption. Put \(\vartheta(x)=\sum_{p\le x}\log p\). Primes between \(\lceil s/2\rceil\) and \(s\) divide \(\binom{s}{\lfloor s/2\rfloor}\), so induction gives \(\vartheta(s)\le3s\log2\). Also
\[
 16N\log2-\log(16N+1)
 \le\log\binom{16N}{8N}
 \le\vartheta(16N)+4\sqrt N\log(16N).
\]
For the last bound, each prime exponent in a binomial coefficient is at most \(\lfloor\log_p(16N)\rfloor\); contributions beyond the first come only from primes at most \(4\sqrt N\). If no prime lay in \((N,16N]\), these inequalities would imply
\(13\log2\le\log(16N+1)/N+4\log(16N)/\sqrt N\).
For \(N\ge8\) the right side decreases and is already less than eight at eight, whereas the left side exceeds nine. For \(2\le N<8\), use primes \(3,3,5,5,7,7\).

[heuristic] Exhaustive floating-point enumeration in `margins.py` gives the prime example
\[
 \delta_2\simeq4.5152348\,10^{-5},\quad
 \delta_4\simeq2.11928\,10^{-8}.
\]
Attainers are \(p_{31,-}-p_{34,-}\) and
\(p_{17,-}-p_{18,-}-p_{35,-}+p_{37,-}\), respectively, both nearest to zero rather than another multiple of \(2\pi\). The search covers every distinct quadruple and eight sign choices, with a separate small-case enumeration check. These decimals are numerical, not interval-certified. On this machine NumPy `longdouble` has only 53-bit significand precision.

[proved] For the designed example's rational scaffold, \(N=80,q=83,A_{79}=13130,d_A=99,B=50,T=26462\). Pair \(p_{a,-}=p_{a-1}\), \(p_{a,+}=p_{a+39}\), and set \(t_a=1\), \(k_a=(p_{a,+}-p_{a,-})/2\), \(\theta_a=(p_{a,+}+p_{a,-})/2\). The exact scaffold margins are
\[
 \delta_2={198\pi\over26462}=0.02350673968\ldots,\qquad
 \delta_4={2\pi\over79386}=0.00007914727165\ldots.
\]
The earlier lower bounds are saturated by
\(p_{10,+}-p_{11,+}\) and
\(p_{1,-}-p_{31,-}-p_{12,+}-p_{38,+}\), the latter equal to
\(-2\pi-2\pi/(3T)\). Moreover
\(k_a\ge\pi(80q-(q-1))/T=6558\pi/26462>\pi/6\), and the analogous upper bound is below \(5\pi/6\), proving \(v_a\ge1\).

[proved] The actual numerical design adds \(10^{-10}\sqrt{r_{i+1}}\) to each scaffold point, where \(r_j\) is the \(j\)-th prime. Thus it does not sit on the rational grid. If \(H=10^{-10}\sqrt{409}<2.023\,10^{-9}\), its margins decrease by at most \(2H\) and \(4H\), respectively. The velocity bound survives, since the added shifts increase with index and change each \(k_a\) by less than \(H\).

[heuristic] Exhaustive enumeration for this perturbed design gives
\(\delta_2\simeq0.0235067396931\), attained by \(p_{10,+}-p_{11,+}\), and
\(\delta_4\simeq0.00007914198998\), attained by
\(p_{26,-}+p_{38,-}+p_{8,+}+p_{9,+}\), just below \(2\pi\).
The minimum velocity is \(1.40705998\). This construction optimizes the asymptotic orders, not the finite-\(40\) constants globally. Exact formulas and computed parameters are in the script and `numerical_results.json`.

## F. Parity-odd spin operators

[derived, bosonization input] In the stipulated half-integer twist sectors,
\[
 \Delta_0(Q)\ge\tfrac12N(1/2)^2=n/4,
\]
with equality exactly when every \(Q_i=\sigma_i/2\). The usual Jordan–Wigner string has neutral half charges \((1/2,-1/2)\) per flavour, or its other branches. Its background phase is \(\pi\rho_a=k_a\), since the occupied interval has length \(2k_a\). Multiplying by local odd fermion products selects
\(\sum_iQ_i\in2\mathbb Z+1\). Thus not every half-integer vector is a local odd spin field. The minimum remains attainable: make one flavour's two signs equal and all other flavour pairs opposite.

[derived, bosonization input] Including that string background, and using the lifts \(p_{a,\pm}=\theta_a\pm k_a\), the lightest momenta are
\[
 P_\sigma=\tfrac12\sum_a[(\sigma_{a,+}+\sigma_{a,-})\theta_a
 +(\sigma_{a,+}-\sigma_{a,-})k_a]\pmod{2\pi},
 \quad \tfrac12\sum_i\sigma_i\ \text{odd}.
\]
For example, at \(n=1\) the local odd lightest fields have momentum \(\pm\theta_1\), not \(\pm k_1\). Applying the integer-charge rule blindly to all half vectors would include disorder sectors. The half-charge dimension formula and string identification here use the twist-sector assumption explicitly added in F; section 2 alone does not specify twist locality or its translation character.

[derived, bosonization input] All odd fields are irrelevant by dimension for \(n\ge9\); \(n=8\) is marginal and \(n\le7\) has relevant twist fields before translation selection. Uniform irrelevance on the Gaussian ball follows if
\[
 b<1-8/n,\qquad
 \Delta_G-2\ge n(1-b)/4-2>0.
\]
Thus \(g<(n-8)/(nC_M)\), within \(g_*\), suffices. Generic momenta can forbid the lightest twists even for smaller \(n\); then all other candidate low-dimension twist charges must also be checked. No translation-based necessity of \(n\ge9\) is claimed.

## G. Coincident Majorana Fermi points

[proved] Normalize \(\chi_R^a(z)\chi_R^b(0)\sim\delta_{ab}/z\) and similarly on the left, with velocity scaled to one. Currents \(J^{ab}=i\chi^a\chi^b\), \(a<b\), obey
\(J^A(z)J^B(0)\sim\delta^{AB}/z^2+if^{ABC}J^C(0)/z\).
Translation fixes right movers and negates left movers. Thus all mixed marginal terms are \(g_{AB}J_R^AJ_L^B\), with real unrestricted \(g_{AB}\). Chiral dimension-two terms are allowed too: chiral four-fermion products and derivative bilinears, modulo total derivatives. They include stresses and velocity changes. Mixed bilinear masses carry momentum \(\pi\) and are forbidden; chiral currents of dimension one are allowed.

[proved] For \(S_{\rm int}=(2\pi)^{-1}\int g_{AB}J_R^AJ_L^B\), the one-loop equation is
\[
 {dg_{EF}\over d\ell}=\tfrac12\sum_{A,B,C,D}
 f^{ACE}f^{BDF}g_{AB}g_{CD}+O(g^3).
\]
The product OPE has coefficient \(-f^{ACE}f^{BDF}/|z|^2\); the negative second cumulant and shell integral \(2\pi d\ell\) give the displayed sign and factor. Since \(\sum_{AC}f^{ACE}f^{ACF}=2(m-2)\delta_{EF}\), the isotropic equation is \(\dot g=(m-2)g^2\). Its irrelevant sign is \(g<0\).

[proved] For \(m\ge4\), instability is explicit in the one-loop ODE. Split the Majoranas into groups of sizes two and \(r=m-2\). Let diagonal current couplings be \(a,b,c\) within the first group, within the second, and across groups. Direct substitution gives
\[
 \dot a=rc^2,\quad \dot b=(r-2)b^2+2c^2,\quad
 \dot c=[a+(r-1)b]c.
\]
Start at \(a=g+h,b=c=g<0\), with arbitrarily small \(0<h<|g|\). While \(a<0\), put \(A=-a,B=-b,C=-c\). The region \(B-A\ge h,C\ge B\) is invariant because
\[
 (B-A)'=(r-2)(C^2-B^2)\ge0,\qquad
 (\log(C/B))'=-A-B+2C^2/B>0.
\]
Also \(B'\le-rh^2\), so \(a\) must cross zero in finite time.

[proved] If \(r>2\) and \(b\) never crosses zero, monotonicity implies \(b\to0\), while \(b'\ge2c^2\) bounds \(\int c^2\) and hence \(a\). But eventually \(a+(r-1)b>0\), forcing \(|c|\) to grow, a contradiction. Once \(a,b>0\), the variables \(a,b,|c|\) have positive derivatives and their minimum satisfies \(\dot z\ge rz^2\), giving finite-time blowup. For \(r=2\), put \(s=(a+b)/2\); then \(\dot s=2c^2,\dot c=2sc\) and
\(s^2-c^2=gh+h^2/4<0\), the explicit runaway Kosterlitz–Thouless trajectory. This proves non-interiority for the one-loop basin. For \(m=2\), the algebra is abelian, the beta function vanishes, and the requested universal assertion is false.

[heuristic] Interpreting runaway as a physical strong-coupling instability, and extending this basin statement beyond perturbation theory, requires additional RG control.

[derived, bosonization input] After detuning, decompose each quartic term into chiral charges. Its neutral components are density products by fermionic reordering. Every nonneutral component has one of the momenta in C and oscillates if the corresponding margins are positive. It drops out below its own mismatch scale, not necessarily the single-particle detuning scale. Distinct Fermi points alone do not suffice: equal \(k_a=k_b\) with distinct \(\theta_a,\theta_b\) leaves flavour exchange momentum \(2(k_a-k_b)=0\).

## H. Operator summary

| Class | Charge; free dimension | Momentum | Exclusion or role | Status |
|---|---|---|---|---|
| Identity | \(0;0\) | 0 | Energy offset | [derived, bosonization input] |
| Currents | \(0;1\) | 0 | Allowed; shift Fermi points | [derived, bosonization input] |
| Neutral level two | \(0;2\) | 0 | Gaussian quadratic forms; derivative redundancies | [derived, bosonization input] |
| Bilinears | \(\pm e_i\pm e_j;1\) | \(\pm p_i\pm p_j\) | \(\delta_2>0\) | [derived, bosonization input] |
| Bilinear descendants | Same charge; 2 | Same | \(\delta_2>0\) | [derived, bosonization input] |
| Same-branch pairs | \(\pm2e_i;2\) | \(\pm2p_i\) | \(\delta_2>0\); survive Pauli exclusion | [derived, bosonization input] |
| Four distinct branches | \(\sum_I\pm e_i;2\) | \(\sum_I\pm p_i\) | \(\delta_4>0\) | [derived, bosonization input] |
| Higher even sectors | Other charges/levels; \(\ge3\) | \(Q\cdot p\) | Dimension bound in B | [derived, bosonization input] |
| Odd spin twists | Half-integer, odd local sector; minimum \(n/4\) | \(P_\sigma\) above | Exact parity, or momentum, or dimension | [derived, bosonization input] |
| Symmetric mixed currents | \(J_R^AJ_L^B;2\) | 0 | Allowed at coincidence; detuning needs C | [proved] |
