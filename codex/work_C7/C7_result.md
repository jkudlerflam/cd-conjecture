C7. Does fixed microscopic dimension bound the conformal central charge?

Closed-book analysis, October 27, 2026.

The unrestricted question remains unresolved in this analysis. I do not have either an unconditional function f(d) or a conformal counterexample at fixed d. There are three substantive results below.

1. A conditional bound c_T <= d^2 - 1 holds when the stress tensor is exhausted by affine currents of microscopic on-site continuous symmetries. No bound on the affine levels or on a conformal onset scale is required.
2. Modular invariance gives an explicit exponential lower bound on the number of CFT states below dimension c_T/4 + 2. This yields an f(d) bound under a stated, uniform finite-volume realization hypothesis. It avoids assuming a Cardy asymptotic that is uniform across theories.
3. An explicit four-state, one-site translation-invariant, nearest-neighbor Hamiltonian has extremal ground states of any prescribed period b >= 5, with couplings independent of ring length. Adding a spectator qubit gives block ground spaces of dimension 2^b at fixed microscopic dimension eight. This solves the block-selection part of a possible counterexample, but not the coupling part. A separate rank theorem rules out obtaining the desired SU(N) permutation interaction by direct projection of one boundary bond when N > d.

The conditional hypotheses in results 1 and 2 are incomparable. I cannot identify a globally weakest natural extra hypothesis. For each proof, I identify the additional input it uses and the term that remains uncontrolled without it.

Conventions and assumptions.

The central charge c below is the usual nonchiral CFT central charge: c_L = c_R = c, not c_L + c_R. Thus c = c_T when several sectors are combined. Logarithms are natural. The lattice spacing is one. The supplied facts about ULS, entropic and thermodynamic obstructions, Lieb–Robinson bounds, and universal Hamiltonians are used as inputs and are not re-derived.

For the modular counting result only, I assume a unitary CFT with a normalizable vacuum, discrete spectrum with finite degeneracies, and a modular-invariant torus partition function finite for every positive aspect ratio. I also assume the absence of a gravitational anomaly, as supplied in the question. A noncompact theory with a divergent torus trace is not covered by that result. For several velocities, the counting result requires that the total scaling-dimension spectrum has such a modular-invariant partition function after forgetting the velocities. Independent compact CFT factors satisfy this condition. A collection of numbers c_a and v_a, without a specified full Hilbert space and modular completion, is insufficient.

The current-algebra result has different hypotheses, stated in Theorem 1. The long-period construction and the boundary-rank obstruction are finite-dimensional Hamiltonian statements and do not assume a CFT.

Theorem 1 (An affine-current completeness bound).

Suppose the independent affine current algebras that exhaust the infrared stress tensor arise from continuous microscopic on-site symmetries. More explicitly, assume the following.

(i) After removing generators that act as a scalar or act trivially in the infrared, their joint compact Lie algebra g embeds into pu(d).

(ii) In each chirality the full stress tensor is the sum of the Sugawara stress tensors of these currents. Equivalently, in the unitary coset decomposition, the commuting coset has central charge zero.

Then

\[
c_T\leq \dim\mathfrak g\leq d^2-1.
\]

Different velocities are allowed for independent factors, provided their joint symmetry algebra satisfies (i) and their Sugawara tensors exhaust all the sectors in (ii). Repeated copies carrying only a diagonal microscopic symmetry do not automatically satisfy these hypotheses.

Proof.

A microscopic on-site symmetry acts on one site through U(d). Its scalar subgroup acts trivially on observables and cannot supply a fluctuating conserved density. Removing it gives a Lie subalgebra of pu(d), of dimension at most d^2 - 1. Passing to the image in the infrared can only decrease the dimension.

A compact Lie algebra is reductive. Write

\[
\mathfrak g=\mathfrak u(1)^r\oplus\bigoplus_j\mathfrak g_j,
\]

with simple compact factors g_j. The supplied Sugawara formula and the abelian current algebra give

\[
c_{\rm Sug}=r+\sum_j\frac{k_j\dim\mathfrak g_j}{k_j+h_j^\vee}
\leq r+\sum_j\dim\mathfrak g_j=\dim\mathfrak g.
\]

Here the nontrivial unitary compact affine levels are positive. Hypothesis (ii) says c_T = c_Sug. Combining the inequalities proves the result. No step uses k_j <= ind(R). In particular, this bound survives arbitrarily large affine levels. QED.

Where the hypothesis enters.

Without (ii), the unitary decomposition instead gives

\[
T=T_{\rm Sug}+T_{\rm coset},\qquad
c_T=c_{\rm Sug}+c_{\rm coset},\qquad c_{\rm coset}\geq0.
\]

The last inequality is the coset positivity input already present in the supplied background. It is a lower-bound mechanism, not an upper bound on the coset. A neutral tensor factor can increase c_coset without changing the specified current algebra. This is a statement about continuum theories, not a construction of such lattice models at fixed d.

For this proof, the weakest quantitative replacement for (ii) is a bound c_coset <= B(d), giving c_T <= d^2 - 1 + B(d). Setting B = 0 is a natural structural condition, current-algebra completeness. A uniform bound on a nonzero coset would be a further unproved hypothesis. On-site symmetry alone, rationality alone, and a statement that the theory contains a WZW factor do not establish (ii).

If the required independent symmetry is on-site only after blocking at most b(d) sites, the same proof gives c_T <= d^{2b(d)} - 1. Allowing an unbounded blocking length loses the conclusion.

Lemma 2 (A universal supply of light CFT states).

Under the modular assumptions above, let N_c(D) count states of total scaling dimension Delta <= D, with multiplicity. If c >= 1, then

\[
N_c(c/4+2)\geq 1+\frac{\pi e}{12}\,c\,e^{\pi c/6}.
\tag{1}
\]

The cutoff proportional to c is important. This is not a statement about a fixed number of conformal levels.

Proof.

Set C = c/12 and E = Delta - C. The torus partition function is

\[
Z(t)=\sum_E e^{-2\pi tE},\qquad Z(t)=Z(1/t).
\]

Differentiation at the self-dual point, and evaluation at t = 2, respectively give

\[
\sum_E E e^{-2\pi E}=0,\qquad
\sum_E(e^{-4\pi E}-e^{-\pi E})=0.
\tag{2}
\]

Trace convergence justifies these operations. Define

\[
A=\frac{C}{2(e^{2\pi C}-e^{-\pi C})},\qquad
F(E)=Ee^{-2\pi E}+A(e^{-4\pi E}-e^{-\pi E}).
\tag{3}
\]

Equation (2) implies sum_E F(E) = 0. The vacuum has

\[
F(-C)=-\frac C2 e^{2\pi C}.
\tag{4}
\]

We need the signs and an upper bound on F.

First, F(E) <= 0 for -C <= E <= 0. To see this, put x = -E and define

\[
r(x)=\frac{x}{e^{2\pi x}-e^{-\pi x}}.
\]

The denominator D(x) vanishes at zero and is strictly convex for x >= 0, since

\[
D''(x)=4\pi^2e^{2\pi x}-\pi^2e^{-\pi x}>0.
\]

Consequently D(x)/x increases and r(x) decreases. For 0 < x <= C, r(x) >= r(C) = 2A, and

\[
F(-x)=(e^{4\pi x}-e^{\pi x})(A-r(x))\leq0.
\]

Second, for E >= 0 the term proportional to A is nonpositive. Therefore

\[
F(E)\leq Ee^{-2\pi E}\leq\frac{1}{2\pi e}.
\tag{5}
\]

Third, F(E) <= 0 for E >= 2C + 2. For positive E, this is equivalent to

\[
A\geq q(E):=\frac{E}{e^{\pi E}-e^{-2\pi E}}.
\]

We have A >= (C/2)e^{-2pi C}. On E >= 2C + 2, the numerator E e^{-pi E} decreases and the denominator 1 - e^{-3pi E} increases. Hence

\[
q(E)\leq
\frac{(2C+2)e^{-2\pi C-2\pi}}{1-e^{-6\pi}}
\leq \frac C2e^{-2\pi C}.
\tag{6}
\]

For the last inequality, C >= 1/12 gives

\[
\frac{4(1+1/C)e^{-2\pi}}{1-e^{-6\pi}}
\leq\frac{52e^{-2\pi}}{1-e^{-6\pi}}<1.
\]

Thus any positive contributions to sum F come from 0 < E < 2C + 2, equivalently C < Delta < 3C + 2. By (4), their sum must be at least (C/2)e^{2pi C}. Each state contributes at most 1/(2pi e) by (5). There must be at least pi e C e^{2pi C} such states, in addition to the vacuum. Substituting C = c/12 proves (1). Extra vacua only strengthen this argument. QED.

Theorem 3 (A finite-volume conditional f(d) bound).

Assume the modular hypotheses of Lemma 2. In addition, suppose there is an integer K(d), independent of the interaction h, such that for every conformal chain one can choose a compatible ring length L <= K(d) on which all CFT states with

\[
\Delta\leq c_T/4+2
\tag{7}
\]

are represented by linearly independent lattice states. This is a finite-volume completeness assumption. It is implied by a sufficiently accurate, multiplicity-preserving spectral realization through this cutoff. The proof only needs the injectivity, not an estimate of the eigenvalue errors.

Then a valid bound is

\[
\boxed{\quad c_T\leq 1+\frac6\pi K(d)\log d.\quad}
\tag{8}
\]

A tighter form for d >= 2 is

\[
c_T\leq\max\left\{1,\frac6\pi
W\!\left(\frac{2(d^{K(d)}-1)}e\right)\right\},
\tag{9}
\]

where W is the nonnegative Lambert W branch. The d = 1 chain has no nontrivial conformal degrees of freedom.

Proof.

For c_T >= 1, injectivity and Lemma 2 give

\[
d^L\geq 1+\frac{\pi e}{12}c_Te^{\pi c_T/6}.
\tag{10}
\]

Writing y = pi c_T/6 converts this to y e^y <= 2(d^L - 1)/e, proving (9). Alternatively, taking logarithms after dropping the 1 in (10), and using c_T >= 1, gives

\[
\frac\pi6c_T\leq K(d)\log d+\log\frac{12}{\pi e}.
\]

The second term multiplied by 6/pi is less than 1. This proves (8). Values c_T < 1 already satisfy it. QED.

What the finite-volume hypothesis means physically.

For a common velocity v, the cutoff (7) corresponds to excitation energies of order

\[
\frac{2\pi v}{L}\left(\frac{c_T}{4}+2\right).
\tag{11}
\]

For independent sectors with different velocities, a state with total Delta <= D has energy at most 2pi v_max D/L. Requiring the continuum spectrum through that physical energy is sufficient for the counting argument. An assertion that only the first few states or only Delta <= D_0 with D_0 independent of c_T have converged does not imply (7).

Preserving a fixed positive fraction of the states counted in Lemma 2 would also suffice, with an additive logarithmic change in the bound. More generally, an exponential loss smaller than exp[(pi/6)c_T] can be tolerated. Thus full spectral equality is stronger than the proof needs. Its indispensable input is a uniform ring length at which an entropy proportional to c_T is already represented.

This is a conditional result, not a derivation of K(d). No supplied fact controls K(d). If the first suitable ring has L proportional to c_T, or much larger, (10) is consistent with unbounded c_T at fixed d.

For comparison, modular invariance also gives a short canonical statement at the self-dual temperature. The CFT entropy at beta = L/v is

\[
S_{\rm CFT}=(1-t\partial_t)\log Z(t)|_{t=1}
=\log Z(1)\geq\pi c_T/6.
\]

Transporting this entropy to a bounded lattice ring would give a similar conditional bound. Without a uniform finite-volume transport assumption, it has the same missing input. This does not repair the thermodynamic obstruction supplied in the question.

Why the four proposed upper-bound routes do not close unconditionally.

(i) Onset scales. The supplied estimates already identify the missing input, a bound on the conformal window in units of the lattice spacing or velocity. A finite interaction norm gives neither. A positive gap in auxiliary degrees of freedom that depends on h is also insufficient without a uniform lower bound relative to the effective critical scale. The explicit construction below illustrates arbitrarily large structural scales at fixed d and bounded interaction norm. It is not itself conformal, so it does not prove that conformal onset scales are unbounded within the critical class.

(ii) Ring counting. The ordered vacuum descendants are spanned by products of L_{-n} with n >= 2. Their number at each level is bounded by a partition count independent of c. For c > 1, the generic Virasoro vacuum module has the character

\[
\chi_{\rm vac}(q)=q^{-c/24}\prod_{n=2}^{\infty}(1-q^n)^{-1}.
\]

Its level multiplicities do not grow with c. The location of the Casimir energy depends on c, but counting those descendants alone cannot bound c. Cardy growth requires control of a high enough spectrum and of its onset across theories. Lemma 2 replaces the latter issue on the CFT side by a modular argument, but Theorem 3 still needs lattice realization through Delta of order c. Ordinary convergence of each fixed low-lying level as L tends to infinity is pointwise convergence; it supplies no bound on the ring length required for this c-dependent collection of states.

(iii) Blocking and renormalization. Exact blocking sends d to D = d^b and rescales the lattice unit. It does not reduce the infrared central charge. A bound involving D at an uncontrolled b does not give a bound involving d. The continuum c-theorem concerns flows between suitable Lorentz-invariant continuum theories; it does not assign the microscopic site an ultraviolet central charge bounded by log d. The supplied ULS examples already rule out the proposed lattice-scale entropic monotonicity argument.

A precise additional assumption that would make this route work is a bounded-entanglement renormalization representation. Suppose the ground state has an exact tensor network, or one with entropy error o(log ell), such that an interval of ell sites can be separated by cutting at most q(d) log_b ell + O(1) bonds of dimension at most chi(d). The cut bounds Schmidt rank, so

\[
S(\ell)\leq\frac{q(d)\log\chi(d)}{\log b}\log\ell+O(1).
\]

Comparing with the CFT interval coefficient S(ell) = (c_T/3) log ell + O(1) gives

\[
c_T\leq\frac{3q(d)\log\chi(d)}{\log b}.
\]

The missing theorem is the bound on chi and the causal width q. A fixed approximation error in trace norm alone does not imply the required o(log ell) entropy error as interval Hilbert spaces grow.

(iv) Operator algebras and compactness. A single-site algebra has dimension d^2, but continuum fields can be limits of operators on arbitrarily many sites, and different low-energy momenta can appear in the same microscopic operator. There is no established injection of the whole continuum field algebra into the single-site algebra. Theorem 1 provides an injection for microscopic on-site symmetry generators, and the coset is the explicit remainder.

After subtracting a scalar and normalizing the interaction, the space of h is finite-dimensional and compact. This alone does not bound a discontinuous infrared invariant on the conformal subset. A sequence of critical interactions could converge to an interaction with vanishing velocities, a collapsed conformal window, or a nonconformal low-energy limit. Upper semicontinuity and finiteness on all accumulation points would be additional assumptions. A finite number of real coupling parameters is not a parameter-counting proof either: a real coefficient can encode unbounded information, and symmetries or exceptional loci can relate many apparent tuning conditions.

Theorem 4 (An explicit arbitrary-period selector at d = 4).

For every integer b >= 5 there is a one-site translation-invariant nearest-neighbor Hamiltonian with local dimension four, interaction norm bounded independently of b, and couplings independent of ring length, whose ground space on every ring L divisible by b has exactly b states. Its extremal infinite-chain ground states have marker period b. For each fixed b the spectral gap above this ground space is bounded below uniformly in these ring lengths.

Construction.

Use the alphabet {#, L, X, R}. The allowed neighboring pairs are

\[
\#L,\ \#X,\ LL,\ LX,\ XR,\ X\#,\ RR,\ R\#.
\tag{12}
\]

Let Q_i be the diagonal projector onto forbidden pairs on bond (i,i+1), and M_i = |#><#|_i. Fix t > 0 and K > 4t, for example K = 8t, and set

\[
H_L=K\sum_i Q_i+mu\sum_i M_i
-t\sum_i\left(|LX\rangle\langle XR|+|XR\rangle\langle LX|\right)_{i,i+1}.
\tag{13}
\]

This has the required form sum_i h_{i,i+1}: assign half of each on-site marker term to each adjacent bond. All indices are periodic. We will take mu = 2t - delta with 0 < delta < 2t. In particular,

\[
\|h\|\leq K+\mu+t<K+3t,
\tag{14}
\]

uniformly in the selected period.

Step 1. The local sectors are invariant.

Represent the four letters as arrows between two auxiliary labels,

\[
\#=(R,L),\quad L=(L,L),\quad X=(L,R),\quad R=(R,R).
\]

The allowed pairs concatenate arrows. The move LX <-> XR changes only the internal label of a two-arrow path with fixed external endpoints. Thus it preserves the allowed/forbidden status of its own bond and both neighboring bonds. All Q_i commute with H, individually. Every marker projector M_i also commutes with H because the hopping move contains no marker. This establishes exact superselection by marker positions and forbidden-bond positions, not a large-K approximation.

Step 2. Solve a legal interval between markers.

Suppose two successive markers are separated by a cell of b sites, counting the first marker but not the second. The b - 1 internal sites must have the form

\[
L^j X R^{b-2-j},\qquad j=0,\ldots,b-2.
\tag{15}
\]

There is exactly one X. The hopping term is the adjacency Hamiltonian of an open path with b - 1 vertices. Its eigenvectors have components proportional to sin[m pi(j+1)/b] and energies -2t cos(m pi/b), m = 1,...,b-1. Therefore the lowest energy of this cell, including its first marker, is

\[
E_b=\mu-2t\cos(\pi/b)
=a_b-\delta,
\qquad a_b=2t(1-\cos(\pi/b)).
\tag{16}
\]

Different cells are independent because no hopping term touches a marker. Their energy densities are

\[
e_b=\frac{a_b-\delta}{b},\qquad b\geq2.
\tag{17}
\]

Step 3. Select an arbitrary integer b uniquely.

Define the crossing values

\[
\delta_b=(b+1)a_b-b a_{b+1}.
\tag{18}
\]

The function a(x) = 2t(1 - cos(pi/x)) is decreasing and strictly convex for x >= 2. Indeed,

\[
a''(x)=\frac{4t\pi}{x^3}\sin(\pi/x)
+\frac{2t\pi^2}{x^4}\cos(\pi/x)>0.
\]

It follows that

\[
\delta_b-\delta_{b+1}
=(b+1)(a_b-2a_{b+1}+a_{b+2})>0,
\]

and delta_b tends to zero. Also delta_b > a_b > 0. Choose

\[
\delta_b<\delta<\delta_{b-1}.
\tag{19}
\]

The adjacent energy-density difference is

\[
e_j-e_{j+1}=\frac{\delta_j-\delta}{j(j+1)}.
\tag{20}
\]

It is positive for j < b and negative for j >= b. Thus e_b is the unique global minimum over all integer cell lengths, and e_b < 0 because delta > a_b.

For b >= 5 the interval (19) lies below 2t. One convenient bound is

\[
\delta_j=a_j+j(a_j-a_{j+1})\leq\frac{3t\pi^2}{j^2},
\]

obtained from 1 - cos u <= u^2/2 and -a'(x) <= 2t pi^2/x^3. With j = b - 1 >= 4 this is less than 2t. Hence mu is positive as promised. The midpoint of (19) is a concrete choice of the coupling for each b.

Step 4. Rule out illegal sectors at finite K.

Let q >= 1 forbidden bonds be present. Their positions are conserved. Cutting at these bonds decomposes the ring into q legal open strings. Inside each string, complete intervals from one marker to the next have energy at least e_b times their length. There are at most two remaining fragments per string, the prefix and suffix; if there are no markers there is only one fragment.

Each remaining fragment has at most one X. Its hopping Hamiltonian is a path, possibly with fixed endpoints or missing moves, and is bounded below by -2t. Any unassigned marker energy is nonnegative because mu > 0. Since e_b < 0, a fragment of length s has energy at least

\[
-2t\geq e_b s-2t.
\]

Adding all complete cells and fragments gives the operator lower bound in this sector

\[
H_L\geq e_bL+(K-4t)q.
\tag{21}
\]

A hopping term on a cut bond vanishes because LX and XR are allowed pairs, so no omitted coupling crosses the cuts. Since K > 4t, illegal sectors cannot compete with a legal tiling that attains e_bL.

In a legal periodic sector with no markers, the only words are all L and all R. They have zero energy and also lose to e_bL.

Step 5. Ground states and the gap.

For L divisible by b, a tiling by cells of length b attains e_bL. Equality requires every cell to have length b and every path Hamiltonian to be in its unique lowest state. There are b distinct translations of the marker pattern. These are orthogonal and exhaust the ground space.

Define

\[
\gamma_b=\min_{k\geq2,\ k\ne b} k(e_k-e_b)>0.
\tag{22}
\]

Positivity follows from unique minimization and the fact that k(e_k - e_b) tends to infinity as k tends to infinity. The first internal excitation of a length-b cell costs

\[
\omega_b=2t\left(\cos(\pi/b)-\cos(2\pi/b)\right)>0.
\tag{23}
\]

Any legal marker sector with a different cell length costs at least gamma_b. Marker-free sectors cost |e_b|L >= |e_b|b, and illegal sectors cost at least K - 4t. Thus on compatible rings

\[
\operatorname{gap}(H_L)\geq
\min\{\gamma_b,\omega_b,|e_b|b,K-4t\}>0.
\tag{24}
\]

This lower bound depends on b but not L. In infinite volume, choosing one marker pattern gives an extremal ground state that breaks translation to period b. Translation eigenstates on finite rings are superpositions of these patterns; the claim is not that such superpositions fail to be translation invariant. QED.

Consistency checks and the limiting interaction.

Every energy above has units of t. The construction is exactly nearest-neighbor and one-site translation invariant. The parameter delta is fixed once b is selected and never depends on L. No perturbative elimination was used. K can be fixed at 8t for all b, so the interaction norm is bounded throughout the family.

There is no conflict with the supplied classical period bound. The term -2t cos(pi/b) is the ground energy of a coherent quantum walk across a cell. It is not the sum of fixed classical edge costs along a word.

For small delta, Taylor expansion gives a_b = t pi^2/b^2 + O(b^{-4}), delta_b = 3t pi^2/b^2 + O(b^{-3}), and the minimizing period grows as pi sqrt(3t/delta). At the midpoint of (19), neighboring-cell excess energies are of order t/b^4, while omega_b is of order t/b^2. There is no uniform auxiliary gap as b increases.

The accumulation point delta = 0 makes the compactness obstruction concrete. Now every legal cell has positive energy a_b, while all-L and all-R have energy zero. A ring with one marker and one length-L cell has excitation energy

\[
a_L=2t(1-\cos(\pi/L))\sim t\pi^2/L^2.
\]

Every legal sector with markers has energy at least a_L, since its cell lengths are at most L. The same fragment argument bounds illegal sectors below by K - 4t. Hence the gap is a_L for sufficiently large L. The limiting Hamiltonian has a quadratic finite-size scale, not the required ordinary finite-velocity conformal spectrum. A sequence in a compact interaction space can therefore develop arbitrarily large structural scales by approaching a nonconformal point. This example makes no assertion about a sequence of conformal central charges.

Corollary 5 (Large block ground manifolds at fixed d = 8).

Tensor every site of Theorem 4 with a qubit and let H act as the identity on the qubits. The microscopic dimension is eight. In a chosen period-b sector, each cell has its unique clock ground vector tensored with b free qubits, giving a cell ground manifold

\[
\mathcal C_b\simeq(\mathbb C^2)^{\otimes b},\qquad N=2^b.
\]

The gap above the entire ground manifold is still bounded below by (24). The manifold is exponentially degenerate in the ring length. It is not a CFT, and these free qubits cannot be counted as massless conformal modes. This corollary provides only an explicit supply of large block spaces with L-independent couplings.

Theorem 6 (A direct-projection obstruction to the proposed ULS encoding).

Partition a nearest-neighbor d-state chain into disjoint contiguous blocks, with at least three blocks on a periodic ring so that each neighboring pair shares one microscopic bond. In each block choose an N-dimensional code with isometry V_i : C^N -> (C^d)^{tensor b_i}. Suppose the desired effective interaction is obtained by directly compressing H to the product of these block codes. If the interaction between two neighboring code spaces is a nonzero multiple of P_N, up to arbitrary one-block terms and a scalar, then N <= d.

Proof.

Choose a Hermitian orthonormal single-site operator basis consisting of the identity and traceless matrices F_alpha, alpha = 1,...,d^2 - 1. Modulo one-site operators and scalars, the microscopic bond is

\[
h_{\rm conn}=\sum_{\alpha,\beta=1}^{d^2-1}
J_{\alpha\beta}F_\alpha\otimes F_\beta.
\]

Compression onto the neighboring blocks maps each boundary operator to V_i^dagger F_alpha V_i. Discard its scalar part, since scalar parts contribute only one-block terms. In traceless operator bases on the two N-dimensional codes, the connected coefficient matrix therefore has the form

\[
J_{\rm eff}=A J B^{\mathsf T}.
\]

Its rank is at most d^2 - 1, regardless of block length or entanglement inside each code.

For an orthonormal Hermitian basis T_A of traceless N by N matrices, the swap satisfies

\[
P_N=\frac1N I\otimes I+\sum_{A=1}^{N^2-1}T_A\otimes T_A.
\]

Its connected coefficient matrix has rank N^2 - 1. A nonzero coefficient multiplying it does not change the rank. Thus N^2 - 1 <= d^2 - 1, proving N <= d. Internal block terms compress to one-block terms and cannot alter this argument. QED.

Scope of the obstruction.

Theorem 6 applies to direct projection and first-order degenerate perturbation theory, or to an exact invariant product code whose effective Hamiltonian is the direct restriction. It does not exclude an encoding dressed across block boundaries or higher-order virtual processes.

In particular, I do not claim that the Schmidt rank at perturbative order k is always bounded by (d^2)^k. A resolvent such as (H_B + H_C - E)^{-1} has denominators depending jointly on many internal excitation energies. Its operator Schmidt rank need not be bounded by the number of bare boundary factors. Additional assumptions on the penalty spectrum would be needed for that stronger claim.

There is a related symmetry obstruction to the most direct protection scheme. A nontrivial microscopic on-site SU(N) action at fixed d induces a nonzero homomorphism su(N) -> pu(d). Simplicity of su(N) makes it injective, so N^2 - 1 <= d^2 - 1. Thus SU(N) with N > d cannot be protected by an ordinary nontrivial on-site action on the original sites. An encoded symmetry supported on long blocks, a dressed non-on-site symmetry, or an emergent symmetry is not excluded.

What is still required for an actual unbounded-c construction.

Starting from Corollary 5, an attractive target would be an effective antiferromagnetic ULS chain on N = 2^b colors. Its central charge would be N - 1 by the supplied background. Achieving this requires all the following ingredients in one fixed Hamiltonian h_b for each b.

Such a target would be consistent with the supplied entropic obstruction. Its conformal onset need not be the block length b. The ULS onset estimate in effective block units requires at least b(N - 1)/(3 log N) physical sites, of order N when N = 2^b. Also, the question permits choosing the effective exchange much smaller than every auxiliary gap for each fixed b. The shrinking gaps in the selector are not themselves a no-go theorem.

(i) A uniform-in-volume effective description. One would need, below the auxiliary gap, an effective interaction of the form

\[
H_{\rm eff}=J_b\sum_jP_{N,j,j+1}+\sum_X R_{b,X},\qquad J_b>0,
\]

with a local bound such as

\[
\sup_j\sum_{X\ni j}e^{\kappa_b\operatorname{diam}X}
\|R_{b,X}\|\leq\epsilon_b J_b,
\tag{25}
\]

where the parameters may depend on b but not on the number of blocks. A finite-size simulation theorem with a penalty or precision that must grow with L would not meet the question. Nor does a total operator-norm error O(L epsilon_b) resolve gaps of order J_b/L.

A naive global Schrieffer–Wolff argument is also insufficient. The retained data space has an extensive interacting bandwidth, while one auxiliary excitation costs only a fixed amount. The full retained band and excited bands can overlap as L grows. A local effective-Hamiltonian theorem in the relevant low-energy regime, an exact encoding, or another volume-uniform argument is needed. The block-selector proof does not supply it.

(ii) Protection of the infrared theory. Small local errors in (25) are not by themselves a proof of the same asymptotic CFT. A nonzero relevant perturbation can grow at arbitrarily large distance even if its coefficient is tiny. Marginally relevant perturbations can postpone the departure to a very large but finite scale. The limit L -> infinity is taken at fixed h_b, not with the simulation error sent to zero as L grows.

It would be too strong to assert that every SU(N)-breaking error is a relevant scalar at SU(N)_1. For clarity, the usual target-CFT identification assigns the level-one primary labeled by the p-th antisymmetric fundamental weight the total dimension

\[
\Delta_p=\frac{p(N-p)}N,
\]

and the effective one-cell translation phase exp(2pi i p/N). Here I use the standard affine-primary weight relation h_R = C_2(R)/(k + h^vee) as additional target-CFT input. The needed Casimir can be checked directly: with roots of length squared two, (Lambda_p, Lambda_p) = p(N-p)/N and (Lambda_p, rho) = p(N-p)/2, so C_2(Lambda_p) = (Lambda_p, Lambda_p + 2rho)/2 = p(N-p)(N+1)/(2N). At level one, h_p = C_2/(N+1), giving the stated Delta_p = 2h_p. The translation assignment is an additional standard identification for the fundamental ULS chain; it has not been derived for any proposed microscopic encoding here. With that translation action, the nontrivial primary perturbations are excluded by translation. Current perturbations and marginal current-current anisotropies still require analysis. Exact encoded SU(N), the correct translation action, and control of the marginal basin would be useful sufficient protection data. I have not established them for a fixed-d realization with N growing.

(iii) Commensurability. The selector is solved on L divisible by b. A ULS target can be studied on the further subsequence for which the number of blocks is divisible by N, avoiding a ground-state representation mismatch. Couplings can remain independent of L throughout such a sequence. On compatible rings, the b broken-translation vacua would give a direct sum of finitely many CFT sectors, not b times the central charge in a single sector.

This does not automatically settle other ring lengths. Incompatible rings require changes in cell lengths or defects. If added couplings make forced defects mobile, they can have nonrelativistic finite-size levels. If they remain immobile, large degeneracies also need treatment. Under a definition requiring the complete low-energy spectrum on every sufficiently large periodic ring to be purely conformal, merely constructing the commensurate subsequence is incomplete. A successful counterexample must specify and control these sectors.

(iv) Additional modes. For each fixed b the undisturbed selector has a positive gap above its clock and pattern ground spaces, so its finite clocks do not themselves force a long-wavelength history-state mode. It remains necessary to prove that the perturbations preserve the auxiliary gap and do not introduce mobile gapless defects. Conversely, a construction using a clock whose length grows with the full ring can leave quadratic history-state modes. Such modes cannot be relabeled as an extra unitary CFT. The supplied universality and undecidability facts do not establish the absence of these modes.

(v) The actual interaction. No fixed-dimensional nearest-neighbor matrix implementing the required higher-order exchanges and protection has been produced here. Finite real-parameter precision is not an obstruction imposed by the question; h_b may encode b to arbitrary precision. The unresolved problem is an L-independent local Hamiltonian with a controlled conformal infrared limit, not the storage of b in a coefficient.

Thus the explicit construction removes arbitrary block size and large block dimension as obstacles. The remaining obstacles are the effective exchange, its volume-uniform derivation, its infrared stability, and the complete periodic-sector spectrum. The direct-projection rank bound explains why the obvious first-order attempt stops.

Numerical checks.

The accompanying c7.py uses only numpy and scipy and reads no input files. It checks the exact local commutators of the marker and forbidden-bond projectors, the selected energy-density minimum for periods 5, 8, 16, 32, and 64 over a range of candidate lengths, the modular functional signs at c = 1, 2, 10, and 100, and the connected swap rank for N = 2, 3, and 4. These are consistency checks; the proofs above do not depend on sampling.

For t = 1, K = 8 and delta at the midpoint of (19), the period-five ring has ground energy -0.794593112989 with degeneracy five; the period-six ring has ground energy -0.551711241165 with degeneracy six. Both agree with (16). Their next gaps are 0.794593112989 and 0.551711241165, respectively, because the marker-free zero-energy states lie below the first internal clock excitation. An initial unshifted sparse eigensolver omitted these isolated zero eigenvalues. The final code shifts the matrix before Lanczos iteration and explicitly verifies both zero-energy product states. All final checks pass. The scripted checks took less than one second in the local environment.

Verdict.

For the full class in the question, the outcome is (c), a located obstruction to both proposed conclusions, together with two proved conditional bounds. The algebraic bound c_T <= d^2 - 1 requires current-algebra completeness; the finite-volume bound requires a controlled ring that already contains exponentially many states in c_T. Neither condition follows from fixed d and nearest-neighbor locality in this analysis.

The counterexample route advances to an explicit fixed-d block selector and large block ground manifolds. It stops before a conformal Hamiltonian. The next concrete target is a volume-uniform higher-order or dressed encoding of the SU(2^b) permutation interaction with proven symmetry and infrared stability. Until that step and the periodic-sector issue are resolved, large block spaces and Hamiltonian universality do not establish unbounded central charge.
