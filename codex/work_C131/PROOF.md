# C131. Compactness and the missing local bounds

2026-11-07. This note proves the topological reduction. It does not prove C(2) finite. Claims are tagged [Proved], [Numerical], [Assumed], [Heuristic], or [Open]. Literature results are distinguished from proposed consequences. The numerical record is in NUMERICS.md.

## 1. Space, normalization, and lengths

[Proved, C107] Let V be the real space of Hermitian two-qubit bonds and

    G = {alpha I + o tensor I - I tensor o}.
    M([h]) = inf_{g in G} ||h+g||.
    S = {[h] in V/G : M([h])=1}.

This is a norm on V/G and S is compact. The proof and existence of a minimizing representative are in ../work_C107/PROOF.md, Section 1. We use that result without a new compactness assumption. Scalar shifts and coboundaries cancel in E_(2L)-2E_L. Positive rescaling cancels in Q_L. Negative rescaling is not a gauge equivalence.

[Proved] V/G has dimension 12, so S has dimension 11. Simultaneous on-site SO(3) rotations preserve M and every ring spectrum. Their compact quotient T=S/SO(3) is also compact and has generic dimension 8. It is not justified to call T an eight-sphere. At symmetric bonds the stabilizer grows and the quotient has singular strata. One may prove the cover on S or on T and lift it. The numerical nine-coefficient normal form parametrizes T redundantly at degeneracies, not S before rotations.

We use S_a=sigma_a/2 and

    h = sum_a J_a S_a tensor S_a
        + sum_(a,b,c cyclic) D_a (S_b tensor S_c-S_c tensor S_b)
        + (1/2) sum_a b_a (S_a tensor I+I tensor S_a).

Thus the physical uniform field in H_N is +b.S. The field B of C129 is b_z=-B. At N=2 both oriented bonds are included, H_2=h_12+h_21. The DM term therefore cancels at N=2. No Jordan-Wigner sector restriction is imposed on E_N, the lowest eigenvalue of the full periodic spin Hamiltonian.

[Definition] The common-length target of this note is

    C_even(2) = sup_{[h] in S} sup_{L=2,4,6,...} Q_L(h),
    Q_L(h) = L(E_(2L)(h)-2E_L(h))/M(h).

All even rings are admitted in this upper-bound target. This is an explicit operational convention, not a claim that every even ring carries the continuum vacuum. It removes the familiar period-two frustration and makes both L and 2L even. We also distinguish the stronger envelope C_all(2), with every integer L>=2. C129 asks for C_U1,all on its subclass. The lemma below applies to each of these choices separately. No equivalence between finiteness of C_even and C_all is proved here.

[Assumed, only for identifying the physical Casimir coefficient] A vacuum-compatible sequence L_j for a given bond must make both L_j and 2L_j realize the untwisted conformal vacuum to order 1/L_j. Concretely, with that bond's physical A, require

    E_(L_j) = e_0 L_j-pi A/(6L_j)+o(1/L_j),
    E_(2L_j) = 2e_0 L_j-pi A/(12L_j)+o(1/L_j).

Then Q_(L_j) -> pi A/(4M). Evenness alone does not imply these expansions for arbitrary twists, charge filling, or an unclassified translation-breaking phase. A bound on all even rings bounds any such even vacuum subsequence. If only a different admissibility set is wanted, use the graph version of the lemma below, without silently replacing it by all even rings.

[Proved, check against C107] C107 uses a different observable, (6/pi)L(Le_0-E_L), not the doubling Q_L of this note. Its ferromagnetic XX example with hopping t has

    E_N=-2t csc(pi/N), N even;
    E_N=-t csc(pi/(2N)), N odd;   M=t.

Consequently the doubling observable here obeys

    Q_L=2L[2 csc(pi/L)-csc(pi/(2L))], L even,
    Q_L=0, L odd.

The even limit is pi/2=pi A/(4M), with A=2t. The odd value is exactly zero, since E_(2L)=2E_L for odd L. Thus all lengths can be allowed for an upper bound, but the raw full-length limit cannot be identified with A. Classical antiferromagnetic Ising similarly has Q_L=0 for even L and Q_L=-4L for odd L in Pauli normalization h=Z tensor Z, M=1. Odd forced defects do not by themselves threaten an upper bound.

## 2. Reduction lemma

[Proved] Fix once and for all an index set A subset {2,3,...}. For every fixed L, Q_L is continuous on S and descends continuously to T. The following are equivalent.

(i) C_A := sup_{h in S,L in A} Q_L(h) is finite.

(ii) For every h_0 in S there exist a relatively open neighborhood U of h_0 and a finite constant C_0 such that Q_L(h)<=C_0 for every h in U and every L in A.

(iii) S has a finite open cover U_1,...,U_m with finite constants C_1,...,C_m satisfying the same uniform estimate on each U_i.

Proof. Let delta=M([h]-[k]). Choose a minimizing representative of the difference, including its scalar shift, of operator norm delta. On an N-site ring its coboundary telescopes, and its scalar part shifts the energy by N times a constant. The variational principle gives

    |E_N(h)-E_N(k)-N alpha| <= N delta

for a suitable alpha (its sign depends on the representative). The alpha terms cancel between sizes 2L and L. Since M(h)=M(k)=1 on S,

    |Q_L(h)-Q_L(k)| <= L(2L delta+2L delta)=4L^2 delta.

This proves continuity. Rotation invariance proves descent. The factor L^2 is not bounded uniformly in L.

If (i) holds, take U=S and C_0=C_A. If (ii) holds, choose one such neighborhood for each h_0. Compactness gives a finite subcover U_1,...,U_m, so (iii) holds. If (iii) holds, then Q_L(h)<=max_i C_i everywhere, proving (i). No interchange of a limit with a supremum occurs. In fact the equivalence of local and global boundedness only uses compactness, not continuity. Continuity is what ensures that the problem cannot fail at a bounded set of lengths. QED.

[Proved, prescribed admissibility graph] More generally, fix a gauge- and rotation-invariant assignment A(h) of admissible lengths before taking the supremum. Set Gamma={(h,L): L in A(h)}. The same equivalence holds with all local and global suprema restricted to Gamma. No continuity of A(h) is needed. It is essential that the same Gamma occur on both sides. An ambient bound on all even lengths is sufficient, but not logically necessary, for a smaller physically selected Gamma. A sector projection changing with h must not be substituted for the continuous full-ring E_L without a separate continuity argument.

[Proved, finite-length patch] The elementary bound established in C98 is

    E_(2L)-2E_L <= 4M, hence Q_L<=4L.

Indeed, after centering an optimized bond to norm M, replacing the two closing bonds of two disjoint L rings by the two cross bonds changes the Hamiltonian by four terms, of total norm at most 4M. The product of the two ground states is a valid trial state. Thus a local lemma for L>=L_0 extends to every admissible L by taking max(C_0,4L_0). Alternatively, continuity and compactness bound the finitely many small lengths without any gluing constant.

[Open source discrepancy] specs/_Cd_target.md quotes 2M for product gluing, but ../work_C98/PROOF.md, Section 3, and RESULT.md explicitly prove 4M. I did not find a proof of the sharper general constant in the requested sources. Nothing here needs it. In particular, I do not promote the numerical Q_2 maximum of 4 to a theorem using that unverified sharper statement.

[Proved, localization of a possible failure] Let O be the set of points admitting (ii). O is open. The bad set B=S\O is closed, hence compact. If C_A is infinite, choose Q_(L_j)(h_j)>j. Fixed-L continuity and compactness imply L_j is unbounded, and we can pass to a subsequence with L_j->infinity and h_j->h_*. If h_* were in O, its neighborhood would eventually contain h_j and contradict its finite bound. Therefore h_* belongs to B. This locates a failure at some accumulation bond, but does not identify its physical phase.

## 3. Useful sufficient forms of a local lemma

[Proved] Write r_N(h)=E_N(h)-N e_0(h). On normalized bonds, any estimate

    |r_N(h)| <= K/N,  h in U, N>=N_0

gives

    Q_L <= L(|r_(2L)|+2|r_L|) <= 5K/2

for admissible L>=N_0. Uniformity in both the bond and N is essential. An expansion at each separate bond does not provide it.

[Proved] The stronger gapped estimate

    |r_N(h)| <= K exp(-N/xi),  h in U, N>=N_0

gives

    Q_L <= K L[exp(-2L/xi)+2exp(-L/xi)]
         <= (5K xi)/(2e).

The last inequality bounds the two terms separately using sup_(x>=0) x exp(-ax)=1/(ae). Add max(4N_0,5K xi/(2e)) for the small lengths. Constants K and xi may depend on the center h_0; compactness does not require one preliminary radius or one preliminary xi for all centers.

[Proved] Strictly exponential convergence is more than necessary. A uniform remainder bound |r_N|<=K f(N) suffices whenever sup_N N f(N)<infinity. For example f(N)=exp(-a N^theta), a,theta>0, works with Q_L<=3K(a theta)^(-1/theta) exp(-1/theta). A uniform O(1/N) bound also works. This matters because spectral-flow locality estimates need not be strictly exponential.

[Proved, phase envelopes] If finitely many phase branches on U satisfy

    E_N(h)=min_a [N e_a(h)+r_(a,N)(h)],
    |r_(a,N)(h)|<=K exp(-N/xi),

uniformly in a,h,N, then e_0(h)=min_a e_a(h) and the same bound holds for |r_N|. To see this, every branch is at least N min_a e_a-K exp(-N/xi), while a minimizing bulk branch gives the reverse inequality. This formulation tolerates crossings and varying ground-state multiplicities. It does not require differentiability of e_0.

## 4. Strata and their local obligations

[Definition] In this table B(h_0;r,C) means the ambient statement Q_L(h)<=C for every normalized bond within quotient distance r of h_0 and every even L>=2. Distance on T minimizes M over simultaneous rotations. Relative neighborhoods within a symmetry class do not count as ambient neighborhoods.

[Open classification] The labels below overlap. Assigning priority to special loci can make a disjoint set partition, but no theorem in TOOLBOX classifies the entire gapless locus of the general bond. The last row is necessary for a logically exhaustive reduction. Without it, this is a list of proposed local lemmas, not a proved finite stratification by physical universality classes.

| Class or locus | Precise statement needed | Relevant tool and its actual scope | Status |
| :--- | :--- | :--- | :--- |
| (a) U(1), Jx=Jy, Dx=Dy=bx=by=0 in some frame | First sup_(h in U1,L>=2) Q_L(h)<infinity. Then B(h_0;r,C) for each center in this class, allowing all transverse perturbations. | C129, Yang-Yang sector ground states, Kozlowski 1508.05741 Theorems 3.4, 3.5, 4.1, twist variational bounds. These do not give transverse control. | [Open] Handed to C129 for the relative class. No result from the parallel run is presumed. Transverse ambient lemma remains separate. |
| (b0) Quasi-free bonds | sup_(h in F,L even) Q_L(h)<=C_F, or local uniform versions covering F. F means Jz=0 and Dx=Dy=bx=by=0 in a JW frame, up to rotations. Jz=0 alone is insufficient. | Exact fermion sums including both JW parities; C118/C119/C128 for recorded subfamilies; C130 free leg for the full family. | [Proved] Recorded free subfamilies only. [Open] Whole free finite-size theorem as specified in C130; no ambient conclusion. |
| (b1) Nondegenerate free critical centers and nearby interacting bonds | B(h_0;r,C), including detuning across the critical surface and both spin parity sectors, with uniform control through mL of order one. | Giuliani-Mastropietro 1211.3030 Theorem 1.1 gives a classical Ising cylinder coefficient; Benfatto-Mastropietro 0907.2837 gives quantum-chain universality relations. C130 proposes the quantum-ring adaptation for parity-even perturbations. | [Open] C130 is drafted, not launched. Even its proposed result is relative to fermion parity. Generic bx,by,Dx,Dy perturbations are outside it. |
| (c1) A stable single gapped phase | G-FS and G-ST below, or directly a uniform O(1/N) remainder on an ambient U. | Hastings-Koma Theorem 2.8 proves clustering within a fixed system. Quasi-adiabatic continuation compares a specified uniformly gapped path. Neither statement as read proves G-FS from a ring gap alone. | [Open] General frustrated-chain local lemma. The exact missing inputs are separated below. |
| (c2) Classical Peierls centers, and ordered coexistence near them | Uniform phase-envelope remainders from Section 3, on lengths compatible with the reference phases, plus control of any additional even lengths if claiming C_even. | Borgs-Kotecky 1990 Theorem 4.1 is a classical cubic-torus estimate. Borgs-Kotecky-Ueltschi 1996 Theorem 2.1 assumes spatial dimension >=2. A one-dimensional quantum-cylinder adaptation uniform as beta->infinity is still needed. | [Proved] Classical C<=0 from the inherited classical theorem. [Assumed] The quantum-ring corollary advertised in the survey is not verified by the stated classical theorem alone. |
| (c3) Frustration-free gapped centers satisfying Local-TQO and Local-Gap | Stable low-energy gap plus a uniform ring-energy comparison, or a direct B(h_0;r,C). | Michalakis-Zwolak 1109.1588 Theorem 1 proves the gap part under both named local conditions. Ring-frustration-free h itself has E_N=N e_0. | [Proved] Pointwise zero defect and the cited stability theorem under its hypotheses. [Open] Requested perturbed energy estimate unless a further comparison theorem is supplied. |
| (d) Other non-integrable critical surfaces, away from controlled free neighborhoods | B(h_0;r,C) for every such center, including both sides of the surface. | No general finite-size theorem found in TOOLBOX. Ising universality is a physical expectation, not the required uniform remainder. | [Open] No analytical attempt in this run. Coverage fraction by C130 balls is not presently defined, see Section 7. |
| (e1) Heisenberg and its symmetry/commensurate-twist images | B(h_H;r,C) with marginal flow and arbitrary allowed transverse perturbations included. | Bethe ansatz at the integrable center does not give this ambient statement. SU(2)_1 with a marginally irrelevant interaction describes the expected IR. | [Open] No analytical attempt. [Numerical] Center values decrease, 4, 2.79125273, 2.61792430, 2.55824732. |
| (e2) Ferromagnetic endpoint Delta=-1 | B(h_F;r,C) uniformly through the z=2 crossover. | C129 addresses its U(1) directions; simple-Fermi-point RG excludes the degenerate center. | [Open] No analytical attempt. |
| (e3) Free Lifshitz/saturation junctions | B(h_lif;r,C), including interaction, pairing, and field detunings through the quadratic touching. | Exact free sums can treat the point; C130's nonzero-velocity hypothesis fails. | [Open] No analytical attempt. |
| (e4) Other boundaries with diverging crossover length | B(h_0;r,C) uniform in length and detuning. | No theorem from a gap or a pointwise CFT expansion alone. | [Open] No analytical attempt. |
| Residual R | Either prove R is empty, or B(h_0;r,C) for every h_0 in R. Here R is the complement of the chosen classified loci. | A rigorous phase-classification/exhaustion theorem would remove this row. None is supplied by TOOLBOX Section 8.4. | [Open] Required for an unconditional exhaustive cover. |

[Proved conditional assembly] If the ambient local statement in every row holds and R is empty or satisfies its stated local lemma, Section 2 supplies a finite subcover and C_even(2)<infinity. A finite number of types of lemma is not yet a finite cover by neighborhoods. The latter exists only after the local statements have been proved at every point. The argument is qualitative unless the neighborhoods, constants, and a covering certificate are made quantitative.

## 5. The gapped problem stated without an implicit theorem

[Open, G-FS] Suppose a single translation-invariant bond h_0 has, for every admissible ring N>=N_0, an exactly q-dimensional lowest eigenspace and a spectral gap at least gamma>0 above it. Prove an estimate

    |E_N(h_0)-N e_0(h_0)| <= K(h_0) exp(-N/xi(h_0)).

For ordered systems with tunneling, replace the exactly degenerate ground space by a specified q-dimensional low-energy band, with width tending to zero and a fixed positive gap to its complement. Those are different hypotheses. A common gapped ordered chain has a split doublet, and fails the literal exact-ground-space premise with q=1.

[Open, G-ST, restricted to stable phases] Establish a neighborhood U in which the appropriate spectral separation persists uniformly in N and h. To combine it with G-FS, the latter's constants must be uniform under the resulting common gap, interaction-norm, phase, and locality hypotheses. Pointwise finiteness of K(h) and xi(h), even combined with pointwise gap stability, is not enough. A theorem bounding the constants in terms of common controlled data, or a direct uniform proof, is required.

[Proved limitation on the formulation] A single uniform positive gap above the exact ground space cannot persist on an unrestricted neighborhood of every degenerate gapped point. For the classical ferromagnet

    H_N(b)=-sum_j Z_j Z_(j+1)-b sum_j Z_j,

b=0 has two ground states and gap 4 (N>=3). For b>0 the all-up state is the unique ground state and the all-down state's energy lies 2bN above it. Thus the first gap is at most 2bN. Taking b arbitrarily small at any fixed N>=N_0 defeats any common positive gap over the neighborhood. Rescaling these bonds to M=1 is continuous near b=0 and does not remove this failure. Nevertheless

    E_N(b)=-N(1+|b|),  Q_L(b)=0

for every b and L. This is not a counterexample to C(2) or to the energy conclusion. It shows that the proposed two-input route must use phase branches, restrict to a stable phase, or formulate a more appropriate band condition at coexistence. Keeping a fixed number q of lowest eigenvalues separated uniformly also needs care when an extensive branch splitting crosses droplet excitations.

[Proved literature theorem, scope checked] Hastings-Koma, [math-ph/0507008](https://arxiv.org/pdf/math-ph/0507008), Theorem 2.8, bounds the connected correlation after subtracting the ground-sector projection term, exponentially in the separation for finite-range interactions under a uniform spectral gap. For a unique ground state this is ordinary connected clustering. For a degenerate space it is not a statement that every ground-state mixture has decaying ordinary connected correlations. In either case it compares observables in a system, not ground energies on two different ring sizes.

[Assumed missing implication] To obtain G-FS from locality one would need, for example, a comparison between the translation-averaged local bond expectation on the N ring and its infinite-volume value, with error K exp(-aN), uniformly in h. Multiplying by N would give |r_N|<=KN exp(-aN), hence an exponential bound with a smaller exponent. Clustering alone does not establish this comparison. Cutting the ring and appealing to a gapped spectral flow also assumes a gapped interpolation that need not exist, especially when the cut produces boundary modes. Quasi-adiabatic continuation assumes the spectral separation along the path, rather than proving it for a generic perturbation.

[Proved literature theorem, scope checked] Michalakis-Zwolak, [arXiv:1109.1588](https://arxiv.org/pdf/1109.1588), Theorem 1, assumes frustration-freeness, a uniform gap, Local-TQO and Local-Gap. It proves stability for sufficiently weak quasi-local perturbations. This is a source for G-ST on its class, not a theorem for arbitrary frustrated chains. Moon-Nachtergaele, [arXiv:1804.04750](https://arxiv.org/abs/1804.04750), concerns frustration-free one-dimensional systems with open boundaries and local-order assumptions. The boundary setting matters for a ring-energy application.

[Assessment] No general frustrated-chain theorem implying both G-FS and G-ST with the required hypotheses was verified. The next tractable statement is a uniform ring-energy comparison within a stable, unique-ground-state phase with a known gapped path. At classical coexistence, the branch-envelope formulation is the appropriate alternative.

## 6. What the contour theorem actually supplies

[Proved literature theorem, checked] Borgs-Kotecky, J. Stat. Phys. 61 (1990) 79, [preprint, Theorem 4.1(i)](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/11/A-Rigorous-Theory-of-Finite-Size-Scaling-at-First-Order-Phase-Transitions.pdf), proves for its classical contour model on a periodic cube T of volume L^d

    |Z_per(T)-sum_q exp(-f'_q L^d)|
         <= exp(-f L^d-b_0 tau L),
    f=min_q f'_q,

under the Peierls/contour assumptions and sufficiently large tau. This is the theorem quoted in TOOLBOX. It is not, as stated, a spectral-gap theorem for a quantum chain or a bound on a cylinder of fixed circumference and unbounded imaginary time.

[Proved literature scope restriction] Borgs-Kotecky-Ueltschi, CMP 181 (1996) 409, [Theorem 2.1, preprint page 10](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/11/Low-Temperature-Phase-Diagrams-for-Quantum-Perturbations-of-Classical-Spin-Systems-.pdf), explicitly assumes spatial dimension d>=2. Its d is spatial dimension, not on-site Hilbert-space dimension. Theorem 2.2 inherits that assumption. Thus these stated theorems do not apply directly to the spin chain here. This is an additional obstruction to treating the survey claim as a corollary. Datta-Fernandez-Frohlich, J. Stat. Phys. 84 (1996) 455, [DOI](https://doi.org/10.1007/BF02179651), remains a recalled candidate reference; its fixed-circumference zero-temperature theorem was not verified. The claim that BK Theorem 4.1 alone supplies both the quantum gap and ring-energy conclusion should not be used as a proved input.

[Proved sufficient replacement] A usable quantum-cylinder estimate would be

    |log Z_N(beta,h)-log sum_a exp[-beta N e_a(h)]|
          <= beta K exp(-N/xi)+R_N(beta,h),
    lim_(beta->infinity) R_N(beta,h)/beta=0,

with K,xi uniform on U and over compatible N. Divide by beta and take beta->infinity at fixed N. The left terms tend to -E_N and -N min_a e_a, respectively, giving the desired ring-energy estimate. This is the concrete contour estimate to request next.

[Proved diagnostic] A formula Z_N=sum_a exp(-beta N e_a)[1+O(exp(-aN))+O(exp(-b beta))] with a relative error bounded independently of beta at fixed N would imply E_N=N min_a e_a exactly after beta->infinity. Generic weakly quantum chains have nonzero exponentially small ring corrections. Thus the survey's displayed relative-error formula cannot be used in that uniform sense. Wrapping contributions must be allowed to accumulate linearly in beta in log Z. A cubic torus estimate does not settle that limit.

## 7. Handoffs and the requested coverage estimate

[Open, C129 handoff] Prove the all-length U(1) estimate with its twist and charge minimization. Its finite-size bound is stronger than the even-length restriction on this subspace. It does not assert that every length has the conformal vacuum expansion. No waiting for or reading a result of the parallel run was needed here.

[Open, C130 handoff] Three distinctions must survive the adaptation. First, parity-even perturbations occupy a proper subset of the full bond space. Transverse fields and off-axis DM are generally nonlocal under the chosen JW map. In fact b cross D != 0 is a rotation-invariant obstruction to the standard on-site parity frame, since parity requires both vectors to lie along its axis. Such perturbations occur arbitrarily close to the parity-preserving class. Thus allowing a change of spin axes does not turn C130's proposed relative neighborhood into an ambient one. Second, a Casimir expansion on a tuned critical hypersurface does not by itself control a neighborhood across the mass crossover. Third, the free finite-size bound should be in terms of M or controlled band norms. The draft's bound solely by sum of velocities at Fermi points cannot hold for every quasi-free bond: a gapped free chain has no Fermi points but generally has nonzero finite-size defect. The requested full class bound can still be true with M on the right.

[Proved, free consistency check] For h=Sx tensor Sx+(Sz tensor I+I tensor Sz)/2, the chain is gapped (the usual transverse-field Ising chain with field above the Ising critical value). Its quasiparticle energies obey epsilon(k)^2=5/4+cos(k), up to a momentum shift, so min epsilon=1/2. Direct two- and four-site diagonalization, equivalently the parity-resolved free sums, gives

    E_2=-sqrt(5)/2,
    E_4=-sqrt((5+sqrt(17))/2),
    E_4-2E_2=sqrt(5)-sqrt((5+sqrt(17))/2)>0.

For the antiperiodic sector, N=2 samples k=pi/2,3pi/2, giving E_2=-sqrt(5)/2. At N=4 the two distinct values of epsilon are sqrt(5/4+1/sqrt(2)) and sqrt(5/4-1/sqrt(2)), each twice. Their sum has square (5+sqrt(17))/2, giving E_4. The Bogoliubov vacuum in this sector has even fermion parity: starting at zero Sx coupling it is the polarized state on an even ring, and no single-particle energy crosses zero as that coupling is increased to one. The unconstrained periodic-sector energies are -1 and -(1+sqrt(5)/2), respectively, both higher than the displayed antiperiodic energies; imposing its parity constraint cannot lower them. These checks identify the full spin ground energies, not just an unprojected fermion sum. The last displayed positive-defect inequality follows from sqrt(17)<5. There are no Fermi points. Thus a uniform free bound with only the sum of Fermi-point velocities on the right really fails already at L=2. The numerical spin matrices reproduce these energies to floating-point accuracy. This objection concerns the draft's stated free-leg normalization, not its proposed M-normalized target.

[Proved literature theorem, scope checked] Giuliani-Mastropietro, [arXiv:1211.3030](https://arxiv.org/pdf/1211.3030), Theorem 1.1 and (1.4), proves c=1/2 in the finite-size free energy of its weakly perturbed classical Ising model, taking the long-direction limit before the circumference limit. It is not the general quantum-ring neighborhood theorem. Benfatto-Mastropietro, [arXiv:0907.2837](https://arxiv.org/abs/0907.2837), concerns universality relations in weakly interacting spin chains. The stronger quantum-ring claim remains C130's proposed work.

[Open, numerical coverage] There is no supplied numerical value of lambda_0(h_free), no established ambient C130 ball, and no specified sampling measure on the unclassified critical surfaces. A numerical fraction lying inside such balls is therefore undefined, not zero. No percentage is reported. A distance-to-free statistic with an arbitrarily chosen threshold would not estimate rigorous coverage. Once radii and a sampled critical surface are available, the test is whether min_(free centers) dist(h,h_free)/lambda_0(h_free)<1 in the stated quotient norm, with the perturbation's symmetry restrictions also checked. Small-ring optimization is not a phase-surface map and cannot supply those missing data.

[Assessment, order of attack] Complete C129 on its relative class and C130's finite-size, crossover-uniform parity-preserving statement. In parallel as a future project, formulate the gapped ring-energy comparison with known stability hypotheses, and a phase-envelope version near classical coexistence. The remaining ambient transverse directions, multicritical neighborhoods, and an exhaustion of the residual gapless set must be retained as explicit obligations. Neither continuity nor compactness eliminates any of them.
