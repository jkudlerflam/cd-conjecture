# C52. Exact-cell redesign and its obstruction

[proved, outcome] No fixed-local-dimension nearest-neighbor redesign closing S was obtained. There is a finite-dimensional obstruction to the most direct proposed redesign. An exactly invariant product-cell code with separated physical endpoints cannot have C39's full-rank nearest-cell fermion hopping, even with additional number-conserving interactions. The unchanged C44 sine-head wrapper independently fails exact invariance under nonzero marker motion. Neither statement refutes the original perturbative construction or all possible redesigns.

[proved, provenance and scope] The separated-endpoint theorem is already proved in [C8 RESULT Section 10](../work_C8/RESULT.md). Section 2 below states its hypotheses, gives its finite proof for completeness, and attributes it to C8. The new application uses the one-particle block of the commutator for C39's general complex hopping, extends the obstruction to the direct-sum mobile code, and supplies an explicit C44 head-leakage formula and a boundary-local repair obstruction. The independent rank and determinant bounds clarify which other exact-fermion designs could evade the first theorem. These are finite algebra results. No constructive RG theorem is invoked.

[proved, claim convention] A tag covers the following paragraphs and formulas until the next tag or heading. “Imported” identifies a result supplied by the cited project notes. “Assumed” identifies an unproved input or an unrealized construction. Numerical statements are separate from the proofs.

## 1. The exact object that would need to be constructed

[proved, definitions] Fix n before taking the physical ring length to infinity. A physical site has dimension d independent of n. For a specified tiling, a cell of length b has an isometry

    G_b : F_n -> (C^d)^(tensor b),       F_n = Fock(C^n),
    P_b = G_b G_b^dagger.

A frozen product encoding has G = tensor_j G_b and P = G G^dagger. In the C44 mobile version there is a sum over marker patterns tau, with cell lengths b or b+1,

    P_mobile = direct_sum_tau P_tau,
    P_tau = tensor_(cells j of tau) P_(length j).

Different marker patterns are orthogonal. C44 has b=n+3 and d=4464. A product encoding may contain arbitrarily entangled states inside each cell. The hypothesis excludes entanglement of the encoding itself across different cells.

The required exact invariance is Q H P=0, equivalently [H,P]=0 for Hermitian H. The compression P H P is always defined, including when this equality fails. Its existence says nothing about exact invariance.

[proved, imported target] C39 Sections 2-3 and C42 give a reference logical hopping matrix

    T_0 = -W [Gamma(w) + i eta Gamma(wv)],                (1)

with W>0, eta>0 and both Gram matrices positive definite. The field is an onsite matrix and does not change (1). The reference has n orbitals and 2n distinct simple Fermi roots, with positive root, quartic and speed margins at every fixed n. The moments are specified by a terminating finite-parameter procedure. Those settled results are used here, not recomputed.

The prospective exact fermion Hamiltonian would have a bond of the form

    K_(j,j+1) = sum_(a,b) [T_ab f_(j,a)^dagger f_(j+1,b) + h.c.]
                  + number-conserving interaction terms,

and arbitrary onsite terms. Fermions on different cells are represented in the usual ordered tensor product, including the parity factor on the left cell in an intercell hopping operator. All arguments below use actual matrix elements and retain this factor.

[proved, convention] A genuine vacuum-normal-ordered monomial with at least two annihilation operators vanishes on the global zero- and one-particle sectors. Thus a number-conserving quartic interaction cannot change the one-particle hopping matrix. If an interaction was written with subtractions or contractions, collect its quadratic part into T and its onsite part into A first. No argument depends on naming an actual quadratic term an “interaction.”

## 2. The inherited exact-product theorem

[proved, theorem, C8 Section 10] Partition a physical ring into M>=4 cells. Write

    H = sum_j B_j + sum_j v_(j,j+1).

B_j acts inside cell j. Each v acts only on boundary regions of its two cells. Within any cell the left and right boundary regions are disjoint. Let P=tensor_j P_j be any product code, and assume [H,P]=0. Define the connected compression K_(j,j+1) by subtracting the normalized partial traces and adding back the scalar trace. Then

    [K_(j-1,j), K_(j,j+1)] = 0.                         (2)

No gap, frustration freeness, translation invariance or termwise code preservation is assumed.

[proved, proof recalled with attribution] Put Q_j=I-P_j and apply Q_j on the output of cell j and P_i on all other outputs to H G. Exact invariance gives zero. Terms not incident on j vanish. After identifying the neighboring code spaces, the remaining equation has the form

    Q_j B_j G_j tensor I
      + sum_alpha Q_j A_alpha^R G_j tensor C_alpha^(j+1)
      + sum_beta D_beta^(j-1) tensor Q_j A_beta^L G_j = 0. (3)

Choose a Hermitian linear functional phi on cell j+1 with phi(I)=0. Apply phi on that cell and normalized traces on the other cells except j. This kills the internal and left-bond terms. It gives

    Q_j A_phi^R G_j = 0,

where A_phi^R is a Hermitian physical operator supported on the right boundary region. Hermiticity gives the reverse off-code block zero too. Therefore [A_phi^R,P_j]=0. Its compression, after subtracting a scalar multiple of I, is the phi-contraction of the connected right bond. These contractions span its active Hermitian Schmidt space. Apply the same argument to the left neighbor to obtain reducing physical lifts of the left Schmidt space.

The two families of physical lifts commute because their supports are disjoint. Restriction to a reducing subspace is multiplicative. Their logical restrictions therefore commute. Expanding the two connected bonds in their Schmidt spaces proves (2).

This proof is why internal counterterms or cancellation between different physical terms do not evade the theorem. Equation (3) used invariance only for the full sum.

## 3. Application to full-rank hopping, with interactions

[proved, theorem] Under the hypotheses of Section 2, suppose the compressed Hamiltonian is number conserving in logical fermions. Let T_j be the hopping block of K_(j,j+1) between one-particle states on its two cells. Then

    T_(j-1) T_j = 0.                                    (4)

In particular a repeated invertible T is impossible. This remains true with arbitrary number-conserving two-cell interactions that vanish on the vacuum and one-particle sectors, including arbitrary vacuum-normal-ordered quartics.

[proved, proof] For three consecutive cells, write |1,a>, |2,a>, |3,a> for the state with one particle of flavor a in the indicated cell and vacuum elsewhere. On this sector the two connected bonds have the block form

    K_12 = [[A_1,T_1,0],[T_1^dagger,A_2,0],[0,0,c I]],
    K_23 = [[c' I,0,0],[0,B_2,T_2],[0,T_2^dagger,B_3]],

where scalar and onsite blocks may be arbitrary. Consequently

    <1,a|[K_12,K_23]|3,b> = (T_1 T_2)_ab.               (5)

For the forward product the only path from cell 3 to cell 1 goes through cell 2. The reverse product cannot make that path. The normalized partial-trace subtractions defining connected bonds change onsite and scalar blocks only, so they do not change (5). Number conservation ensures that the one-particle subspace is invariant throughout this calculation. Higher normal-ordered interactions vanish there. Equation (2) now proves (4).

For any invertible T,

    ||[K_12,K_23]|| >= ||T^2|| >= sigma_min(T)^2 > 0,     (6)

which contradicts exact product-code invariance. The second inequality follows from ||Tx||>=sigma_min(T)||x|| twice, then restricting to a unit vector. In fact ||T^2|| can be much larger than this lower bound.

[proved, C39 hypotheses checked] For a complex vector x, positivity of the first Gram matrix gives

    Re x^dagger(-T_0/W)x = x^dagger Gamma(w)x
                         >= gamma_n ||x||^2,
    gamma_n = lambda_min(Gamma(w)) > 0.

Cauchy-Schwarz implies

    ||T_0 x|| >= W gamma_n ||x||,
    sigma_min(T_0) >= W gamma_n.                         (7)

Thus T_0 is invertible at every fixed n. The complex bus phase is fully included. No assumption of simultaneous diagonalization is needed.

If a proposed exact logical bond differs from the reference by an operator error of norm delta, its one-particle hopping changes by at most delta. For delta<W gamma_n/2, the actual hopping remains invertible with smallest singular value at least W gamma_n/2. C42 allows its finite-parameter error tolerances to be made arbitrarily small at each n. One may add this strictly positive tolerance to its existing momentum and speed tolerances. This is an additional requirement, not an assertion that every recorded diagnostic parameter already meets it.

Therefore replacing the C39 parent by commuting projectors, or by another frustration-free parent, cannot produce an exactly invariant product-cell realization of the same reference with sufficiently small corrections. The theorem is independent of the choice and spectrum of the parent. Enlarging d does not remove the separated-endpoint obstruction.

[proved, qualification] Equation (4) does not exclude all hopping in product codes. For example, different left and right orbitals can give a nonzero nilpotent intercell T, while intracell dynamics propagates between them. It excludes the invertible n-by-n hopping required here. Nor does it bound the central charge of an arbitrary interacting spin chain.

## 4. The mobile direct sum does not remove this obstruction

[proved, theorem] Suppose P_mobile from Section 1 is exactly invariant under a local Hamiltonian. Fix an admissible marker pattern tau and let Pi_tau project onto the full physical space with that marker pattern, retaining all internal excitations. Its diagonal Hamiltonian block is

    H_tau = Pi_tau H Pi_tau.

Then the product code P_tau is invariant under H_tau.

[proved, proof] For a code vector psi in P_tau, invariance gives H psi in the sum of the P_sigma. Projecting its output to marker pattern tau gives H_tau psi in P_tau. Hence (Pi_tau-P_tau)H_tau P_tau=0. Hermiticity gives the reverse block. Marker-changing terms cannot cancel leakage into the original marker pattern because their outputs are orthogonal to it.

In H_tau, physical terms are still internal-cell terms or boundary terms. Apply Section 2 whenever the two boundary regions in each cell are disjoint. Thus a common C39 payload bond in the diagonal length blocks is impossible in an exactly invariant direct-sum product code, by Section 3. A weak length-dependent perturbation does not help when its effect on each hopping matrix is smaller than the bound in (7). Products of neighboring invertible matrices remain invertible even when the two length labels differ.

[proved, scope] This statement allows arbitrary off-diagonal marker dynamics. It assumes that the diagonal length block retains C44's common payload hopping, or a sufficiently close version. It does not cover an entirely different encoded fermion model whose apparent payload propagation arises solely from marker-changing processes and which has no such diagonal hopping. No such alternative with the C39/C42 margins has been constructed here.

## 5. The literal C44 marker move has nonzero leakage

[proved, imported code] C44's selector factor is independent of its charge, helix and payload factors. A length-l cell has normalized head ground state

    phi_l(i) = sqrt(2/l) sin(pi i/l),   i=1,...,l-1.

The directed move #L -> R# grows a length-b cell to b+1 and shrinks the adjacent length-(b+1) cell to b. It leaves the left head position i unchanged and removes the first possible head position in the right cell. C44's formula is the compressed overlap, not the norm of the complete output.

[proved, theorem] Define

    A_b = sum_(i=1)^(b-1) phi_b(i) phi_(b+1)(i),
    a_b^2 = phi_(b+1)(1)^2
          = 2 sin(pi/(b+1))^2/(b+1),
    s_b = 1-a_b^2.

For unit directed selector hopping on the initial product head ground state, the output vector w_b obeys

    ||w_b||^2 = s_b,
    <phi_(b+1) tensor phi_b | w_b> = A_b^2,
    ||(I-P_head,target)w_b||^2 = s_b-A_b^4
                              >= s_b(1-s_b)>0.           (8)

This holds in particular for every b=n+3 in C44.

[proved, proof] In the target head spaces define

    u_i = phi_b(i) for i<b,       u_b=0,
    v_j = phi_(b+1)(j+1),         j=1,...,b-1.

The output is u tensor v. Its squared norm is s_b because ||u||=1 and ||v||^2=1-phi_(b+1)(1)^2. The left overlap is A_b. Reflection of both sine waves makes the right overlap A_b too. Orthogonal projection gives the first equality for the leakage.

Cauchy-Schwarz applied to the first b-1 components of phi_(b+1) gives A_b^2<=1-phi_(b+1)(b)^2=s_b. Hence A_b^4<=s_b^2, proving the lower bound. Since 0<a_b^2<1, it is strictly positive. Geometrically, the output has no amplitude at the newly available left head position i=b, whereas the new ground state has positive amplitude there.

[proved, inclusion of the other C44 tracks] On an initial product of complete selected-cell code states, the directed marker move produces w_b tensor xi in the target marker pattern, where xi contains the typed-hole transport and the local helix rotations. It is nonzero. For example C44's strictly positive beta_v implies that its projection onto the target wrapper ground factors is already nonzero. The final code is contained in the target head ground subspace. Consequently

    ||Q_mobile T_marker psi||^2
       >= (s_b-A_b^4) ||xi||^2 > 0                     (9)

for the output component in this target pattern. Here T_marker denotes this directed local move, with its nonzero coefficient removed. In a ring with sufficiently many cells, shifting a specified marker gives a pattern distinct from that reached by any other single marker shift. The other terms in C44 either preserve the marker pattern or shift a different marker. They cannot cancel this component. Thus the unchanged C44 Hamiltonian has Q_mobile H P_mobile != 0 whenever u_mov is nonzero, even if its payload ports are switched off.

[proved, a stronger local-repair statement] Keep C44's sine-head code but replace the marker transition by any operator supported in a boundary region. Suppose the growing cell has two untouched interior head positions i_1<i_2, and the selector letters in the operator's support are the same for both positions. The literal endpoint move has this property. On those two input components the operation acts identically on every other factor. A nonzero exact map into the new code would therefore require

    phi_(b+1)(i_1)/phi_b(i_1)
       = phi_(b+1)(i_2)/phi_b(i_2).                     (10)

These ratios are strictly increasing in i, so (10) is impossible.

To prove monotonicity, omit the constant normalization factor and set alpha=pi/(b+1), beta=pi/b. For 0<x<b,

    d/dx log[sin(alpha x)/sin(beta x)]
       = alpha cot(alpha x)-beta cot(beta x)>0.

The function y cot y decreases on (0,pi), since its derivative is (sin y cos y-y)/sin(y)^2<0. With alpha x<beta x this proves the strict inequality. All integer head positions used above lie in this interval.

This argument includes a sum of local counterterms reaching the same target marker pattern, provided their combined support leaves those two head positions untouched. For a fixed-size boundary region it applies once b is large enough. To obtain exact nonzero transport one must change the head ground-state family, act over a growing fraction of the cell, or abandon the product-head description. A payload-only redesign cannot do it.

## 6. Rank bounds without exact invariance

[proved, proposition] Let a physical interaction across a cell cut act on one d-dimensional site on either side. For arbitrary product-cell isometries, its compressed bond has operator Schmidt rank at most d^2. If its logical compression is number conserving and has one-particle hopping matrix T, then

    2 rank(T) <= d^2.                                  (11)

This bound needs no exact invariance or separated-endpoint hypothesis.

[proved, proof] Expand the physical bond in a basis of the d^2 operators on its left site,

    h = sum_(alpha=1)^(d^2) A_alpha tensor B_alpha.

Compression is a tensor product of linear maps on operator spaces. It cannot increase the Schmidt rank. On each logical side apply the linear map taking an operator X to its 2n matrix elements

    (<a|X|0>, <0|X|a>)_(a=1,...,n).

The resulting 2n-by-2n coefficient matrix of the compressed bond is

    [[0,T],[conjugate(T),0]],

up to interchange of row or column conventions. Pair creation and annihilation entries vanish by total number conservation. Its rank is 2 rank(T). Applying linear maps cannot increase a coefficient matrix's rank, which proves (11). Onsite and scalar terms vanish under both extraction maps. Genuine normal-ordered interactions of degree at least four also vanish in this extraction.

[proved, charge-refined bound] If the physical bond conserves an additive integer charge and the code identifies it with logical particle number up to a constant, let d_q be the dimensions of the charge-q eigenspaces of one boundary site. A charge-raising physical operator belongs to a space of dimension

    m_1 = sum_q d_(q+1) d_q.

Only those operators contribute to <a|G^dagger X G|0>. Decomposing the charge-preserving bond into transfer channels therefore gives

    rank(T) <= m_1 <= floor(d^2/4).                     (12)

For the last inequality group even and odd q. The sum includes some of the products between the two groups, so it is at most d_even d_odd, which is at most floor(d^2/4). This argument uses an actual integer-charge representation, as in C39/C44, not parity alone.

[proved, interpretation] Even if endpoint separation were removed, these bounds prohibit a full-rank n-channel family for all n at fixed d in this single-bond product compression class. At d=4464, (12) gives n<=4,981,824. This weak rank bound by itself says nothing against a finite-n violation of c<=d-1. The stronger separated-endpoint theorem is the obstruction relevant to the actual cells.

## 7. A deficient-rank free hopping cannot keep the same root count

[proved, theorem] Consider an n-orbital number-conserving nearest-neighbor Bloch Hamiltonian

    h(z)=A+zT+z^(-1)T^dagger,       |z|=1,

with rank(T)=r and det h not identically zero. Then det h is a Laurent polynomial with exponents between -r and r. There are at most 2r distinct Fermi roots, counted with multiplicity, on the unit circle. Consequently 2n distinct simple roots require r=n.

[proved, proof] Expand the determinant multilinearly in its columns. A term using more than r columns of T vanishes by linear dependence, regardless of the other columns. The same holds for T^dagger. A nonzero term containing p columns from T and q from T^dagger has exponent p-q, with p,q<=r. Hence all exponents lie in [-r,r]. Multiplication by z^r gives a polynomial of degree at most 2r. A nonzero polynomial has at most that many roots. At a simple Fermi crossing with the other eigenvalues nonzero, det h has a simple zero. The identically-zero determinant case has a zero eigenvalue at every momentum and fails the stated finite simple-root hypothesis.

[proved, implication for redesign] Reducing the nearest-cell hopping rank to satisfy (4), while changing only onsite fields, cannot retain C39/C42's 2n simple roots in an n-orbital Gaussian reference. For a fully Gaussian microscopic bus design, the project STATE and STATUS already record the settled quasi-free nearest-neighbor bound c<=log_2(d), including arbitrary periods (C1 and REPORT Section 4). That existing result rules out the all-n Gaussian replacement and is not re-derived here. Its locality assumption concerns genuinely local physical fermions; a general spin Hamiltonian need not have such a representation.

## 8. Exactly local dressing, with its range specified

[proved, theorem] Suppose a unitary U maps the proposed invariant space to a product-cell code, and conjugation by U enlarges the support of every local operator by at most R physical sites on each side. Suppose all cell lengths exceed 4R+2. If H is physical nearest neighbor, the connected bonds of the compressed UHU^dagger commute. In particular it cannot give the full-rank hopping of Section 3.

[proved, proof] A physical bond originally on sites x,x+1 is supported after conjugation on [x-R,x+1+R]. At a cut between 0 and 1, only original bonds with -R<=x<=R can cross the cut. Their union lies in [-2R,2R+1], at most 2R+1 sites on each side. The left and right such boundary strips inside any cell are disjoint under the length hypothesis. Every remaining transformed term is internal to one cell, and no term reaches across an entire cell. Group the terms into internal and adjacent-cell interactions and apply Section 2.

[proved, scope] A range bound independent of n therefore fails for the growing C39/C44 cells. Within this form of exact dressing, evading the theorem requires at least

    R >= (b_min-2)/4,                                  (13)

or abandoning the product-cell target description. The inequality is a necessary condition for evasion, not a construction. Arbitrary unitary changes confined separately to each cell are absorbed into G_b and never evade Section 2, regardless of their intracell depth.

[assumed, genuine remaining possibility] A cross-boundary dressing of physical range comparable to n is still finite range at each fixed n, and is allowed by the broad wording of the task. It is not excluded here. No such dressing that turns a repeated nearest-neighbor fixed-alphabet bond into the desired exact interacting fermion model, while also supplying the mobile wrapper and defect bounds, was constructed. A quasi-local dressing with nonzero tails also falls outside this exact-disjoint-support proof and returns to the matching question addressed by the other S runs.

## 9. Why the old perturbative construction is not contradicted

[proved, elementary source leakage] C39 has PVP=0 and a nonzero two-cell source. From its three orthogonal charge-source channels in each of the two buses,

    (V G_2)^dagger (V G_2)
       = 3(1+eta) q^(2n+2)/S_q^2 I,
    S_q = sum_(ell=0)^(n+1) q^(2ell).                  (14)

Here V is the unwrapped unit port operator and G_2=G tensor G. The three source factors are isometries into orthogonal species sectors, the two buses are orthogonal, and both endpoint vacuum amplitudes are nonzero. Thus QVP is nonzero. Making the already frustration-free parent commuting would not remove this algebraic fact while keeping the same source and code.

The readout is not invariant either. For the code state with only data bit j occupied, put p_j=q^(2j+2)/S_q. The readout has values one when the vacuum token is at that data site and zero otherwise. Its squared off-code norm is the variance

    ||Q R G|j>||^2 = p_j(1-p_j)>0.                     (15)

C39/C45 already retain these excursions in their error budget. Equations (14)-(15) identify why those ingredients cannot be promoted to exact invariance by terminology.

[proved, resolvents can increase effective rank] The physical-cut rank estimate applies to direct compression, not to a virtual effective Hamiltonian. For example take a product source from logical basis |i,j> into an auxiliary sector with energy a_i+b_j, where the positive a_i and b_j are distinct. Its second-order return operator is proportional to

    sum_(i,j) [1/(a_i+b_j)] |i><i| tensor |j><j|.

The Cauchy coefficient matrix has full rank. To see this without a determinant formula, if sum_i c_i/(a_i+b_j)=0 at m distinct b_j, multiply the rational function sum_i c_i/(a_i+x) by product_i(a_i+x). The numerator has degree at most m-1 and m distinct zeros, hence vanishes identically. Its residues then give every c_i=0. A product source has thus generated m independent Schmidt directions through the non-product resolvent. C39 uses virtual propagation to obtain its full-rank Gram matrices. The exact-product theorem does not apply to that mechanism.

[proved, an invalid shortcut made explicit] At the rooted-cell level one can define G_j G_(j+1) K_target G_j^dagger G_(j+1)^dagger, extended by zero off the two code spaces, and add positive code-complement penalties. This is an exactly invariant Gaussian reference by construction and keeps the finite momentum margins. One can similarly insert the exact intertwiners between G_b tensor G_(b+1) and G_(b+1) tensor G_b. As physical operators these act on whole adjacent cells. Their support grows with n, and the marker intertwiner must change untouched head amplitudes. They are not the required repeated two-site physical bond. Calling them nearest neighbor in cell coordinates would change the problem. No physical theorem is claimed for this shortcut.

## 10. Relation to published RG theorems and final ledger

[assumed, literature descriptions only] The supplied literature notes identify Mastropietro-Porta 2108.13708 and Porta-Scola-Singh 2509.08665 as nearby rigorous fermionic frameworks. C46 explains the mismatch between a microscopic finite-range fermion Hamiltonian and an unproved retarded, arbitrary-degree action, and separates correlation conclusions from the central-charge conclusion. C50 separately audited the primary RG route and still leaves its all-scale extension open. This run uses the supplied local descriptions only and makes no theorem-number application from them. No network was used.

[proved, application verdict] There is no redesigned microscopic bond here to which those primary theorem hypotheses have been checked. The failure occurs before RG. Exact invariance with separated product-cell endpoints contradicts the required hopping, and the unchanged mobile wrapper leaks. Writing down a logical finite-range quartic Hamiltonian alone does not provide its exact fixed-d physical realization.

| Requirement | Status |
| :-- | :-- |
| C39/C42 all-n reference, U(1), simple-root and nonresonance margins | [proved, imported] Retained as the target; their finite construction is not changed. |
| Commuting-projector or other exact invariant product-cell redesign with full-rank hopping | [proved impossible in this class] Sections 2-4. |
| Arbitrarily weak interacting correction to that reference in the same invariant product code | [proved impossible in this class] Sections 3 and 7, with the actual hopping remaining invertible. |
| Unchanged C44 marker move exactly preserving its selected code | [proved false] Section 5. |
| Boundary-local repair retaining the sine heads and two untouched interior positions | [proved impossible] Equation (10). |
| Exact cross-boundary dressing with physical range comparable to the cell length | [assumed, open construction] Not excluded by Section 8. |
| An alternative non-product code or different mobile selector | [assumed, open construction] No complete bond and exclusion audit supplied. |
| Published theorem implying S or c=n+1 for a new physical model | [not established] No applicable exact physical model was obtained. |

[proved, conclusion] This run closes the direct exact-product-cell proposal negatively, with a theorem whose finite hypotheses are checked for C39/C44. It does not close S. A successful redesign would need to leave this class, construct the cross-boundary encoded dynamics explicitly, and repair or replace the sine-head mobility while re-establishing the all-length exclusion and physical cyclic-gluing statements. These requirements are not supplied by the existing code gap or the all-n momentum margins.

## 11. Consistency checks and files

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` completes locally in under one second. It checks operator Schmidt rank 2n for full-rank hopping at n=1,2,3; invariance of the extracted transfer matrix under onsite and density-quartic additions; rank-one hopping after a charge-preserving two-qubit cell rotation; the determinant Fourier-degree bound; the identity (5) with arbitrary number-conserving higher-particle matrix blocks at n=1,2; (8) at b=4,5,8,16; and rank growth in a virtual Cauchy resolvent. The largest Hilbert matrix used is 64 by 64.

The commutator identities and quartic extraction errors are zero in the recorded run. The forbidden determinant Fourier coefficients are at most 2.13e-15. The head identity errors are below 4.8e-16. At b=4, the head leakage squared is 0.1701573216462289 and its analytic lower bound is 0.11909830056250527. These are controls of the algebra, not evidence for or against an interacting infrared central charge.

[proved, record] RESULT.md states the outcome, PROGRESS.md records the investigation, STATE.md is the local handoff, and checks.py/checks.json contain the finite controls. No project-level STATE update, edits outside work_C52, subagents, cluster computation, network access, or microscopic many-cell simulation was used. The project STATE is outside the permitted write scope.
