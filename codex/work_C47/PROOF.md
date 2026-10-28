# C47. What can be proved for the ten-letter demonstrator

[proved, convention] A tag governs its paragraph and the displays and continuation paragraphs following it, until the next tag or heading. A conditional implication tagged proved includes a proof of the implication, not a proof of its stated assumptions. Imported results not re-established here are tagged assumed. Numerical checks do not establish an infinite-volume statement. References to source availability describe this offline audit, not an exhaustive literature search.

[proved, outcome and scope] The results below do not establish c=1 for the ten-letter chain or c=2 for the n=1 wrapper. They establish the exact second-order XXZ Hamiltonian, a finite-order local unitary with a uniform spatial remainder bound, and a symmetry-preserving spectral embedding with a uniform comparison of the enlarged and original models. Section 5 verifies a conditional application of the recent local finite-cell matching result. Sections 6 and 7 explain why this still does not invoke a published XXZ infrared theorem. The microscopic calculations in Sections 1-5 concern the invariant period-three sector. The distinction from selection of that sector is addressed in Section 8.

## 1. Cell gap, code and symmetries

[proved] Use the literal bond of `../../contrib/routeA_demonstrator/cell.py`, with D0=D1=D=1, kappa=3/5, Jp=2, Jz=-1, mu=11/20 and K=2. Site type is conserved because every bond is block diagonal in P or B at each endpoint. In a period-three P B P sector, subtract the constant -mu M and block into M cells of dimension 32. The Hamiltonian is

    H_M = A + lambda V,
    A = sum_j A_j,
    V = sum_j V_(j,j+1),
    V_(j,j+1) = Jp(sigma+_(pR,j) sigma-_(pL,j+1) + h.c.)
                   + Jz Z_(pR,j) Z_(pL,j+1).

The physical ring length is 3M. The statements below hold at fixed cell parameters and fixed lambda, whenever the displayed smallness conditions hold. No parameter is chosen as a function of M or inverse temperature beta.

[proved] Write the pair basis as g=Phi+, a=Phi-, u=|01>, d=|10>. On the pair and its logical qubit, h_pair vanishes on g tensor C^2. On its six-dimensional complement the only off-diagonal entries are

    <a,0|h_pair|u,1> = kappa,
    <a,1|h_pair|d,0> = -kappa,

and their transposes; every diagonal entry is D. Its eigenvalues on that complement are D-kappa twice, D+kappa twice, and D twice. Therefore

    h_pair >= (D-|kappa|) q_pair,
    q_pair = (1-|g><g|) tensor I_c.

The two q_pair projections in a cell commute, since they act on disjoint port-partner pairs. Each h_pair preserves both projections, including the other side's projection, even though the checks share the logical qubit. With

    P_j = |g><g|_(pL,sL) tensor I_c tensor |g><g|_(sR,pR),
    Q_j = 1-P_j,

one obtains

    A_j >= (D-|kappa|)(q_L+q_R) >= (2/5) Q_j.        (1)

The kernel is exactly ran P_j and has dimension two. Putting one side in its D-kappa eigenvector and the other side in g attains 2/5, so this is the exact cell gap. Also ||A_j|| <= 2(D+|kappa|)=16/5.

[proved] The reduced port state in g is I/2. Consequently P_j O_port P_j=(tr O/2)P_j. Every single port Pauli matrix has zero compression, and P V P=0 for P=product_j P_j. One application of a port Pauli to g lies in its pair complement. Thus a bond V sends the global code into a state with both incident cells excited.

[proved] Define the cell charge and conjugation by

    Qcell = (Z_pL-Z_sL+Z_c-Z_sR+Z_pR)/2,
    Ccell = X_pL X_sL X_c X_sR X_pR.

The a-to-u conversion raises the pair charge by one in its reverse direction and lowers logical charge by one, and the a-to-d conversion does the converse. Each conversion is neutral. Port hopping is neutral, as is ZZ. Thus the full Hamiltonian commutes with total charge. Under Ccell, g maps to g, a to -a, u to d, and sigma+_c to sigma-_c. The relative minus sign in h_pair makes the two conversions exchange, proving charge-conjugation invariance. It sends Qcell to -Qcell. In the code these operations are logical Z/2 and X. The on-cell grading is

    Pi = (-1)^(1/2-Qcell) = Z_pL Z_sL Z_c Z_sR Z_pR.

There are 16 states of each parity, and one code state of each parity. All bond terms are even. These algebraic facts, rather than a numerical symmetry test, establish the selection rules.

## 2. Exact second order

[proved] Only the same interface can close a two-cell excitation at second order. For two different interfaces on a ring with at least three cells, some endpoint excited by the first interface is untouched by the second. A preserves the associated syndrome projection. Projection back to P therefore kills that contribution. Thus

    -P V Q A^(-1) Q V P = sum_j h2_(j,j+1),             (2)

where each h2 is an isolated-interface calculation. The inverse is on Q. The unused outer pair in either cell stays in g, so only two eight-dimensional half-cells are needed. The intermediate excited-excited subspace has dimension 36.

[proved] On these half-cells the source from logical bits c,e is

    V |g,c;g,e> = Jz |a,c;a,e>
                    + (Jp/2)(|u,c;d,e>+|d,c;u,e>).     (3)

The factors 1/2 follow from sigma+ g=u/sqrt(2), sigma- g=d/sqrt(2), and Zg=a. All formulas that follow allow D>|kappa| and real Jp,Jz.

[proved] Three entries of the excited resolvent suffice. Set

    A0 = (2D^2-kappa^2)/(4D(D^2-kappa^2)),
    B0 = 2D/(4D^2-kappa^2),
    C0 = -kappa^2/(4D(D^2-kappa^2)).                    (4)

For two active two-state blocks the Hamiltonian is 2D+kappa X_1+s kappa X_2, s=+1 or -1. In the X basis, its inverse has diagonal computational-basis entry

    (1/4)[1/(2D+2kappa)+2/(2D)+1/(2D-2kappa)] = A0.

For opposite signs s=-1, the matrix entry flipping both bits is

    (1/4)[2/(2D)-1/(2D+2kappa)-1/(2D-2kappa)] = C0.

One active block and one inert state have inverse diagonal entry

    (1/2)[1/(2D+kappa)+1/(2D-kappa)] = B0.

Two inert states have inverse 1/(2D). These identities follow by diagonalizing X_1 and X_2 and averaging their four eigenvalues. They account for every source sector in (3), because each half-cell charge is separately conserved by A.

[proved] For parallel logical bits the a-a source uses A0 and the two charged sources each use B0. For antiparallel bits the a-a source uses A0, one charged source has two inert states and the other has two active states. The diagonal entries are consequently

    Epar  = -Jz^2 A0 - (Jp^2/2) B0,
    Eanti = -Jz^2 A0 - (Jp^2/4)[1/(2D)+A0].             (5)

For the off-diagonal matrix element from 01 to 10, there are two overlaps of an a-a source with a charged source in the same conserved-charge block. Their combined source coefficient is Jp Jz. The effective Hamiltonian carries an overall minus sign, so the hopping is -Jp Jz C0. Therefore

    h2 = a0 I + t(sigma+ sigma- + sigma- sigma+) + u ZZ,
    a0 = (Epar+Eanti)/2,
    t = Jp Jz kappa^2/[4D(D^2-kappa^2)],
    u = (Epar-Eanti)/2
      = 3 Jp^2 kappa^4/[32D(D^2-kappa^2)(4D^2-kappa^2)]. (6)

For the last equality, putting 1/(2D)+A0-2B0 over a common denominator gives numerator 3kappa^4. This supplies the diagonal coefficient as well as the hopping, including the scalar term.

[proved] At the demonstrator parameters the exact result is

    a0 = -20505/11648,
    t = -9/32,
    u = 243/11648.                                     (7)

Staggering the logical sigma+ by (-1)^j changes the hopping to |t|. On an even ring it is periodic; on an odd ring it leaves a single twisted boundary bond. Writing S=Pauli/2 gives, after the staggering,

    lambda^2 h2 = lambda^2 a0 I
      + J [Sx Sx + Sy Sy + Delta_XXZ Sz Sz],
    J = 2|t|lambda^2 = (9/16)lambda^2,
    Delta_XXZ = 2u/|t| = 27/182.                        (8)

The boundary twist does not change this bulk algebraic identification. No thermodynamic claim has been inferred from it.

[proved, checks of the formula] At kappa=0, t=u=0 and a0=-(Jz^2+Jp^2/2)/(2D), giving -3/2 for D=1,Jp=2,Jz=-1. Every coefficient scales as inverse energy if D and kappa are scaled together with Jp,Jz held fixed. Conjugation forbids a field; U(1) forbids a pairing term. The coefficients in (6) satisfy both restrictions. These are exact consistency checks.

[proved, isolated-bond qualification] A Hermitian two-qubit operator commuting with Z_1+Z_2 and X_1 X_2 has the form in (6): charge conservation gives two one-dimensional charge blocks and a central two-dimensional block; conjugation equates the outer diagonal entries, equates the inner diagonal entries and makes the inner off-diagonal entry real. For a sufficiently small isolated two-cell coupling the four low eigenvalues are separated from the rest. Indeed ||V_bond||=3 and Weyl's variational eigenvalue inequalities give a separating gap at least 2/5-6|lambda|. The spectral projection and its polar identification with the code commute with both symmetries whenever the projection overlap is invertible, which holds in a neighborhood of zero by continuity. Its exact isolated-bond effective matrix is therefore XXZ plus a scalar. This does not identify the many-cell effective interaction with a sum of those isolated matrices.

[numerical] `checks.py` independently constructs the 36-dimensional rational resolvent using Fraction Gaussian elimination. Three exact fixtures reproduce (6), including kappa=0. It also computes the isolated third-order coefficient X^T V_QQ X, X=A_Q^(-1)VP, obtaining a nonzero hopping coefficient 1845/2048 at the demonstrator point. The third-order numbers in checks.json are finite-algebra controls; no infinite-series convergence is inferred. An independent floating-point 64-dimensional construction agrees with (7) within 1.8e-15.

## 3. A fully proved finite-order unitary, with its remaining leakage

[proved, theorem U] For rings with M>=3 and every fixed a>=0 define an interaction norm, using connected cell intervals as supports, by

    ||Phi||_a = sup_j sum_(X contains j) exp(a|X|)||Phi_X||.

For the model above put

    Delta = 2/5,
    v = 6 exp(2a),             aA = (16/5) exp(a),
    s1 = 2v/Delta,             t2 = 24v^2/Delta,
    s2 = 48v^2/Delta^2,
    lambda_* = min(1,Delta/(64v)),
    C_a = 12s2 v + 96(2s1 s2+s2^2)aA
                      +384s1^2 v+16384s1^3(aA+v).       (9)

For |lambda|<=lambda_*, there is an anti-Hermitian generator S=lambda S1+lambda^2 S2, of range at most three cells, for which the exact finite-ring identity is

    exp(S) H_M exp(-S) = D_M+R_M,
    D_M = A+lambda B1+lambda^2 B2,
    [Bq_X,P_X]=0,
    P D_M P = lambda^2 sum_j h2_(j,j+1),
    ||R||_a <= C_a |lambda|^3.                          (10)

Here P_X=product_(j in X)P_j. Constants do not depend on M or beta. The identity applies to the full 32^M-dimensional frozen-sector space. It is not a unitary equivalence to a 2^M-dimensional spin chain.

[proved, local construction] For a term O_X define its diagonal and off-diagonal parts relative to P_X by

    D_X(O)=P_X O P_X+(1-P_X)O(1-P_X),
    O_X^off=O_X-D_X(O_X).

Let A_X=sum_(j in X) A_j and let R_X be its inverse on 1-P_X, zero on P_X. By (1), ||R_X||<=1/Delta. Define

    L_X(O)=R_X O P_X-P_X O R_X.

For Hermitian O it is anti-Hermitian, has norm at most 2||O||/Delta, and direct multiplication gives

    [L_X(O),A]=-O_X^off.                               (11)

Spectral gap inversion here is only on one finite support, never across the extensive logical band. Let S1=sum_X L_X(V_X), B1=D(V), and

    T2=[S1,B1]+(1/2)[S1,V^off],
    S2=sum_X L_X(T2_X),       B2=D(T2).                 (12)

Terms in T2 are grouped on unions of intersecting bonds, so their supports have at most three cells. They are Hermitian. Expanding the conjugation through order lambda^2 and using (11) gives (10) at those orders, since

    [S1,V]+(1/2)[S1,[S1,A]] = T2,
    [S2,A]=-T2^off.

The construction commutes with charge, charge conjugation and cell translation, since P_X, A_X, their resolvents and the original terms transform covariantly. All ensuing terms preserve these symmetries.

[proved, code compression] In the global code, B1 vanishes. The code block of [S1,B1] vanishes because B1 preserves P and B1 P=0. On that block (1/2)[S1,V^off]=-PVQA^(-1)QVP. One may check the sign by using S1 P=A_Q^(-1)VP and P S1=-PV A_Q^(-1). In this equality the inverses are the finite-support inverses in (11); on a source from the global code they coincide with the global inverse because all spectator cells have zero energy. Equation (2) then proves the retained block in (10).

[proved, interaction estimates] If Phi and Psi have maximum support sizes r and s, respectively, their commutator has maximum size r+s-1 and

    ||[Phi,Psi]||_a <= 2(r+s)||Phi||_a||Psi||_a.         (13)

Only intersecting supports contribute. For contributions where the pinned site lies in X, sum the Y terms over the at most r sites of X, using exp(a|X union Y|)<=exp(a|X|)exp(a|Y|). This contributes at most 2r times the product of norms. Pinning instead in Y contributes 2s. Their sum gives (13), allowing harmless double counting. Block diagonal extraction and off-diagonal extraction each have operator norm at most ||O|| for Hermitian O: the former is a compression on orthogonal blocks, and the latter has norm equal to its off-diagonal rectangular block. Thus

    ||S1||_a<=s1,
    ||T2||_a<=12s1 v=t2,
    ||S2||_a<=s2.                                     (14)

The original bond has eigenvalues -1,-1,-1,3, so v=6exp(2a) bounds its incidence norm on the ring.

[proved, convergence and remainder] Set sigma=|lambda|s1+lambda^2s2. The restriction in (9) implies |lambda|s2<=s1 and 4sigma<=1/4. For any interaction Phi of range at most three, repeated use of (13) gives

    ||ad_S^m(Phi)||_a/m!
        <= binomial(m+2,2)(4sigma)^m ||Phi||_a.         (15)

At the kth commutator the support size is at most 3+2k, yielding the product (4sigma)^m(3)(4)...(m+2). The BCH series therefore converges in this interaction norm, uniformly in volume. For x<=1/2,

    sum_(m>=3) binomial(m+2,2)x^m <= 32x^3,

because dividing by x^3 leaves a positive-coefficient series increasing in x, whose value at 1/2 is 32. The BCH terms of degree at least three in the commutator expansion have norm at most 2048sigma^3(aA+v), hence at most the last term of C_a times |lambda|^3.

The pieces not already cancelled in the first two commutators are lambda^3[S2,V], the cross and S2-squared pieces of (1/2)ad_S^2(A), and (lambda/2)ad_S^2(V). Treating all generators and initial interactions as range at most three, (13) bounds a double commutator by 192 times the product of the three norms. Their bounds are, respectively,

    12|lambda|^3 s2 v,
    96|lambda|^3(2s1 s2+s2^2)aA,
    384|lambda|^3 s1^2 v.

Here |lambda|<=1 and sigma<=2|lambda|s1 were used. These are the other terms in (9), and prove the remainder bound. Since the equality is an operator identity, it holds at every temperature; the norm itself does not refer to temperature.

[proved, precise limitation] The construction enforces local block diagonality only in D. It does not establish [R,P]=0, much less R=PRP. Compression of R would discard real transitions, not constitute a unitary change of Hilbert space. Increasing a makes the spatial norm stronger but also changes lambda_* and C_a. Equation (10) is a finite-order reduction, not an all-order elimination theorem.

[assumed, imported thermal theorem] C46 PROOF Section 2 proves that for such a locally block-diagonal D, if J_B=sup_j sum_(X contains j)||B_X|| and delta=Delta-2J_B>0, then

    0 <= f_P(beta,M)-f_D(beta,M)
       <= beta^(-1) log[1+15 exp(-beta delta)].         (16)

Its assumptions allow defect motion and require no logical gap. We use that already proved result as an imported input, rather than rederive it.

[proved, application and limitation of (16)] Equations (9) and (14) imply J_B<=|lambda|v+lambda^2t2<=11Delta/512, so delta>=245Delta/256. Also ||R_M||<=M C_a|lambda|^3. The min-max principle bounds each ordered eigenvalue shift by this norm, hence

    |f_H-f_P| <= C_a|lambda|^3
          + beta^(-1)log[1+15 exp(-245beta Delta/256)]. (17)

The same norm bounds the ground-energy-density shift. Even after subtracting both ground energies, the bound contains 2C_a|lambda|^3. At any fixed nonzero lambda this is not o(beta^(-2)). Thus (17) cannot identify the desired thermal coefficient.

## 4. A charge- and conjugation-preserving massive embedding

[proved, theorem E] The demonstrator cell embeds into the Fock space of one logical mode and 60 auxiliary modes. The auxiliary Hamiltonian is quadratic, every auxiliary energy is at least g=1/5, and the embedding intertwines A_j, charge, grading and charge conjugation. The number 60 is a convenient symmetry-preserving bound, not a minimality claim.

[proved, construction] The 15 even excited eigenstates can be chosen as simultaneous eigenvectors |s,+> of A_j and Qcell. Put |s,->=Ccell|s,+>. They are an orthonormal eigenbasis of the odd excited subspace with the same energies E_s>=2/5 and opposite charges. Since an even state has Qcell=q_s=1/2+2k_s for an integer k_s, its partner has charge -q_s.

Introduce four distinct auxiliary modes a_s1,a_s2,b_s1,b_s2 for each s, all of energy E_s/2. Their charges are 1,2k_s-1,-1,1-2k_s. These are odd integers. With the logical occupation l=0,1, define

    J G|l> = |l> tensor |Omega>,
    J|s,+> = |0> tensor a_s1^dagger a_s2^dagger|Omega>,
    J|s,-> = |1> tensor b_s1^dagger b_s2^dagger|Omega>.   (18)

All image vectors are orthonormal. Their energy is E_s, their parity matches the preimage, and their charge under

    Qtilde = 1/2 - n_log + sum_alpha q_alpha n_alpha

is q_s or -q_s as required. Therefore J is an even isometry and H_mass J=J A_j. It preserves U(1), not merely its parity subgroup.

[proved, conjugation] Let U_aux exchange a_s1 with b_s1 and a_s2 with b_s2 for every s, with the same ordering inside both pairs. Its second quantization is an involution preserving the vacuum and H_mass, and reverses auxiliary charge. On the ordinary logical tensor auxiliary representation take Ctilde=X_log tensor U_aux. It exchanges the last two vectors in (18), acts as X on the auxiliary vacuum, and sends Qtilde to -Qtilde. Thus Ctilde J=J Ccell. In CAR coordinates its action on auxiliary annihilators includes the harmless minus sign from X_log reversing logical parity; it is still a linear mode permutation with signs and preserves the free auxiliary covariance. For a chain the Jordan-Wigner strings supply the usual cell-dependent signs. Working in each total-parity boundary condition retains locality of even interactions and all the bounds below. This verifies that the embedding need not sacrifice the demonstrator's charge-conjugation symmetry.

[proved, coefficient conversion] For an operator O on N fermion modes and field weight h>=1, its normal-monomial coefficient sum is at most

    (8h^2)^N ||O||.                                    (19)

There are 4^N occupation matrix units, with coefficients bounded by ||O||. A matrix unit is a product of creation operators, the vacuum projector product_i(1-n_i), and annihilation operators, up to a sign. Expanding the projector gives at most 2^N monomials. Normal ordering leaves degree at most 2N. Summing absolute coefficients and using the triangle inequality proves (19). Thus a complete-support operator on r enlarged cells costs at most (8h^2)^(61r). This constant is finite even though the embedding is impractical as a numerical representation.

[proved, uniform comparison with the original model] Extend every full-support interaction O_X by J_X O_X J_X^dagger. Let Pvalid_j=J_j J_j^dagger. These projections commute with H_mass and with every extended interaction, and every term touching an invalid cell vanishes there. All auxiliary-vacuum states are valid, so an invalid state has auxiliary energy at least g. Put

    w = sup_j sum_(X contains j)||(B+R)_X||,
    Dinv = 2^61-32.

In a sector with invalid set S, cut all interactions meeting S. Their total norm is at most w|S|. Let Eoutside be the ground energy on the remaining valid cells with these interactions cut. The cut all-valid system also has ground energy Eoutside, by putting every isolated cell in its code. Restoring the cut terms gives Evalid<=Eoutside+w|S|, whereas the invalid-sector ground energy is at least Eoutside+g|S|. Hence

    E_invalid(S) >= Evalid+(g-w)|S|.                    (20)

If w<g, the original and enlarged spectra below Evalid+g-w coincide, including multiplicities. This is a comparison with the full original frozen-cell model, not yet with a logical spin chain.

For the partition function the isolated valid one-cell heat trace is at least two, and the invalid trace is at most Dinv exp(-beta g). The cut all-valid trace factors with the same outside trace as the invalid sector. Adding back the cut interactions changes its logarithm by at most beta w|S|, by the same ordered-eigenvalue bound used for (17). Consequently Z_invalid(S)/Z_valid is at most [(Dinv/2)exp(-beta(g-w))]^|S|. Summing over subsets proves

    0 <= f_original-f_enlarged
       <= beta^(-1) log[1+(Dinv/2)exp(-beta(g-w))].      (21)

This proof is uniform in volume and beta and treats arbitrary summable supports. It also proves that the two models share any low-temperature power coefficient if that coefficient is independently established for either one.

## 5. What the resummed finite-cell theorem supplies

[assumed, explicit imported input M] The local note `../../contrib/direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md`, with the audit in its `../alternatives/IR_AUDIT.md`, supplies the following finite-cell effective-action result. It invokes Pedra and Salmhofer, "Determinant Bounds and the Matsubara UV Problem of Many-Fermion Systems", Theorems 2.4 and 4.5. Its primary text was not present among the supplied local sources and was not fetched. M, including its norm conventions, determinant bound, time-slicing limit and weighted extension, is therefore an imported assumption in this report. No primary-source theorem application is claimed independently of M.

M uses an enlarged Hamiltonian H_mass+H_log+K_D+Rtilde, where K_D annihilates the local auxiliary vacuum on both sides. In a complete anchored instantaneous coefficient norm that also bounds its pinned-field norm, let k and r bound K_D and the nonconstant part of Rtilde with input field weight H'=H+2delta. With alpha=2/(g-kappa_t), omega=2alpha/delta^2 and k+r<=1/(4omega), its auxiliary effective-action map obeys

    ||W_beta(K_D+Rtilde)-W_beta(K_D)||_(H;a_x,kappa_t)
         <= 4r,                                      (22)

uniformly in M and beta. At zero temperature W_infinity(K_D)=0. The norm has exponential spatial and circular-time diameter weights, arbitrary external logical degree, and retains the logical fields without contracting them. Scalars and source insertions are handled separately in M. Its finite-temperature baseline W_beta(K_D) is not asserted to vanish.

[proved, conditional application of M to this cell] Use the transformation of Section 3 with a=768, then the symmetry-preserving embedding of Section 4. Write H_log=PDP as the logical XXZ interaction acting on the logical tensor factors in every auxiliary sector. Define, term by term,

    K_D,X = J_X B_X J_X^dagger - (P_X B_X P_X)_logical,
    Rtilde_X = J_X R_X J_X^dagger,
    B = lambda B1+lambda^2 B2.                         (23)

The proof of M uses an even finite isometry, a quadratic positive auxiliary Hamiltonian, and inclusion of the whole logical auxiliary-vacuum space in the valid subspace. Equations (18)-(21) verify these properties for the symmetry-preserving embedding, so its state-pairing modification does not change that argument. The notation R_X here denotes a remainder interaction term, not the finite-support resolvent in Section 3. Since B_X is reducing for P_X, the first operator in (23) agrees with the logical extension on the auxiliary vacuum and has no off-diagonal matrix element out of that vacuum. Thus K_D,X annihilates that vacuum from both sides. Its norm is at most 2||B_X||. The extended remainder has norm ||R_X||. All charges, conjugation and translation are retained by (18) and the covariant transformation.

Choose

    H=2,          m_aux=60,       delta=2sqrt(60),
    g=1/5,        kappa_t=1/10,   a_x=1,
    alpha=20,     omega=1/6,      H'=2+4sqrt(60).        (24)

The massive energies are at least g, and their covariance is diagonal in cell and species. For 0<t<beta its absolute value is exp(-e t)/(1+exp(-beta e)). Its circular-time-weighted integral is at most 1/(e-kappa_t)<=10. This verifies the decay assumption entering M. The determinant estimate with delta in (24) is the imported part of M, not a consequence of this integral bound alone.

Equation (19) pays the complete coefficient conversion with F=61 modes per cell. There is ample spatial weight for the anchored or pinned multiplicities as specified in M: even allowing a reserve 12r in the exponent, the necessary weight is at most

    a_x + 61 log(8H'^2) + 12 < 768.                    (25)

For an elementary bound, H'<33, 8*33^2<exp(10), so the left side is less than 623. The reserve also bounds factors such as 122r^2 from marking or pinning a density. This follows from 122r^2<=exp(12r), r>=1; differentiating log(122r^2)-12r shows its maximum on this interval is at r=1, where it is negative. With the complete density norm and normalization convention of M, the input bounds can therefore be taken as

    k <= 2b_lambda,
    r <= C_768 |lambda|^3,
    b_lambda = |lambda|v+lambda^2 t2,                  (26)

using v,t2 at a=768. A separate scalar remainder density is bounded by C_768|lambda|^3 as well. Further fixed convention factors can be absorbed in the strict reserve in (25); no inverse logical hopping is used.

Impose the explicit additional restrictions

    |lambda|<=lambda_*,
    2b_lambda+C_768|lambda|^3 <= 3/2,
    b_lambda+C_768|lambda|^3 < 1/5.                    (27)

All left sides tend to zero with lambda. Hence a fixed positive lambda interval exists, independent of volume and temperature. The first two conditions verify M's smallness hypotheses. The third gives w<g in (20)-(21). This is a noncircular schedule: fix the cell and coefficient weight, compute the finite constants, then choose lambda.

Under M, (22) gives a zero-temperature logical action correction of norm at most 4C_768|lambda|^3. Relative to the hopping scale W=(9/32)lambda^2, this is at most

    (128 C_768/9)|lambda| -> 0.                        (28)

For any prescribed positive tolerance epsilon, also choosing |lambda|<=9epsilon/(128C_768) makes this relative bound at most epsilon. This proves the conditional finite-cell matching implication. The enormous constants prevent claiming that lambda=.05,.10,.15 satisfy (27). Nothing here validates those finite numerical couplings by these bounds.

[proved, reason this does not repeat C46's failed raw expansion] The order-one logical dependence of the original excited-cell Hamiltonian has been included in the exact one-cell eigenstates and the free massive energies. It is not Taylor-expanded in external logical fields during an excursion. The interactions K_D and Rtilde now go to zero with lambda in a finite coefficient norm. The zero-temperature vacuum cancellation in M uses two-sided annihilation, not just a vanishing vacuum expectation: every surviving normal monomial has an auxiliary creator and annihilator, so the earliest vertex kills the auxiliary vacuum. The resummation and representation have changed. This is consistent with the divergence of the original raw external-field expansion in C46.

[proved, limits even under M] The output of (22) is an action with kernels at different imaginary times and with arbitrary logical degree. It is not an instantaneous, finite-range spin Hamiltonian and is not a unitary equivalence of Hilbert spaces. At positive temperature W_beta(K_D) is a further term. C46 bounds its effect on the reference pressure but does not identify all its source kernels or its interaction with Rtilde. Equation (22) alone supplies neither a thermal o(T^2) estimate nor a current/stress scaling limit. Embedding an observable also changes its off-code expression; source matching must be included when identifying microscopic correlations.

## 6. Audit of the proposed XXZ theorem application

[assumed, local bibliography only] The supplied `../../contrib/ir_control/README.md` attributes spinless anomalous-exponent and beta-function results to Benfatto, Gallavotti, Procacci and Scoppola, CMP 160 (1994), and Benfatto and Mastropietro, CMP 258 (2005). It cites arXiv:1012.0403 for a response relation. It also summarizes Mastropietro and Porta, arXiv:2108.13708, Assumption 2.1 and Theorem 8.5, and Porta, Scola and Singh, arXiv:2509.08665. These are references reported by that review. No complete primary perturbed-XXZ theorem statement was located in the specified local material. Network access was forbidden and not used.

[proved, theorem-application verdict] No Benfatto-Mastropietro or related infrared theorem is invoked as proving this microscopic model's asymptotics. A theorem number cannot be supplied honestly as an application when its full hypotheses and the necessary input action have not both been verified. The following ledger states the actual matching boundary, rather than silently assuming an extension of a published result.

| Required input or conclusion | What has been established here |
| :-- | :-- |
| One spin-1/2 per retained cell | [proved] Exact code in Section 1. |
| XXZ reference and its normalization | [proved] Equations (6)-(8), Delta=27/182. |
| U(1), conjugation, and cell translation | [proved] Sections 1, 3 and 4; the conditional matching can preserve these symmetries. |
| Small norm relative to the logical hopping | [proved conditional on M] Equation (28) for the matched zero-temperature action. |
| An exact instantaneous small perturbation of XXZ | [assumed, missing] The finite-order Hamiltonian still leaks, and M produces a retarded action. |
| Applicability to arbitrary-degree retarded vertices | [assumed, missing] No primary theorem or proved extension with these precise inputs is supplied. |
| Inclusion of Delta=27/182 in a quantified weak-coupling region | [assumed, missing] A qualitative smallness neighborhood of the free XX chain does not certify this number. |
| A treatment of half filling and umklapp | [assumed, missing] It must be checked in the chosen primary theorem. |
| Physical source/counterterm matching | [assumed, missing] Effective-action norm control alone does not construct the lattice current or stress limit. |
| Low-temperature pressure or a complete c=1 scaling limit | [assumed, missing] Neither is proved by Sections 1-5. |

[proved, a fixed anisotropy is not made smaller by lambda] In Jordan-Wigner variables the nearest-neighbor density coupling is 4u and the hopping magnitude is |t|. Therefore

    V_density/|t| = 4u/|t| = 27/91,

independent of lambda. Taking lambda smaller suppresses (28) but does not reduce the reference interaction. Thus an existential theorem valid for |Delta|<epsilon_0, with no usable bound on epsilon_0, cannot by itself verify Delta=27/182. This is a quantifier issue, not evidence that this anisotropy is physically unstable.

[proved, half filling is a separate check] The staggered free XX reference at zero field has Fermi points +pi/2 and -pi/2. Twice the difference is 2pi, so the repeated-pair umklapp momentum is a reciprocal lattice vector. It cannot be discarded by a nonresonance assumption for generic separated seas. Antisymmetry can force derivative factors in a spinless quartic; establishing their irrelevance in an interacting RG is an additional step, not a nonzero momentum margin. A suitable half-filled spinless theorem may address this, but its statement must actually be checked.

[proved, an optional in-family adjustment] For D=1,Jp=2,Jz=-1 and variable 0<kappa<1, (6) gives

    Delta_XXZ(kappa)=3kappa^2/[2(4-kappa^2)],
    |t(kappa)|=kappa^2/[2(1-kappa^2)].                  (29)

Thus anisotropy can be made arbitrarily small without changing the alphabet, while keeping a positive hopping at each fixed kappa. One must choose kappa first and lambda afterwards, since the relative remainder contains 1/|t(kappa)|. At kappa=0 the entire logical second-order dynamics disappears. This adjustment can meet an existential weak-coupling radius once a precise theorem is available; it does not preserve the requested anisotropy near .15 and does not prove the missing action-to-RG step.

[heuristic] The exact U(1), conjugation and translation, the easy-plane reference and the small conditional matching error are good reasons to expect the observed Luttinger liquid to persist. That expectation is not substituted for the missing theorem in the ledger.

## 7. Unitary equivalence, temperature and central charge

[proved] A spatial exponential bound and a fixed-coupling infrared theorem are different statements. For an interaction R with ||R||_a<=epsilon, truncating intervals longer than ell leaves local norm at most epsilon exp(-a ell). At fixed ell this is a nonzero constant as T tends to zero. The eigenvalue comparison that led to (17) therefore cannot make its thermal error o(T^2). Increasing ell with beta can improve a comparison, but produces a temperature-dependent approximation, not the requested one fixed finite-range Hamiltonian. It also supplies no RG theorem for the remaining interaction.

[proved] An on-cell excitation gap does not isolate the entire logical band. A concrete example is an independent auxiliary two-state system of energy Delta at each cell, tensored with lambda^2 times the XXZ chain (7). Parallel and alternating Z-product states have expectations differing by 2u lambda^2 M on an even ring. Hence the logical spectral width is at least that large. The lowest state with one auxiliary excitation is only Delta above the logical ground. Once 2u lambda^2 M>Delta, the logical band is not a cluster of the lowest 2^M levels separated below all auxiliary levels. This example satisfies an exact on-cell gap and even has an exact tensor factorization. It disproves using that gap as a volume-independent justification for global-band perturbation theory. It does not disprove a more restricted low-energy construction for the demonstrator.

[proved] The elementary identity f(T)=e0-pi c T^2/(6v)+o(T^2), if proved together with a known single speed v and a conformal interpretation, would identify c. Equations (17), (21) and (22) do not supply its needed interacting logical coefficient. Only (21) transfers a coefficient that is independently known for the enlarged problem. For the fixed period-three sector, physical free energy is f_cell/3 and physical velocity is 3v_cell, so the same c would be obtained in physical-site units. This conversion follows from L=3M and does not add a scaling-limit assumption.

[proved] Normalized correlations alone can lose central-charge information. Tensor a system with an independent neutral critical system and take observables acting only in the first factor. Their normalized correlation functions are unchanged, whereas the free energies add. Therefore a collection of charge-correlation asymptotics, without a proof of completeness of the massless degrees of freedom or of the stress tensor, does not logically identify the total central charge. C46 PROOF Section 8 gives a separate conditional current-algebra lower bound under a CFT assumption. It cannot be promoted to equality for this microscopic model by counting one exact U(1).

[assumed, sufficient but unproved replacement] An infrared theorem for the action in Section 5, including its thermal winding terms and microscopic sources, that proves

    f_frozen(T)=e_frozen-pi T^2/(6v_cell)+o(T^2)

with the identified mode speed, or constructs its complete compact-boson stress/current scaling limit, would close the frozen-sector c=1 inference. Such a statement is not derived or invoked here. Exact all-order quasi-local unitary elimination is stronger than the Euclidean route and is likewise not established.

## 8. Selection of the frozen sector

[proved, zero-coupling selection] At lambda=0, write a type sector as B's separated by runs of k P's. For k=0 the run energy is K=2. For k>=2 the two endpoint checks can simultaneously vanish. For k=1, both partners claim the same port. Let P1 and P2 be their Bell projections. Direct Bell contraction gives P1 P2 P1=P1/4 on ran P1. The singular overlap is consequently 1/2, and diagonalizing the two-projection blocks gives ||P1+P2||=3/2. Thus q1+q2=2-P1-P2>=1/2. The lower bound h_pair>=(2/5)q from Section 1 gives a shared-port cost at least 1/5, independently of the spectator logical qubits.

Subtract the period-three trial energy -mu L/3. For one B followed by k P's the lower-bound excess is

    e_k+mu(k-2)/3 >=
       49/30,             k=0,
       1/60,              k=1,
       0,                 k=2,
       11(k-2)/60,        k>=3.                         (30)

An all-P sector has zero energy and is higher than -mu L/3. Summing the positive run bounds proves that on L divisible by three, the ground sectors at zero coupling are precisely the three translates of P B P with their logical codes. This proves a zero-coupling selector statement without numerical run eigenvalues.

[proved, no fixed per-site separation from every other type sector] At any fixed bounded nearest-neighbor bond, start from a ground state in a selected type sector. Trace out one P site and replace it by any B state. All nonincident bond expectations are unchanged, and each of the two incident expectations changes by at most 2||h||. The new density matrix is entirely in a different type sector, so that sector's ground energy is at most E_selected+4||h||. The difference per site is at most 4||h||/L. Thus the README's headline wording of a positive per-site separation from every other type configuration cannot be literal uniformly in L. A positive cost per defect is compatible with this observation.

[proved, limitation of the supplied run comparison] Lower bounds of the form E_sector>=L u_2+delta N_defect and an upper bound E_selected<=L u_2+epsilon L imply exclusion when delta N_defect>epsilon L. They do not exclude a single defect for arbitrary L when epsilon>0. This is the arithmetic reason that comparing all run lower bounds with u_2 does not alone prove a dressed dilute-defect gap. The later demonstrator prose acknowledges the density error, but its assertion about dilute defects requires a separate comparison. Sections 1-7 avoid importing that assertion by working in the exactly invariant period-three sector. No claim that the selected phase is actually unstable follows.

[assumed, full-chain transfer still needed] To pass a frozen-sector infrared theorem to the unrestricted ten-letter chain at nonzero lambda, prove a ground-relative positive dressed cost of type defects and control their thermal contribution. That proof is not supplied here. The reported numerical type selection is evidence for it, not a replacement. The zero-coupling bound (30) is a starting point, not the completed nonzero-coupling result.

## 9. The n=1 member of the 4464-letter wrapper

[assumed, imported finite construction] C44 and C45 give two retained lengths b=4 and b+1=5, one payload qubit and one length qubit, complete auxiliary exclusions, second-order common-bond matching and a finite-order local remainder. C46 supplies reference thermal and physical-quotient comparisons. Those results remain the imported finite and reference inputs; they are not rederived in this file.

[proved, why the finite-cell application above does not cover the wrapper] Its rooted microscopic cell space includes arbitrarily long legal cells. A finite cutoff is not invariant: marker motion transfers length between neighboring cells, so a configuration at the cutoff can leave it. The spectral embedding and (19) cost a factor exponential in the number of excited one-cell states. The auxiliary vacuum projector itself has coefficient sum (1+h^2)^m. As the cutoff increases, m increases without a uniform bound. Hence the finite constant argument (25)-(27) has not produced a uniform infinite-length estimate. Letting lambda depend on the cutoff and tend to zero would not establish the fixed-Hamiltonian statement.

[proved, additional infrared distinction] Even a successful matching would retain two interacting modes, the payload and length fluid. An ordinary single-spin XXZ theorem would at most address one component or a decoupled reference, not the mixed terms present at finite wrapper parameters. An infrared theorem for the full two-mode action and its physical quotient would still be needed. The n=1 specialization removes the large number of payload flavors, but not these two obstructions.

[assumed, outstanding wrapper conclusion] c=2 for the full n=1 wrapper remains conditional on massive matching uniform in the length space, a suitable two-mode RG theorem and a complete central-charge identification. No new wrapper central-charge measurement or proof is claimed.

## 10. Verification and theorem boundary

[numerical] Run `python3 -B checks.py`. The control uses only the standard library, numpy and scipy. Exact rational checks cover the resolvent identities at three fixtures. Numerical matrices have dimension at most 64. They give two cell zero modes, gap 0.3999999999999997, exact floating-point zero charge and conjugation commutators, and second-order discrepancy 1.7763568394002505e-15. The isolated two-cell remainder divided by lambda^3 is 3.1954, 3.2471, 3.2740 at lambda=.02,.01,.005. These include the scalar correction and therefore need not match the traceless residual quoted in the older README. The new order-two unitary control has code-block discrepancies below 5.2e-16 and remainder/|lambda|^3 approaching 10.775 as lambda decreases. The results are finite-cell controls only.

[proved, final boundary] Equations (1)-(15), the embedding (18)-(21), and zero-coupling selection (30) are established here with proofs. Equations (17) and (22)-(28) are applications conditional on the named imported thermal and massive-matching inputs. The exact anisotropy is 27/182. Neither the requested all-order unitary reduction nor an applicable published perturbed-XXZ infrared theorem has been established. Consequently c=1 for the full ten-letter model and c=2 for the n=1 wrapper are not theorems of this work.
