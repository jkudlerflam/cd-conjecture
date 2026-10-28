# C34. The open qubit strata after the second reduction

2026-10-28. No full proof of c <= 1 and no counterexample.

[proved] This run adds an elementary gap bound for an easy-axis family with a misaligned Dzyaloshinskii-Moriya vector. It includes time-reversal-only Hamiltonians, reflection-only Hamiltonians, and Hamiltonians with no uniform internal symmetry. It also proves that arbitrary site-dependent spin rotations cannot enlarge the direct applicability of the supplied Bravyi-Gosset-Liu-Wong gap theorem to a uniform periodic bond and field.

[proved] These results leave bare search regions of maximal dimensions five in each antiunitary branch and eight in the symmetry-free branch. Section 5 gives an exact algebraic partition and the remaining inequalities for the certificates used here. These are dimensions of search regions, not dimensions of unknown critical manifolds.

[heuristic] A proof for the full family is not yet within reach by combining the current exclusions. The missing step is still a dynamical exhaustion argument, not one last isolated junction calculation. In particular, the available results do not exclude a strongly interacting neutral critical sector at finite velocity. The unitary Z2 stratum also remains incompletely classified in C15, so settling the three strata here alone would not finish the full family.

[proved, scope] An exact set of all points left unresolved by every physics argument and numerical test cannot be extracted from the supplied work. The continuum neighborhoods have no specified bare-coupling radii, and finite scans have no uniform crossover bound. This report therefore distinguishes an explicit algebraic remainder from the smaller, presently uncharacterized remainder after all controlled continuum arguments. Claiming that the latter is a finite collection of known low-dimensional critical sets would be unsupported.

Tags apply to the paragraph, derivation, or table introduced by them. [proved] includes finite-dimensional deductions from an identified input theorem. [assumed] records external theorem input or continuum assumptions. [numerical] records actual sizes and diagnostics. [heuristic] includes proposed tests and interpretations not established for the whole family. No outside source was fetched and no thermodynamic simulation was performed here.

## 1. Conventions and inherited exclusions

[proved, parametrization] Write

    H = sum_j [sigma_j^T J sigma_(j+1) + B.sigma_j],
    J = S + A(D),     S^T=S,     A(D)_ab=epsilon_abc D_c.

Thus D_z multiplies XY-YX. Remove constants, bond coboundaries, uniform proper spin rotations, and positive scale. One can fix

    tr(S^2) + 2|D|^2 + |B|^2 = 1.

The generic quotient has dimension 12-3-1=8. Negative scale is not an equivalence of ground-state problems. All finite-ring statements below use periodic L>=4.

[assumed, inherited exact-model input] The quadratic Jordan-Wigner class has c<=1 by the supplied free theorem. Uniform U(1) reduces to twisted XXZ and zero-field symmetric exchange to XYZ, with c<=1 under their supplied complete infrared classifications. The additive-charge classification in `../../contrib/d2_gap_reduction/additive_charge.md` leaves only the additional staggered-charge family

    a(XX-YY)+u ZZ+b Z.

Its still-open two-parameter region, at a=1 and b>0, is

    -1 < u < sqrt(1+b^2/4),     u != 0.

The two outer wings and the quadratic boundary are already settled in that contribution. This family has a uniform unitary pi rotation, so it does not lie in the strata requested here. It remains relevant to the full-family verdict.

[proved, inherited mechanism exclusions] C15-C17 exclude literal independent physical sublattices at range one, a conserved one-site qubit marker with an active walker, a longest frustration-free segment above three sites, and the microscopic SU(2) higher-level route. Exact U(1) rules out independently movable atom and pair minima and reduces the many-body problem to XXZ. The completed Abelian parton count gives c<=1 only when its current-completeness and coset assumptions hold. None of these statements excludes a neutral sector, dressed composites, or frustrated quantum-selected cells in the three present strata.

[assumed, controlled extension] A complete one-boson or at-most-two-Majorana continuum description remains sufficient in its regular scaling regime. The relativistic c-theorem applies only when the complete ultraviolet continuum theory and its hypotheses have been established. Unequal velocities, a new soft branch, a competing vacuum, or a vanishing matching scale require further analysis. A path through bare Hamiltonians is not an RG trajectory.

## 2. Exact symmetry tests

[proved, reproduced from C6 with the relevant argument] Uniform internal symmetries act on Pauli vectors as

    G_u = {R in SO(3): R S R^T=S, R D=D, R B=B},
    G_a = {O in O(3): det O=-1, O S O^T=S, O D=-D, O B=B}.

The extra sign on D follows from its axial-vector transformation. If G_u is trivial, two distinct antiunitary Pauli actions cannot both occur, because their product would be a nonidentity proper action. The surviving improper action squares to the identity. An improper orthogonal involution in three dimensions is either -I or a plane reflection. Consequently the three target strata are

    T: B=0, G_u={I};
    K: B!=0, G_u={I}, and there is n with
       S n parallel n, D parallel n, B.n=0;
    N: G_u={I}, G_a empty.

T is time reversal, squaring to -1 on one spin. K is a reflection of one spin component, squaring to +1. Spatial symmetries are not included in these names. Translation is always present.

[proved] For S with distinct eigenvalues, diagonalize S. Its proper stabilizer consists of the identity and the three principal-axis pi rotations. Its full orthogonal stabilizer consists of eight diagonal sign choices. This is a finite exact test of the displayed equations. With a repeated eigenvalue, use the invariant equations themselves, retaining rotations in the degenerate plane. A choice of principal frame must not create a spurious symmetry.

[proved] The quadratic, uniform-U(1), and zero-field-XYZ loci all have nontrivial G_u. They are disjoint from T, K, N. The supplied classification of arbitrary additive on-site charges implies that no nonzero such charge exists on these target strata, including site-dependent charges on even rings. This is a statement about microscopic additive spin densities, not emergent or multispin currents.

## 3. What the supplied gap theorem can cover

### 3.1 A finite linear-algebra certificate

[assumed, supplied theorem] We use the Bravyi-Gosset-Liu-Wong theorem exactly as quoted in `../../contrib/d2_gap_reduction/rank_one.md`. On an arbitrary finite graph, a bond with ferromagnetic longitudinal coefficient K>=||A_transverse||op, no mixed longitudinal-transverse block, and longitudinal fields all at least mu>0 has a unique ground state and gap at least 2 mu. Transverse fields are unrestricted. The original preprint was not independently accessed in this no-network run.

[proved, conditional on that theorem] Define

    s = ||J||op,
    E_-(J) = ker(J+s I) intersect ker(J^T+s I),
    g(J,B) = 2 || P_E B ||,

where P_E is the Euclidean orthogonal projector onto E_-(J). A nonzero g certifies a unique ground state and gap at least g for every L. If J=0, E_- is all of R^3, giving the correct independent-spin gap 2|B|.

To prove this, a unit axis n of the supplied criterion obeys

    Jn=J^Tn=lambda n,
    lambda <= -||P J P||op,       P=I-nn^T.

Its cross blocks vanish. The largest singular value of the resulting block matrix is max(|lambda|,||PJP||op), so lambda=-s. Conversely, every unit n in E_- has transverse norm at most s and satisfies the criterion. Maximizing |B.n| over this subspace gives ||P_E B||. The certificate is exact; a floating-point nullspace near a singular threshold is only a diagnostic.

[proved] If D!=0, subtracting the two eigenvector equations gives A(D)n=0, so n is parallel to D. A certificate therefore exists precisely when

    S D = lambda D,
    lambda <= -||J restricted to D-perp||op,
    B.D != 0.                                             (1)

For D=0 it exists precisely when the eigenspace of S at eigenvalue -||S||op has a nonzero field projection. This formulation includes all eigenvalue degeneracies.

[proved] In axes with D=d e_z and S=diag(x,y,z), the transverse norm is

    t = || [[x,d],[-d,y]] ||op
      = sqrt(((x+y)/2)^2+d^2) + |x-y|/2.                   (2)

For completeness, put m=(x+y)/2 and q=(x-y)/2. The two eigenvalues of T T^T are (sqrt(m^2+d^2) +/- |q|)^2. Taking the larger square root gives (2). The gap condition is z<=-t and B_z!=0. No diagonalization of a many-body Hamiltonian is needed.

### 3.2 Site-dependent axes do not enlarge this certificate

[proved, new finite-dimensional lemma] Suppose an arbitrary product of on-site spin rotations transforms the uniform periodic H(J,B) into the hypotheses of the supplied theorem, with positive longitudinal fields bounded below by mu. Then a uniform rotation already supplies a certificate with at least that lower bound.

Let n_j be the physical unit vector chosen as the longitudinal axis on site j. Choose its sign so that B.n_j<=-mu. The absence of mixed blocks on bond (j,j+1) implies

    J n_(j+1) = -K_j n_j,
    J^T n_j = -K_j n_(j+1).

Rotations preserve singular values. Dominance of K_j over its transverse block implies K_j=||J||op=s on every bond. Average the equations around the ring, writing nbar=L^(-1)sum_j n_j. Periodicity gives

    J nbar = J^T nbar = -s nbar,
    B.nbar <= -mu.

Thus nbar is nonzero, lies in E_-(J), and has norm at most one. The unit vector nbar/|nbar| gives longitudinal field at least mu/|nbar|>=mu. Section 3.1 applies. This proof also covers repeated singular values, rotating axes, and size-dependent choices of axes.

[proved, scope] This is closure under product spin rotations for direct application of the quoted theorem to the full bond. It is not closure under blocking, entangling circuits, a new comparison theorem, or arbitrary perturbations of gapped models. It prevents another search for a helical or staggered choice of axes from evading (1). It does not imply that Hamiltonians failing (1) are gapless.

### 3.3 Consequences for the requested strata

[proved] The direct theorem gives g=0 everywhere on T, because B=0. On K with D!=0, its only possible axis is D, whereas B.D=0, so again g=0. On N with D!=0 and D not an eigenvector of S, E_- is zero and g=0. These failures hold even when the exchange norm is large and ferromagnetic-looking.

[proved] The theorem does remove regions on the D=0 part of K and N, and on the aligned-D part of N. At a strictly negative dominant eigenvalue it allows arbitrary transverse fields. In particular it closes the mixed-field ferromagnetic Ising problem, as the supplied contribution already established. It does not close the antiferromagnetic mixed-field problem, since staggering the longitudinal axis changes the sign of the uniform field on alternate sites.

[assumed] Suitable stability theorems or controlled gapped perturbation theory can extend certified phases beyond their algebraic slices. Such extensions need their own hypotheses and quantitative neighborhood. No unspecified stability radius is used below to remove bare parameter points.

## 4. A new gap bound beyond the common-axis criterion

### 4.1 Domain-wall inequality

[proved, new] Consider

    H = sum_j [a(XX+YY) + c ZZ + d_z(XY-YX)
                       + d(YZ-ZY) + b Z],
    c<0,
    r=sqrt(a^2+d_z^2),
    kappa=|c|-r-2|d|.                                     (3)

A uniform rotation about z makes the perpendicular component of D point along x, so this includes arbitrary D at fixed uniaxial S. If b=0 and kappa>0, the ground space is exactly the two z-polarized product states and

    gap above that doublet >= 4 kappa.                    (4)

If b!=0 and kappa>=0, the favored polarized state is the unique ground state and

    E1-E0 >= 2|b|.                                        (5)

Both statements hold for every periodic L>=4. They do not use the external gap theorem.

Here are all steps. Define the commuting wall projectors

    w_j=(1-Z_j Z_(j+1))/2,      W=sum_j w_j.

The Ising energy is cL+2|c|W. The longitudinal-DM and transverse symmetric exchange on bond j have an odd-spin block

    [[0,2(a+i d_z)],[2(a-i d_z),0]],

up to reversing the basis order. Its eigenvalues are +/-2r and it vanishes on aligned spins. Therefore its operator lower bound is -2r w_j, and the sum is bounded below by -2r W.

The perpendicular DM contribution can be regrouped without approximation:

    d sum_j (Y_j Z_(j+1)-Z_j Y_(j+1))
       = d sum_j Y_j (Z_(j+1)-Z_(j-1)).

Let Q_j=(1-Z_(j-1)Z_(j+1))/2. The jth summand is supported on Q_j, commutes with Q_j, and has eigenvalues +/-2|d| there. Hence it is bounded below by -2|d| Q_j. The diagonal identity

    Q_j = w_(j-1)+w_j-2 w_(j-1)w_j
         <= w_(j-1)+w_j

gives sum Q_j<=2W. Combining the bounds proves the full operator inequality

    H-cL >= 2 kappa W + b sum_j Z_j.                       (6)

For b=0, both fully polarized vectors are annihilated by every exchange contribution other than c ZZ. On a periodic binary configuration the number of walls is even. It is zero only on these two vectors and at least two on their orthogonal complement. Thus

    W>=2(I-P_ferro),
    H-cL>=4 kappa (I-P_ferro).

This proves the ground-space claim and (4), without assuming that W commutes with the entire Hamiltonian.

For b!=0 subtract E_ref=cL-|b|L. The field term becomes 2|b|N_bad, where N_bad counts spins opposed to the favored polarization. It satisfies N_bad>=I-P_favored. For kappa>=0, (6) yields

    H-E_ref >= 2 kappa W+2|b|N_bad
             >=2|b|(I-P_favored).

The favored product vector has energy E_ref, proving uniqueness and (5).

[proved, limits] At d=0, b=0 the bound is 4(|c|-sqrt(a^2+d_z^2)), consistent with the easy-axis ferromagnetic XXZ wing after the axial twist. At a=d_z=d=0 it equals the exact Ising domain-wall-pair cost 4|c|. At kappa=0 and b=0 it supplies no positive gap and no conclusion about the infrared. The coefficient 2|d| is a sufficient estimate, not a claimed sharp phase boundary.

### 4.2 Which open strata this actually settles

[proved] Let a!=c and d d_z!=0. The proper stabilizer of S=diag(a,a,c) consists of axial rotations and pi rotations about transverse axes. No nonidentity one fixes a vector D with both parallel and perpendicular components. At b=0, the model therefore belongs to T. At b!=0, any antiunitary would have to be a plane reflection with normal parallel to D and an eigenvector of S. D is not such an eigenvector. The model belongs to N. The supplied common-axis certificate fails for both, but (4)-(5) apply.

[proved] If d!=0, d_z=0 and b!=0, the same model is reflection-only. The reflection normal is e_x, D is parallel to it, and B is perpendicular to it. There is no nonidentity proper symmetry preserving both D and B. For parameters satisfying (3), this closes an interacting subfamily of K with D!=0, where the supplied theorem has no certificate.

[proved, dimensions] After rotations and scale, the uniaxial family (3) has four coordinates. Its b=0 time-reversal-only part has dimension three. Its d_z=0, b!=0 reflection-only part has dimension three. Its generic b!=0, d d_z!=0 symmetry-free part has dimension four. The inequalities in (3) contain relatively open regions of these respective families. They are lower-dimensional inside the maximal 5,5,8 target strata.

[proved, extra localization of T] All uniform product eigenstates in the time-reversal-only stratum occur within the uniaxial family just considered, although the gap bound need not hold there. To see this, rotate the proposed product vector to |0>^L. Its adjacent-double-flip amplitude vanishes only if S_xx=S_yy and S_xy=0. Its single-flip amplitude at B=0 is 2(S_xz+i S_yz), requiring S_xz=S_yz=0. Antisymmetric terms cancel in the single-flip sum. Therefore S=diag(a,a,c). If S is scalar, D supplies a surviving continuous proper symmetry, so it is outside T. If D lies entirely along or entirely perpendicular to the unique axis, a proper symmetry also survives. Thus genuine T product eigenstates require a!=c and d d_z!=0. A simple-spectrum S in T has no uniform product eigenvector at all. This does not exclude nonuniform product helices or dressed vacua.

[numerical, convention checks] `checks.py` checks the complete matrix inequality (6), the ground-space energy, and the gap bound for six choices of (a,c,d_z,d,b) at L=4,6,8. This includes both signs of b, the reflection-only case, and kappa=0 with b!=0. The smallest remainder eigenvalue is no less than -3.6e-15, consistent with roundoff. For (a,c,d_z,d,b)=(0.4,-1.4,0.2,0.3,0), the analytic bound is 1.411145618 and the gaps at L=4,6,8 are 2.878407, 2.952049, 2.955272. These are checks of a proof, not central-charge estimates.

## 5. An explicit remainder and its dimensions

### 5.1 The exact algebraic envelope

[proved, definition and partition] Let G={g(J,B)>0} be the certificate set of Section 3, and let W_gap be the union of the globally rotated families in Section 4 satisfying b=0, kappa>0 or b!=0, kappa>=0. Define

    R_alg = (T union K union N) \ (G union W_gap).           (7)

Equation (7) is an exact, checkable remainder for these two certificates. It is an outer search set, not a claim that every member lacks every other proof. The following rows are disjoint and cover it. Every row retains the G_u/G_a restrictions of Section 2. Dimensions are maximal dimensions after uniform rotations and positive scale; they are not conformal codimensions.

| Piece | Bare conditions, in addition to its named symmetry | Dimension | Additional certificate removal | Concrete next test [heuristic] |
|---|---|---:|---|---|
| T_s | B=0, S has three distinct eigenvalues | 5 | Neither G nor W_gap meets it | Track neutral low gaps in all spin and bond channels while varying the three components of D; determine whether a new tower has finite velocity |
| T_d | B=0, S has a double eigenvalue; choose S=diag(a,a,c), D=(d,0,d_z), d d_z!=0 | 3 | Remove c<0 and |c|>sqrt(a^2+d_z^2)+2|d| | Track the first soft wall or composite on leaving the proved ordered wing; test dispersion, conversion amplitudes, and independent towers |
| K_D | D!=0, SD parallel D, B.D=0 | 5 | Remove the three-dimensional reflection-only wall-bound subset of Section 4 | Continue the tilted-field Ising/liquid investigation with nonzero D, then vary both transverse eigenvalues; test independent Ising mass and liquid velocity |
| K_0 | D=0, B!=0, a reflection survives | 4 | Remove P_E B!=0 | Complete the tilted-field Neel boundary first, then vary anisotropy and field direction; entropy plus momentum-resolved towers |
| N_m | D!=0, SD not parallel D | 8 | Remove the four-dimensional wall-bound subset of Section 4 | Locate any extra neutral or composite gap closing along a calibrated liquid; measure its response to all independent bond directions |
| N_a | D!=0, SD=lambda D | 6 | Keep only lambda> -||J restricted to D-perp||op | Tilt the field and vary the transverse exchange through the dominance boundary; test gap closing and each surviving velocity |
| N_0 | D=0, no internal symmetry | 5 | Keep only the cases with no eigenvalue -||S||op | Search the positive-dominant and competing-sign anisotropy regions with all three field components; resolve translation-odd neutral order as well as spin correlations |

[proved, simplifications in the table] In N_a, B.D cannot vanish, since that would give a reflection antiunitary. B cannot be parallel to D, since that would give a proper pi rotation. Hence the field condition of (1) is automatic and only the dominance inequality needs testing. In N_0, S must have three distinct eigenvalues and all three principal components of B must be nonzero. If a component vanished, reflection in its normal would be a symmetry. If S were degenerate, one could choose a reflection normal in its degenerate plane perpendicular to B. Thus its negative-dominant eigenspace, when present, always has nonzero field projection, justifying the simplified condition.

[proved, dimension count] In a simple-spectrum principal frame, T_s has 3 exchange eigenvalues and 3 DM components, minus scale, giving 5. K_D has 3 eigenvalues, 1 axial DM component, and 2 transverse field components, minus scale, giving 5. K_0 has 3 eigenvalues and 2 field components, minus scale, giving 4. N_m has 3+3+3-1=8. N_a has 3+1+3-1=6. N_0 has 3+3-1=5. Discrete axis choices and signed permutations do not change dimension. For T_d, the residual rotation in the double eigenspace fixes the direction of D_perp, leaving a,c,|D_perp|,D_parallel, minus scale, giving 3.

[proved, degenerate cases] Scalar S has no time-reversal-only points at B=0, since D leaves an axial proper symmetry. This proves that T_s and T_d exhaust T. The other rows use invariant conditions and therefore already include repeated eigenvalues and rank-deficient J. For example a double eigenvalue lowers the generic dimension of N_m from 8 to 6, while K_0 then has maximal dimension 3. The maximum dimensions displayed in the table remain valid after the indicated removals: there are open parameter choices in each top-dimensional row for which its inequality fails to certify a gap. Removing finitely many measured points does not alter these dimensions.

[proved, useful diagnostic] Aligned DM in a simple-spectrum frame is codimension two in the full symmetry-free family. D=0 is codimension three. Thus the direct external gap criterion does not reduce the generic eight-dimensional family to a hypersurface or to isolated junctions. Section 4 supplies an additional lower-dimensional exact region; it does not change that conclusion. In particular, a generic misaligned-D perturbation is not covered by the direct theorem even if it is physically still inside the same gapped phase.

### 5.2 What is not specified by (7)

[proved, logical limitation] Let C_reg denote the union of the controlled continuum and stable gapped neighborhoods actually established with all their hypotheses, and let Z_known denote any independently settled special sets. The physical unresolved set would be

    R_phys = R_alg \ (C_reg union Z_known),

with the empirical evidence recorded separately. C15-C17 do not give C_reg as inequalities in J,B, or a uniform algorithm that decides membership. Consequently neither the exact boundary nor the dimension of R_phys, and certainly not of its conformal subset, has been established. It could include finite-c multicritical points away from every known singular parent. The table is the exact algebraic envelope that the present certificates justify, not a substitute for that missing information.

[proved] Numerical evidence is not a parameter-space subtraction with a certified radius. At a point with finite-size evidence for a gap or c<=1, one still needs scale control to exclude a much longer crossover. At the physics standard, converged calculations can be persuasive evidence for that point and tracked neighborhoods. They do not establish a quantified exclusion across an entire five- or eight-parameter stratum. This is why an exact remaining critical set with assigned codimensions cannot honestly be reported here.

## 6. What the junction and two-fluid hunts establish

### 6.1 Ising lines and liquid regions

[numerical, inherited, not rerun] `../../contrib/d2_multicritical/README.md` reports rings through L=20 and VUMPS through chi=64 for the axial interacting-Kitaev slices. At delta=0.5, V=6, h=2.0, the liquid gives c_fit=0.997 with correlation lengths 43 to 164. At V=6,h=3.7, chi=32 to 64 gives a short-correlation-length state, xi=6.7 to 6.9, while the Ising line is near h=3.973. The V=8,h=4.8 check is similarly gapped with xi about 7.8. These data support an intervening gapped strip at the sampled cuts. They do not measure an independent Ising sector inside the liquid.

[proved, inherited free limit; numerical, interacting continuation] At zero interaction the chiral-Ising collision has a cubic slow branch, not an additional finite-velocity sector. The interacting endpoint tracking uses L<=20 and identifies the same first soft momentum approaching zero. That is evidence for a Lifshitz continuation, not a proof for all interactions. The exact quadratic and zero-field XYZ exclusions remain stronger on their respective loci.

[numerical, inherited limitation] The tilted-field XXZ map at anisotropy 2 uses only L=12,16, with h_x<=2 and h<=8. Its largest chord estimate is 1.21. The Neel Ising line and its junction were not resolved. This is a two-dimensional sampled slice inside K_0, itself a four-dimensional part of the full reflection-only branch. The predicted commensurate-incommensurate junction is therefore [heuristic], not an observed exclusion. The quoted hunt contains no comparable TR-only scan and does not cover N_m. It also does not scan the generic nonzero-D reflection branch K_D.

[assumed, controlled mechanism] In a one-species dilute onset, K tends to one. A nonoscillating pair source cos(2 theta) has dimension 1/K and is relevant there, whereas its irrelevance in a floating phase requires K<1/2 in this normalization. This explains why the ordinary band-edge Ising boundary and that floating regime are not independent degrees of freedom. It applies only after identifying the same dilute species and the permitted vertex. A different neutral Ising sector is not eliminated by this argument.

### 6.2 An exact variational improvement to the two-fluid screen

[proved, new optimization of the supplied trial energy] In the exact-vacuum slice of `../../contrib/d2_two_fluid/README.md`, take hopping amplitude tau=exp(i phi), interaction V, chemical potential mu, and conditional creation g_L,g_R, with no direct vacuum source or pairing. Minimizing the phase of a uniform product state of density n gives the energy relative to the empty vacuum

    e(n)=(-2 cos(phi)-mu)n + A n^2
                     -2G n^(3/2) sqrt(1-n),
    A=V+2 cos(phi),      G=|g_L+g_R|.

This expectation is exact for the trial state. For n>0 divide by n. To minimize A n-2G sqrt(n(1-n)), put n=(1-cos t)/2. The expression becomes

    A/2 - (A/2) cos t - G sin t.

Its minimum over 0<=t<=pi is (A-sqrt(A^2+4G^2))/2. Thus the exact uniform-product threshold is

    mu_prod = -2 cos(phi) + (A-sqrt(A^2+4G^2))/2.           (8)

For G>0 the minimizing density at threshold is

    n_* = [1-A/sqrt(A^2+4G^2)]/2.

For mu>mu_prod some product state has energy below the vacuum, so the exact vacuum cannot be the ground state. Equation (8) is an upper bound on the actual vacuum-stability threshold, not a lower bound on any excited band. For G=0, it reduces to min(-2 cos(phi),V), with an endpoint density; when also A=0, all densities have the same energy per particle at threshold. The phase convention is checked directly with a two-site matrix in `checks.py`.

[proved, limits of this witness] A variational state below the vacuum does not exclude a competing dilute state below both, or a finite-density two-fluid phase elsewhere. Comparing (8) with a rigidly filled, number-truncated two-body band cannot promote that band onset to an exact threshold. When g_L+g_R=0 the uniform witness loses the conditional-creation gain; nonuniform product states may still give a stronger bound. Direct vacuum sources and pairing also invalidate this exact-vacuum setup. These are concrete remaining pieces, not already closed by the hunt.

[numerical, inherited] The hunt sampled 10,100 general exact-vacuum bonds and 1,382 near the antisymmetric-creation slice. Its two-body problem is exact only with U(1); otherwise it is a truncation. Full-ring screening used L=12. Product trials used periods through six on selected survivors. Candidate C at mu=-2.30 gives c_fit=1.01 with a seven-site VUMPS cell and chi up to 48. Candidate S3 gives c_fit=1.03 at chi up to 48. The additional xi-versus-chi slopes are consistent with one mode. S1 had no usable estimate and S2 remained unconverged at chi<=16. Rings for C at L=14 to 20 had strongly drifting momenta and did not determine c or z. Thus the evidence is useful against those particular two-fluid constructions, but there is no numerical exclusion of the whole exact-vacuum slice, much less the dressed-vacuum generic family.

[proved, dimension of that slice] A specified uniform product eigenvector imposes four real conditions: two adjacent-pair amplitudes and two single-flip amplitudes. Allowing its orientation adds two variables. From 12 physical J,B coefficients this gives a generic incidence dimension 10 before quotienting, hence dimension 6 after rotations and scale. Equivalently, in its own z frame the surviving bond has equal symmetric xx and yy coefficients, antisymmetric xy, arbitrary four mixed xz/yz coefficients, zz, and B_z, with B_x,B_y fixed by single-flip cancellation; this is eight real coefficients, minus the residual z rotation and scale. Therefore the hunt's exact-vacuum slice is already lower-dimensional than N_m. Its complement is not addressed by replacing the vacuum with an uncontrolled dressed version.

[numerical, new finite check] One hundred random two-site bonds verify the product expectation and threshold (8), with maximum errors 8.9e-16 and 1.8e-15. This checks a variational identity, not a thermodynamic central charge.

## 7. Pushing the continuum argument into T, K, and N

### 7.1 What antiunitarity buys

[proved, finite-size symmetry statement] A uniform antiunitary commutes with translation as an operator. If a translation eigenvector has eigenvalue exp(ik), its antiunitary image has eigenvalue exp(-ik) and the same energy. Thus in an antiunitary-invariant vacuum the excitation spectrum is paired under k -> -k. On odd rings physical time reversal also squares to -1, giving Kramers degeneracy. Neither statement counts independent conformal sectors. If the vacuum breaks the antiunitary, the relation may instead connect excitations above different vacua.

[assumed, conditional dilute-field exclusion] Suppose two hypothetical soft minima at +/-Q are described by local bosonic annihilation fields psi_+,psi_-, exchanged by the antiunitary with A psi_+ A^(-1)=psi_- and vice versa. Their product V=psi_+psi_- has zero momentum. Antiunitarity sends gV+g*V^dagger to g*V+gV^dagger, so one real pairing coefficient is allowed, not forced to vanish. At independent dilute repulsive stiffnesses K_+=K_-=1,

    dim(V)=1/(4K_+)+1/(4K_-)=1/2.

A nonzero such coupling is relevant. This excludes a stable two-fluid continuation under these premises, and explains why opposite-momentum valleys alone are especially weak evidence in T and K. The same conclusion at a non-diagonal Gaussian point requires the full covariance matrix.

[proved, limitation of the inference] The conclusion fails as a general exclusion if the local product is forbidden by statistics or the actual symmetry representation, if its coefficient is tuned to zero, if the stiffness makes it irrelevant, or if the presumed fields do not constitute the complete low-energy theory. A Kramers spinon need not transform like the local boson assumed above. In a symmetry-free vacuum the two centers need not even be opposite, so the displayed pairing can oscillate. This is a mechanism test, not a proof on the five-dimensional antiunitary strata.

[assumed, inherited anomaly input] Time reversal with one Kramers doublet per unit cell gives an LSM-type obstruction to a unique symmetric gapped ground state. The ordered doublet in Section 4 satisfies that obstruction by breaking time reversal. A gapless anomaly-matching sector can coexist with an anomaly-free neutral sector, so the anomaly gives no upper bound on total c. Plane reflection has no Kramers-per-cell obstruction of this kind. The symmetry-free family has no corresponding internal-symmetry anomaly.

### 7.2 The neutral Ising obstruction remains concrete

[assumed, candidate field theory] Suppose a calibrated liquid has an additional independent local Ising sector with finite velocity. Ising times one boson has c=3/2. The Ising energy epsilon has dimension one and its coefficient must vanish. Its order field sigma has dimension 1/8. If sigma is local, uniform, and even under all exact symmetries, its coefficient also has to vanish. If sigma is odd under the antiunitary or has translation momentum pi, that particular magnetic perturbation is forbidden.

[proved, conditional operator analysis] Removing the uniform magnetic perturbation does not remove every mixed term. A local sigma times an allowed dimension-one boson current has dimension 9/8 and is relevant. A current odd under the same antiunitary can make this product even. With translation momentum pi, a boson harmonic of momentum pi can do the same. The terms epsilon times a current have dimension two and require marginal-flow analysis with the actual velocities. Nonderivative boson vertices and their products with sigma or epsilon must be tested for locality, momentum, symmetry and dimension. These are algebraic consequences of the chosen product operator spectrum, not claims about a realized qubit fixed point.

[proved, codimension statement] Let p be coordinates on one of the bare strata of dimension d, and g_alpha(p) the independent relevant coefficients in an identified continuum matching. If their simultaneous zero exists and has locally constant Jacobian rank r, the constant-rank theorem gives local dimension d-r. Without a known zero and its rank, there is no established microscopic codimension. For example, a protected Ising mass alone would require one tuning and could leave dimension four in a five-dimensional branch, or seven in the eight-dimensional branch, if the remaining mixed couplings were harmless. An unprotected uniform Ising field adds at least one further constraint. These are conditional counts. Neither number forces emptiness.

[heuristic, concrete candidate test] At any extra gap closing, first identify the translation momentum and antiunitary action of its putative order field using spin, bond-energy, chirality and longer local-string matrix elements. Perturb in independent tangent directions of J,B and measure their low-energy matrix elements. An integrated local field of scaling dimension Delta has leading matrix elements proportional to L^(1-Delta), which allows estimation of the independent response directions once the tower is identified. Tune the Ising mass and every permitted relevant mixed coefficient. Then establish a second finite velocity and a common scaling window. A peak in a second correlation channel can be a harmonic of the first boson and does not pass this test.

### 7.3 Two fluids, cells, and non-product theories

[assumed, candidate Gaussian theory] For two independent bosons without density-current mixing, the vertex exp[i(l.phi+m.theta)] has dimension

    Delta = (l^T K l + m^T K^(-1) m)/4,

in the C16 convention. One must first identify the physical local charge lattice. At a fixed positive K, only finitely many vectors have Delta<=2, so the resonance test is finite for an identified candidate. There is no uniform finite cutoff near singular stiffnesses, and generic chiral covariance requires a more general quadratic form. C15's constrained gas has K=(1-rho)^2 approaching 1/4, so a universal K>=1/2 argument cannot close this test.

[heuristic, concrete test] A multi-fluid candidate needs two positive stiffness and susceptibility eigenvalues, two independently dispersing towers, stable energy convexity, and absence or irrelevance of every local resonant conversion. Without exact microscopic charges, arbitrary twists about two spin axes do not measure two emergent fluid stiffnesses. The densities and currents must first be identified in the low-energy operator algebra. Few-particle truncations must be checked for leakage into the omitted sectors, especially in N_m.

[proved, inherited limitation] The finite frustration-free witness bounds only the length of an exact zero-energy segment. It does not classify perturbations of an exponentially degenerate parent. C16's blockade projection contains a complex hopping, a chemical potential, and a complex constrained source. Second order generates further interactions. No symmetry or bare Fourier-degree argument here limits the number of neutral modes in that projected dynamics or in a frustrated quantum-selected cell.

[heuristic, concrete test] Near a singular parent, derive the retained space and every leading relevant virtual term, then compare its towers with the finite-penalty microscopic chain. For a spontaneous cell, first establish its period, wall gap, internal multiplets and endpoint matrix elements. Test compatible and incompatible ring lengths. Separate the fixed cat-state entropy from propagating modes. For a candidate non-product CFT, identify towers and local operator products directly instead of forcing an Ising-plus-boson fit.

[proved, logical obstruction] Excluding all Ising-plus-boson products would still not exhaust unitary CFTs with c>1. A single interacting CFT need not split into independent stress tensors or a recognized parton coset. Thus a claim that there is no second visible spin-wave branch is weaker than the desired result. The final dynamical statement has to control the complete infrared theory, including neutral and non-product possibilities.

## 8. Numerical protocol, checks, and remaining work

### 8.1 What was computed here

[numerical] Run from this folder:

    python3 checks.py

The script uses Python, NumPy and SciPy, with one BLAS thread. It writes `checks_results.json`; `checks_output.txt` records the printed run. The complete run took less than one second locally. It performs no network calls and reads no external data.

[numerical] The checks are:

| Check | Size or number | Outcome |
|---|---|---|
| Symmetry sign actions on simple-spectrum examples in T, K and N | Seven 3x3 parameter examples | Expected proper/improper actions and certificate failures |
| External-theorem certificate against complete spectra | Three bonds, L=4,6,8 | All nine gaps exceed their certified lower bounds |
| Averaging lemma with genuinely rotating axes | Six sites, J=-R_z(2pi/6), degenerate singular values | Bond-equation error 5.2e-16; rotating lower bound 0.852154, uniform lower bound 1.4 |
| Transverse singular-value formula (2) | 200 real 2x2 blocks | Maximum error 1.4e-15 |
| Uniform TR product eigenvector | L=6 | Zero leakage, consistent with the algebraic condition |
| New full domain-wall operator inequality | Six bonds, L=4,6,8, including kappa=0,b!=0 | All eighteen inequalities, ground energies and gap bounds pass within roundoff |
| Uniform-product two-fluid witness | 100 two-site bonds | Expectation and minimization errors below 1.8e-15 |

[numerical, scope] These are algebraic and finite-spectral checks. No unknown central charge was estimated here, so there is no new thermodynamic exclusion based on small rings. The finite-dimensional inequalities are proved in the text for all L, independently of these checks. No additional random phase scan or DMRG run was undertaken.

### 8.2 A test that would distinguish progress from another false positive

[assumed, standard CFT diagnostics] For a genuine candidate with independent nonchiral sectors,

    S_L(ell) = (sum_a c_a/3) log[(L/pi)sin(pi ell/L)] + const + ...,
    E0(L) = e_inf L - (pi/6L) sum_a c_a v_a + ... .

Each v_a must be measured from its own tower. With left-right drift, replace c_a v_a by c_a(v_a,R+v_a,L)/2. Twisted or nonvacuum sectors have additional 1/L terms. First verify z=1 and positive limiting velocities, then compare entanglement with the velocity-resolved spectrum. The existing free chiral controls give ring entropy estimates as high as 2.37 at L=20 despite exact c=1, so a large short-ring entropy slope alone is not discriminating.

[heuristic, prioritized work] The next calculations should be organized as follows.

1. Complete the unresolved Neel Ising line in K_0, using larger compatible rings and converged finite-entanglement data, then continue into K_D. The present L<=16 map does not locate its junction. This is a bounded numerical task with an existing calibration, although it will not settle all of K.
2. Use (3)-(6) as an analytic anchor in T_d, and track all soft wall and composite branches as kappa decreases. Then break the double eigenvalue of S to enter T_s. This starts from a proved gapped phase and tests where the one-mode description first fails. The sufficient boundary kappa=0 must not be mistaken for the actual gap-closing surface.
3. Finish S1 and S2 from the two-fluid hunt before calling that tested set exhausted. Use (8) for a cheap exact vacuum witness, retain nonuniform trial states when g_L+g_R vanishes, and distinguish full-ring results from rigid two-body filling. Dressed-vacuum finite-density candidates still require a separate calculation.
4. In N_m, search for an independent neutral tower along a calibrated c=1 region and determine the relevant-coupling response matrix. This addresses the main multicritical obstruction more directly than another broad entropy-only scan.
5. For a global proof, establish an exhaustion principle for the complete scaling theories on the remainder, including singular neighborhoods and isolated non-product fixed points. Additional negative samples cannot supply this step.

[assumed, local physics extension of the new result] Standard stability of the ordered easy-axis phase would extend the strict b=0 wall-bound points to relatively open gapped regions in T. Standard stability of a unique gapped product phase would extend the b!=0 points into open symmetry-free neighborhoods. These are useful local physics consequences when the stability input is supplied. This run gives no explicit radii for such neighborhoods and does not use them to claim an exact smaller remainder.

### 8.3 Verdict

[proved] The direct gap theorem has now been exhausted under arbitrary product spin rotations, and a distinct domain-wall estimate settles additional genuine members of all three requested symmetry classes. The explicit search envelope is (7) and the seven-row table. Its maximal dimensions remain 5,5,8. The stronger claim that all remaining conformal candidates lie on known lower-dimensional junctions has not been established.

[heuristic] The mechanisms and numerical evidence favor a single-fluid account of the familiar critical regions, but a full physics-level proof is not one final finite-size test away. The last essential step would be a microscopic argument that every remaining critical point admits a complete continuum description with c<=1, or a flow from a complete c<=1 continuum theory with all competing soft sectors controlled. This must include the still-open unitary-Z2 cases, the singular-parent neighborhoods, and possible isolated strongly interacting fixed points. Calling this merely the absence of another boson would omit non-product CFTs.

[proved, completion record] The requested sources were read, the finite-dimensional results above were derived and checked, and local progress and memory files were updated. No parent file, paper, repository status file, or external service was modified. The project-wide STATE.md could not be updated without violating the instruction to write only inside this folder; `STATE.md` here records the handoff. The work is unrefereed.

Sources used locally: `../../specs/_d2_background.md`, `../../STATUS.md`, `../work_C15/RESULT.md`, `../work_C16/RESULT.md`, `../work_C17/RESULT.md`, the four notes in `../../contrib/d2_gap_reduction/`, `../../contrib/d2_multicritical/README.md`, `../../contrib/d2_two_fluid/README.md`, and the symmetry formulas in `../work_C6_d2/c6.tex`. The project memory was read before the research files. Source central-charge conclusions and finite-size evidence have not been promoted to new proofs.
