# C76: simple mechanisms outside route A

Best admissible result found: c-(d-1)=0. No counterexample was established for 3<=d<=12. The maximum c/(d-1) in the surveyed Hermitian conformal families is 1. The largest c/d in this range is 11/12, attained by SU(12)_1 and by the eleven-species supersymmetric t-J liquid with one hole state.

The useful conclusions are these. Spin-charge separation does not turn q gapless species densities into q+1 sectors. Pure single-occupancy hopping actually freezes the particle order, rather than supplying q coherent liquids. Interactions can exceed C1's bound c<=log_2 d, very economically at d=4 with c=3, but the explicit examples still saturate d-1. Doped higher-spin integrable chains really can acquire neutral sectors, yet their local-state accounting makes them inefficient.

Status labels used below are “literature” for continuum identifications read in primary sources, “proved” for deductions in PROOF.md, “numerical” for this folder's finite checks, and “conditional” for claims requiring a complete Gaussian infrared description. This is a survey of the specified mechanisms, not a classification of all nearest-neighbor bonds or all integrable representations.

## Survey and local-dimension accounting

| Family and regime | Local d | Gapless c | c-(d-1) | Status |
| :-- | --: | :-- | :-- | :-- |
| Fundamental AF SU(d) permutation, all colors present | 3 through 12 | d-1 | 0 | Literature, established reference family |
| Supersymmetric t-J, q flavors and one hole, interior compressible phase | q+1, q=2 through 11 | 1+(q-1)=q | 0 | Literature, ordinary Hermitian graded permutation |
| Repulsive spinful Hubbard, away from half filling and saturation | 4 | 1 charge + 1 spin = 2 | -1 | Literature, Bethe ansatz |
| Attractive Hubbard, paired balanced phase | 4 | 1 | -2 | Spin gap, charge liquid |
| Attractive Hubbard, partially polarized coexistence phase | 4 | 2 | -1 | Pair and excess-particle liquids, not an additional spin copy |
| Ordinary repulsive q-component Fermi-Hubbard liquid | 2^q, hence 4 or 8 for q>=2 in range | q | q-2^q+1 = -1 or -4 | Literature/RG; ordinary q>2 lattice Hubbard is not assumed integrable |
| Projected degenerate Hubbard, at most two fermions per site | 1+q+C(q,2)=4,7,11 for q=2,3,4 | q | -1,-3,-6 | Literature Bethe-ansatz identification, not independently re-proved |
| Takhtajan-Babujian spin-S chain, S=1,3/2,...,11/2 | 2S+1 | 3S/(S+1) | 3S/(S+1)-2S <0 | Literature; strings are not independently countable bosons |
| Integrable doped spin-S chain, S=1/2,1,...,5/2 | 4S+1=3,5,7,9,11 | 5-6/(2S+1)=2,3,3.5,3.8,4 | 0,-1,-2.5,-4.2,-6 | Literature low-energy/TBA count; neutral-sector decomposition has qualifications below |
| Rational AF SO(d) vector chain, d>=3 | d | d/2 | 1-d/2 | Literature SO(d)_1; SO(3) is SU(2)_2 |
| AF Sp(2r) fundamental chain | 2r=4,6,8,10,12 | r in Martins's multivelocity identification | 1-r | Published analytical/numerical identification, not a generic Sp(2r)_1 assertion |
| Strongly interacting Fock-space ordinary swap | 4 or 8 | 3 or 7 | 0 | Same SU(d) reference, explicit fermion bond below |
| Pure hard-core q-species hopping with equal amplitudes | q+1 | Charge c=1 plus frozen color degeneracy | Not a complete CFT | Proved order conservation, numerical full-spectrum check |
| Mixed-grading rational OSp(n|2m), n,m>0 | n+2m | Do not enter an admissible c | Excluded | Standard bond fails ordinary Hermiticity |

For each integer d in the requested range, the saturating benchmark is

| d | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
| --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| c | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| c/d | .666667 | .75 | .8 | .833333 | .857143 | .875 | .888889 | .9 | .909091 | .916667 |

This maximum is over the surveyed admissible examples, not a proved maximum over all bonds.

## Explicit Hamiltonians and mechanisms

Write X^{ab}=|a><b| and n=sum_{a>0}X^{aa}. The basic hard-core mixture is

    H = sum_i [-sum_a t_a(X_i^{a0}X_{i+1}^{0a}+h.c.)
               +sum_{a,b>0} J_ab X_i^{ab}X_{i+1}^{ba}
               +sum_{a,b>=0} V_ab X_i^{aa}X_{i+1}^{bb}]
        -sum_{i,a>0} μ_a X_i^{aa},

with J_ba=J_ab* and real V_ab. This is a literal one-site invariant two-site model on d=q+1 states.

At equal t_a=t and J_ab=t, choose the diagonal term -t n_i n_{i+1}. This gives the supersymmetric generalized t-J bond

    h=-t T+t(P_occ-n⊗n),
    T=sum_a(X^{a0}⊗X^{0a}+h.c.),
    P_occ=sum_{a,b>0}X^{ab}⊗X^{ba}.

The interior phase has c_charge=1 and c_flavor=q-1. q=2 gives d=3,c=2. The spin-orbital q=4 example gives d=5,c=4, with the thermal coefficient checked against the two distinct velocities in [Sirker, Damerau and Klümper, Sec. II](https://www.jensdamerau.de/Su41_v2.pdf). The ordinary permutation P_d=sum_{a,b}X^{ab}⊗X^{ba} provides the complementary ungraded saturator. Its fundamental SU(d)_1 identification is also the project benchmark; [Sutherland's nearest-neighbor multicomponent model](https://doi.org/10.1103/PhysRevB.12.3795) is the original integrable construction.

The Hubbard representative is

    H=-sum_{i,a}t_a(c†_{i,a}c_{i+1,a}+h.c.)
      +sum_{i,a<b}U_ab n_{i,a}n_{i,b}-sum_{i,a}μ_a n_{i,a}.

Here the actual Fock dimension is 2^q, not q+1. For q=2, repulsive generic filling has two modes, as resolved by the dressed charge matrix in [Frahm and Korepin](https://doi.org/10.1103/PhysRevB.42.10553). A q-component repulsive metal similarly has q bosons when its charge and flavor sectors remain gapless; see the bosonization and numerical analysis of [Szirmai, Legeza and Sólyom](https://arxiv.org/abs/0708.3799). This is saturation of C1, not an excess. At a Mott transition the charge mode can be removed. An attractive balanced gas can lose flavor modes to pairing or larger bound complexes. The attractive partially polarized case has a bound-pair fluid and an excess-fermion fluid, as derived in [Cheng et al.](https://arxiv.org/abs/1708.07776). Counting a bound complex and all its constituent liquids simultaneously is unjustified.

A different model projects the same hopping onto at most two occupants per site. With P_<=2 the projector, its bond is the projected hopping plus the same on-site repulsion. Its d is 1+q+C(q,2). The Bethe-ansatz analysis in [Frahm and Schadschneider, Eqs. (1), (17)-(21)](https://arxiv.org/abs/cond-mat/9207026) identifies q gapless sectors. Thus the interesting small cases are (d,c)=(7,3),(11,4). It must not be confused with ordinary SU(q) Hubbard on 2^q states or with the one-occupant t-J chain.

For spin S let Π_j project the bond onto total spin j, and H_j=sum_{r=1}^j 1/r, H_0=0. An explicit TB bond, up to scale and constant, is

    h_TB=sum_{j=0}^{2S} H_j Π_j.

It realizes SU(2)_{2S}, c=3S/(S+1). At S=1 this is proportional to S_i·S_j-(S_i·S_j)^2. The neutral parafermion piece in SU(2)_{2S}/U(1) contributes c=2(2S-1)/(2S+2), already included in the quoted c. Counting Bethe string types instead would overcount. [Alcaraz and Martins, finite-size analysis](https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Physics%20A%20Mathematical%20and%20Theoretical/1988%20Volume%2021/0305-4470_21_23_021.pdf). The spin-S realization of the SU(2S+1) permutation can have 2S occupied Bethe seas under Cartan fields, but gives only c=2S=d-1. Several seas in staggered chains or models deformed by higher conserved charges require a range audit; higher charges generally introduce three-site or longer interactions and do not furnish a nearest-neighbor exception.

The particularly relevant neutral-sector example is the doped spin-1 chain. Each site is V_1⊕V_{1/2}, d=5. With s_i the local spin magnitude, x=S_i·S_j, and P exchanging the full local states, its explicit bond is

    h = (1/2)[x/(s_i s_j)-1+δ_{s_i s_j,1}(1-x²)]
        -(1-δ_{s_i,s_j})P x.

It has c=1+3/2+1/2=3 in the small-doping field-theory description. The added singlet Majorana is real physics, but the hole carries a doublet and costs two local states. [Frahm, Pfannmüller and Tsvelik, Eq. (1)](https://arxiv.org/abs/cond-mat/9803145).

For the higher-spin generalization, d=(2S+1)+2S=4S+1. At small positive doping the continuum content is a charge boson, SU(2)_{2S}, and the minimal model M_{2S+1}, giving

    c=1+6S/(2S+2)+1-6/[(2S+1)(2S+2)]
     =5-6/(2S+1).

The microscopic fused nearest-neighbor construction and TBA sum rule are in [Frahm, Secs. 3 and 6](https://arxiv.org/abs/cond-mat/9904157). At generic doping the two magnetic thermal coefficients vary with doping while their sum remains fixed; they should not individually be called fixed minimal-model central charges. The small-doping decomposition is the controlled interpretation used here. Even the full sum never reaches d-1 except for S=1/2.

For the orthogonal family a convenient bond is

    h_SO=P-2K/(d-2), K=|Ω><Ω|, |Ω>=sum_a|aa>.

It is Hermitian and has c=d/2. In contrast the standard mixed OSp bond uses the graded contraction E rather than the positive rank-one K. The explicit Frobenius norm calculation in PROOF.md shows its non-Hermiticity for every n,m>0 in the nonsingular rational family. Published sector-dependent effective charges, for example z/2 or z-1 with z=n-2m in [Frahm and Martins, Eq. (2.6)](https://arxiv.org/abs/2202.13405), are not counts of positive-norm liquids. The d=3 OSp(1|2) example is explicitly identified as non-Hermitian in [Tavares and Ribeiro](https://doi.org/10.1016/j.nuclphysb.2017.05.020).

The pure Sp(2r) AF bond is P+K_J/(r+1), where |J>=sum_{a,b}J_ab|ab> and K_J=|J><J| for the antisymmetric invariant form. This positive-K sign follows by differentiating R(u)=u(u+r+1)I+(u+r+1)P-uK_J and using PK_J=-K_J. The Sp(4) R-matrix and transfer construction are given in [Li et al., Sec. 2](https://arxiv.org/abs/1812.03618). [Martins](https://arxiv.org/abs/cond-mat/0012159) identifies r elementary sectors, c=r, whereas the other r-1 apparent branches are composites. For Sp(4), using the slower velocity in a one-velocity Casimir fit gives a coefficient approaching 4 although the proposed total c is 2. This is a useful concrete warning. The published c=r interpretation is reported at its stated analytical/numerical status; even the alternative Sp(2r)_1 value r(2r+1)/(r+2) would remain below 2r-1 for r>=2.

General SU(N|M) notation alone is also insufficient. The Hermitian rational bond is ±P_g. In the -P_g convention, the all-even-color subspace is a ferromagnetic permutation chain and can have quadratic color waves. Indeed the bare -P_g chain has the frustration-free ground energy -L whenever an even color is available, and its one-particle threshold has z=2. A chemical potential is essential for selecting the finite-density conformal phase. Thus it is false that every grading, sign, and chemical potential automatically gives c=N+M-1. The conformal one-even-hole, q-odd-particle t-J subfamily above does saturate. Extra ferromagnetic multiplets or nested root levels are not additional relativistic sectors. No all-gradings classification is asserted here.

Haldane-Shastry and other inverse-square chains were excluded. As a further screen, the Izergin-Korepin noncompact regime has the formal black-hole charge c_BH=2+6/(k-2), but normalizable continuum states start above zero and its untwisted effective charge is 2. [Vernier, Jacobsen and Saleur, Eqs. (77)-(78)](https://arxiv.org/abs/1404.4497). The relevant standard vertex-chain realization is among the non-Hermitian chains discussed in [their later analysis](https://inspirehep.net/files/a308db51ed5a23e757b47d719fad53f3). It is not an admissible Hermitian d=3 counterexample. Formal c, finite-size effective c, and unitary vacuum-sector c must be distinguished.

## Interacting bands can beat C1, but not the target in this survey

A particularly simple answer to part (ii) is

    h=P_4=(2S_i·S_j+1/2)(2T_i·T_j+1/2), d=4,c=3.

Identify the four states with |0>,|a>,|b>,|ab> of two spinless fermionic orbitals. PROOF.md gives a literal nearest-neighbor canonical-fermion polynomial and audit.py verifies it equals P_4. It is parity-even and conserves both orbital charges. It includes correlated exchanges beyond ordinary quartic Hubbard interactions. This exceeds the quadratic bound 2, with no new local states. The analogous three-orbital configuration swap has d=8,c=7>3. In both cases c-(d-1)=0. These are strong-coupling examples; no claim is made that infinitesimal interactions create an additional boson around an ordinary free-band metal.

For unprojected spinful two-band models, local d=16 already lies outside the requested range. Counting only the one-electron states while retaining a charge liquid is a dimension error. Projection to one occupant changes both the space and dynamics, returning to the t-J accounting above.

## Answer to the q-species question

For q coherent single-sea density liquids, c=q. Spin-charge separation changes the basis to one total density and q-1 relative densities. The hole density is constrained by n_0=1-sum_a n_a. Species-dependent hopping changes velocities and stiffnesses, not this count. [Proved within the complete Gaussian hypothesis; standard RG stability at generic incommensurate densities.]

With only vacancy hopping and diagonal potentials there is a stronger microscopic obstruction to the proposed picture: particles cannot pass, and the ordered species word is conserved. Equal hopping gives one mobile spinless liquid plus static flavor degeneracy. Unequal hopping does not remove the order conservation. Explicit color exchange is needed for the usual coherent flavor liquid.

Therefore the proposed mechanism does not produce c=q+1. I did not prove that every two-site hard-core color Hamiltonian forbids c>=q+1. Such a theorem would have to exclude an extra neutral critical sector or independently stable emergent pockets. None was found here. The doped-spin example illustrates why charge-rank counting alone cannot settle that stronger question, while also showing the actual dimension cost in the known integrable realization.

## Small exact diagonalizations

All calculations use numpy/scipy, sparse symmetry sectors, and a fixed random seed. They are small-ring checks, not thermodynamic determinations for unfamiliar candidates. The largest symmetry block has 48,400 states. No cluster calculation was needed.

| Model | L | E_0 | Interval c fit, ell>=2 | Casimir c using reference v,e∞ |
| :-- | --: | --: | --: | --: |
| SU(3) swap | 6 | -4.605551275464 | 2.564803 | 2.113464 |
| SU(3) swap | 9 | -6.579435759086 | 2.325552 | 2.056080 |
| SU(3) swap | 12 | -8.624526804645 | 2.239152 | 2.035142 |
| SU(4) swap | 4 | -4.000000000000 | Too few cuts | 3.402078 |
| SU(4) swap | 8 | -6.921009740038 | 3.736949 | 3.113154 |
| q=2 t-J, J=2t, particle density 2/3 | 6 | -6.991015313210 | 2.586349 | Not extracted, two velocities |
| Same | 12 | -13.542600178410 | 2.245299 | Not extracted, two velocities |

The half-chain size estimates are c=1.982403 for SU(3), L=6→12; c=3.089079 for SU(4), L=4→8; and c=1.981850 for t-J, L=6→12. These calibrate the large upward bias in short-ring interval fits. In particular, SU(4) can appear to have c≈3.74 at L=8 even though it is a c=3 reference. The t-J L=9 ground level is degenerate; its arbitrary real eigenvector entropy is saved but excluded from the estimators.

A separate full-spectrum open-chain check at q=2, L=8, counts (holes,a,b)=(4,2,2), verifies six copies of the 70-state spinless spectrum, with maximum error 1.155e-14. Unequal hoppings (1,0.7) retain six disconnected color-word sectors. These directly test the frozen-order mechanism.

Hubbard spectra at U=0,4 and L=4,8,12 provide additional finite-bond checks. Shell effects, including free ground-state degeneracies at L=4,12 in the tensor-periodic convention, prevent combining those sizes into a reliable c estimate. They are not used to infer an excess. At the nondegenerate free L=8 point, all interval entropies match the independent correlation-matrix result to 2.23e-15. The interacting U=4,L=8 energy is -6.245317025012. This tests the implementation, not a new thermodynamic claim.

For Sp(4), the L=4 and L=8 Casimir coefficients obtained by dividing by the slow velocity are 4.549254 and 4.150369. They approach the published velocity-weighted value 4, not a central-charge violation. L=6 gives -1.903640 in the same subtraction and is not in that vacuum size sequence; it is explicitly excluded from an extrapolation. The positive-K bond passes the L=4 conserved-charge check ||[H,sum_i[h_i,h_{i+1}]]||=8.54e-15. Reversing only the K sign gives 82.62 and fails. This is a consistency test, not a proof of integrability.

The four mixed-OSp bond tests reproduce ||h-h†||_F=8sqrt(nm)/|2-n+2m|, with E²=(n-2m)E to machine precision. The canonical-fermion P_4 identity has zero numerical error. Eigenpair residuals in the main checks are at most 1.2e-12, including the excluded degenerate state.

Reproduce in this order:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 ed_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 extra_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 audit.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 symplectic_check.py

Files are PROOF.md, ed_checks.py, extra_checks.py, audit.py, ed_results.json, extra_results.json, audit_results.json, symplectic_check.py, and symplectic_results.json. No files outside work_C76 were modified. Project-level STATE.md was left untouched under the task's explicit confinement rule; this folder has its own STATE.md.

The remaining issue is narrow. A new small-d counterexample would need an additional propagating sector whose independence and stability are demonstrated microscopically. Adding the names “hole”, “charge”, or “spin” to the q density modes does not provide one.
