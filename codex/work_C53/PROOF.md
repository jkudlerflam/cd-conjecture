# C53. A theorem-covered logical chain at fixed alphabet

[proved, scope] There are two different modifications. A diagonal counterterm and field put the second-order reference in a published weak-coupling class, but leave C47's leakage problem. Replacing the port bond by a direct logical bond gives an exact application in the invariant PBP sector. The second modification keeps d=10, the original cell checks and one-site translation invariance, but uses four-site interactions. It does not solve the nearest-neighbor problem or Input S for the wrapper.

[read] The requested S background and C47 RESULT.md and PROOF.md were read in full. The collaborator's [demonstrator note](../../contrib/decisive_attack_20261028/demonstrator/README.md) and its selector audit were compared with the literal [cell builder](../../contrib/routeA_demonstrator/cell.py). Primary sources were checked through the browsing tool. The session's higher-priority source-verification instruction required this despite the background's offline instruction. All filesystem writes and Python calculations are confined to work_C53. No subagents were used.

[proved, conventions] A tag governs its paragraph and ensuing displays. Imported published results are labeled imported, not proved here. Units are per cell unless stated otherwise. The microscopic coupling is lambda, the freely chosen anisotropy is g, and the cell gap is Delta_cell. This avoids using Delta for two unrelated quantities. All parameters are fixed before volume and temperature limits.

## 1. Imported finite algebra and an explicit diagonal modification

[imported, proved in C47 Sections 1-2] The physical site is C^10 = P direct-sum B, with dim P=2 and dim B=8. B contains qubits sL,c,sR. In the PBP sector the cell qubits are pL,sL,c,sR,pR. At D=1, kappa=3/5, Jp=2, Jz=-1,

    A_j >= (2/5)(1-P_j),
    P_j = |Phi+><Phi+|_(pL,sL) tensor I_c
                     tensor |Phi+><Phi+|_(sR,pR).

The unperturbed cell has dimension 32, code dimension 2 and exact gap 2/5. Write V_j for the original port bond without lambda. C47 proves PVP=0 and

    -PVQA_Q^(-1)QVP = sum_j [a0 I + t T_j + u Z_j Z_(j+1)],
    T_j = sigma+_j sigma-_(j+1)+sigma-_j sigma+_(j+1),
    a0=-20505/11648,   t=-9/32,   u=243/11648.

Cross-interface terms vanish. Those results are used, not claimed as new derivations.

[proved] Set

    alpha=9/16=2|t|,
    q(g)=alpha g/4-u = 9g/64-243/11648,
    W(g,h)=q(g) sum_j Z_c,j Z_c,j+1
                         -(alpha h/2) sum_j Z_c,j,
    H_tune=A+lambda V+lambda^2 W(g,h).                 (1)

Here Z_c is the bare logical-qubit Pauli matrix, including identity on both partner qubits. On the physical alphabet the diagonal bond is the single, translation-invariant four-site rule

    q(g) lambda^2 Z_B,i nP,i+1 nP,i+2 Z_B,i+3,
    Z_B = 0 on P,   Z_B = I_sL tensor Z_c tensor I_sR on B.

Sum over every physical i. The field is the on-site rule

    -(alpha lambda^2 h/2) sum_i Z_B,i.

Thus the change uses no new letters. It has physical support length four, not nearest-neighbor support. Every term is diagonal in the letter basis, preserves the type sectors and commutes with the original total U(1) charge. The field breaks charge conjugation, as it must to move away from zero magnetization.

[proved] Since P_j is identity on c, W preserves the code and has the stated compression. Degenerate perturbation theory through degree two, or C45's local identity with F=W, gives

    H_log,2 = lambda^2 a0 M + s H_BM(g,h),
    s=alpha lambda^2,
    H_BM(g,h)=sum_j[-(Sx_j Sx_j+1+Sy_j Sy_j+1)
                            +g Sz_j Sz_j+1-h Sz_j],
    Sa = sigma_a/2.                                   (2)

Indeed SxSx+SySy=T/2, so the hopping is -alpha/2=t, the ZZ coefficient is u+q=alpha g/4, and the field is -alpha h Z/2. After the staggered transverse rotation the exchange is positive and the standard antiferromagnetic XXZ anisotropy is g. At g=0 the ZZ coefficient cancels without reducing the nonzero hopping. This improves on tuning Jp towards zero, which also reduces the hopping.

[proved] In Jordan-Wigner conventions n=(1+Z)/2, the same dimensionless reference is

    H_BM=-1/2 sum_j (a†_j a_j+1+h.c.)
          +g sum_j (n_j-1/2)(n_j+1-1/2)
          -h sum_j (n_j-1/2).                         (3)

The density interaction is instantaneous and nearest-neighbor. For the double-sum convention in the imported theorem,

    v(r)=1/2 if r=+1 or -1, and 0 otherwise,
    vhat(k)=cos k,    lambda_paper=g.

There is no factor-of-two ambiguity: the two ordered pairs in the double sum reproduce one bond of coefficient g.

[proved] At g=0 and h0=1/3 the one-particle energy is

    epsilon(k)=-cos k-1/3,
    k*=arccos(-1/3),    v*=sin k*=2sqrt(2)/3,
    rho*=k*/pi=0.6081734479693927... .                 (4)

There are two simple Fermi points. The quartic umklapp distance is

    dist(4k*,2pi Z)=4 arcsin(1/3)>0.

The bilinear distance dist(2k*,2pi Z) is also positive. The theorem below can fix these Fermi points by its chemical-potential counterterm. We do not assert that changing the field removes every possible higher-degree commensurability in a general action. The relevant quartic obstruction at half filling is absent here.

## 2. The remaining matching error for H_tune

[proved, using C47's local construction] W commutes with P_X on each of its cell supports, although it need not commute with A. Use C47's S1,S2 without change. Through order lambda^2,

    exp(lambda S1+lambda^2 S2) H_tune exp(-lambda S1-lambda^2 S2)
       = D+lambda^2 W+R_new,
    P(D+lambda^2 W)P=H_log,2.                        (5)

To make the new bound explicit, use C47's constants at spatial weight a and put

    w_a = 2 exp(2a)|q(g)| + exp(a) alpha |h|/2.

This bounds the incidence norm of W. C47 has sigma=|lambda|s1+lambda^2s2 <= 2s1|lambda| and 4sigma<=1/4. Its repeated-commutator estimate, also valid for range-two W, implies

    ||exp(ad_S)W-W||_a
      <= [(1-4sigma)^(-3)-1] w_a
      <= 24 sigma w_a.

For the last inequality expand (1-x)^(-3); the quotient by x increases on [0,1/4] and is less than 6 there. Consequently

    ||R_new||_a <= (C_a+48s1 w_a)|lambda|^3.          (6)

This is uniform in cell count. The logical scale s is fixed and nonzero at every fixed lambda>0.

[proved, limitation] Equation (6) does not give an exact Hamiltonian of the form (3). At third and later orders it still contains leakage and other interactions. Conditional massive integration gives small retarded kernels with higher logical degree. Neither the counterterm nor the field makes these kernels instantaneous density-density terms. A bound of order lambda^3 on free-energy density remains a constant when T tends to zero at fixed lambda. It does not determine a T^2 coefficient. Accordingly the published theorem in Section 4 applies to the reference (2), not to H_tune itself.

[numerical] The 64-dimensional isolated-interface ED retains the outer Bell pairs and independently reconstructs the virtual resolvent. At g=0,+1/50,-1/50 and h=1/3, the canonical low-space matrix minus (2) scales as lambda^3, with norm/lambda^3 tending towards approximately 3.30. At lambda=.0025 the measured ratios are 3.28773, 3.28774 and 3.28772. Cancellation of the degree-two ZZ term has not falsely been reported as cancellation of all interactions.

## 3. An exact finite-range redesign

[proved, definition] To obtain an actual theorem application, replace the port-port coupling, rather than only adding a diagonal term. Keep all original positive cell checks, the BB penalty K=2 and the B chemical potential mu=11/20. Delete lambda V on PP. Add the following four-site rule, with identities on the partner qubits understood:

    s nB,i nP,i+1 nP,i+2 nB,i+3
      [-1/2(sigma+_c,i sigma-_c,i+3+h.c.)
                                      +(g/4) Z_c,i Z_c,i+3],

and the on-site field -s h Z_B/2. Sum the rules over all physical sites. The factors nB at the endpoints mean that the displayed c operators act only on B. Call this H_exact. It has d=10, one-site translation invariance, range at most four physical sites and the original U(1). It is a different microscopic Hamiltonian. In particular, deleting lambda V is an order-lambda change, not a small order-lambda^3 correction to the original model.

[proved, exact matching theorem] In a frozen PBP sector, after subtracting -mu M,

    H_exact = A + s H_BM(g,h)_c.                     (7)

The suffix c means that the spin chain acts on the bare c qubits and as identity on every port and partner. Reorder factors as

    H_cell = H_syndrome tensor C^2_c,
    dim H_syndrome=16,
    P_j = |Omega><Omega|_syndrome tensor I_c.

This is a permutation of tensor factors, not a spectral embedding into a larger alphabet. Let G insert |Omega> in every syndrome factor. Since A G=0 and every c operator commutes with G's insertion,

    H_exact G=G s H_BM,
    G† H_exact G=s H_BM,
    [H_exact,P_j]=0 for every j.                     (8)

The commutator with A vanishes because A_j has the exact Bell code as a reducing subspace, and all other cells commute with P_j. The c interactions commute with each P_j directly. There is no Schrieffer-Wolff remainder and no retarded effective action in the all-code sector. The matching is an operator identity at every s, not a second-order approximation.

[proved, uniform excitation threshold] The inequality imported from C47 now gives the especially useful comparison

    H_exact >= I_syndrome tensor s H_BM
                              +(2/5) sum_j Q_j.      (9)

The two operators on the right commute. In a sector with r excited syndrome cells, the energy is at least E_log+(2/5)r, where E_log is the ground energy of s H_BM on the same ring and with the same boundary convention. The code attains E_log. Thus every frozen-sector ground state lies in the code, and all states below E_log+2/5 are logical states. This statement holds uniformly in M without requiring the entire extensive logical band to lie below the excited sector.

[proved, thermal matching] For any beta>0 the code contributes Z_log to the trace. The variational eigenvalue principle applied to (9), followed by factorization of the right-hand trace, gives

    Z_log <= Z_exact <= Z_log [1+15 exp(-2beta/5)]^M,
    0 <= f_log(T)-f_exact(T)
                   <= T log[1+15 exp(-2/(5T))].       (10)

The factor 15 counts the excited syndrome states. No operator monotonicity of the exponential is used. If X>=Y, their ordered eigenvalues satisfy E_n(X)>=E_n(Y), which implies Tr exp(-beta X)<=Tr exp(-beta Y). Both ground energies agree, so (10) transfers every existing low-temperature power coefficient of the logical chain. Unlike C47's comparison, it has no fixed additive lambda^3 error.

[proved, observable matching] Bare c operators preserve the code. In a code ground state, their real- or imaginary-time correlations under (7) equal those of s H_BM, since time evolution on the invariant code is exactly (8). This includes Sz_c, the logical bond current and their connected functions. The microscopic conserved charge reduces to sum Sz_c because each Bell pair is neutral. There is no source-dressing error in these identifications. Fermionic single-particle operators still carry the ordinary Jordan-Wigner string; they are not claimed to be single-site spin observables.

[proved, physical scope] These conclusions concern one invariant period-three sector and its infinite-volume state. We have not proved that this sector minimizes the unrestricted ten-letter Hamiltonian at nonzero s, or compared its Gibbs trace with every type pattern. The original zero-coupling selector still holds because H_exact agrees with the original model at lambda=0. Its extension requires a dressed defect argument. We do not infer that extension from a comparison of extensive run-energy lower bounds. In physical units f_site=f_cell/3 and v_site=3v_cell.

## 4. Imported theorem, with hypotheses and conventions exposed

[imported, primary source] Benfatto and Mastropietro, “Universality relations in non-solvable quantum spin chains,” Theorem 1.1, [arXiv:0907.2837](https://arxiv.org/html/0907.2837), published DOI [10.1007/s10955-010-9920-0](https://doi.org/10.1007/s10955-010-9920-0). Use their equations (1)-(3), (12)-(18), (21)-(23). The following is its statement in our symbols, with the gapless counterterm convention from the proof made explicit.

[imported, theorem BM] For

    H=-sum_j(Sx_j Sx_j+1+Sy_j Sy_j+1)-h sum_j Sz_j
                         +g sum_(j,k) v(j-k) Sz_j Sz_k,

assume real even exponentially decaying v, periodic Jordan-Wigner fermions, two interior simple Fermi points ±pF, and sufficiently small |g|. Fixing pF uses h=-cos pF-g-nu(g,pF) for vhat(0)=1, with the theorem's counterterm nu=O(g). There is analytic K(g) with

    X_plus=K,  X_minus=1/K,
    nu_bar=1/(2-1/K),  2eta=K+1/K-2,
    K=1-g[vhat(0)-vhat(2pF)]/(pi sin pF)+O(g^2),
    Omega(p0,p)=K/(pi vs) [vs^2 p^2/(p0^2+vs^2 p^2)]+R(p0,p),
    R continuous, R(0,0)=0,   susceptibility=K/(pi vs),
    vs=sin pF+O(g)>0.                                (11)

Here X_plus and X_minus govern the oscillating density and pairing correlations, eta is the fermion anomalous exponent, and nu_bar governs the anisotropy-induced inverse correlation-length scale. Omega is the Fourier transform of the connected Sz correlation in the zero-temperature thermodynamic limit. The smallness radius is existential, not numerically supplied.

[proved, parameter choice] Choose pF=k* from (4), v(r)=delta_(|r|,1)/2, and

    h(g)=1/3-g-nu(g,k*).                              (12)

Take |g| below the theorem's positive smallness radius and, if necessary, smaller so that |h(g)-1/3|<1/12 and vs>v*/2. Choose any fixed lambda>0 and set s=(9/16)lambda^2. This is an existential family of fixed Hamiltonians with no volume-dependent coupling. At g=0, every coefficient is explicit, h=1/3 and nu=0. No decimal nonzero g is claimed to be certified by an unspecified radius. The numerical g=±.02 controls use h=1/3 to check matching; they do not compute the RG counterterm (12).

[proved, hypothesis ledger] Each entry follows from the explicit Hamiltonian, not from an assumed universality extension.

| Hypothesis | H_exact in the frozen sector | H_tune at fixed lambda |
| :-- | :-- | :-- |
| Spin 1/2 at each logical site | Exact tensor factor c in (7) | Code only before removing leakage |
| XY coefficient normalized to -1 | Divide by s>0 | Same at degree two |
| Real even exponentially decaying density potential | v(+1)=v(-1)=1/2, all others zero | Same only in reference |
| Instantaneous two-body density interaction | Exact equation (3) | Not established for matched action |
| Small dimensionless interaction | g is chosen within BM's radius | Reference g is small; arbitrary-degree corrections remain |
| Interior simple Fermi points | Fixed pF=k*, sin k*>0 | Reference only; no microscopic theorem invoked |
| Half-filling quartic resonance absent | Positive margin after (4), fixed by (12) | Reference detuned; this alone does not prove RG |
| Sources covered by theorem | Bare c density and logical current, equation (8) | Dressed-source matching remains necessary |
| No extra low-energy cell sectors | Threshold (9) and pressure estimate (10) | C47's finite-order leakage remains |

[proved, boundary qualification] BM chooses the spin seam so that the Jordan-Wigner fermions are periodic. Our exact algebra and thermal comparison hold with that same seam, so this is the convention for the theorem application. For the ordinary local periodic spin ring, Jordan-Wigner gives periodic or antiperiodic fermions according to total parity. We keep that twist in the numerical checks. Changing one seam changes the Hamiltonian norm by O(s), so free energies per cell agree in the thermodynamic limit at fixed temperature, and their ground-energy densities agree. The O(s/M) estimate alone is not a proof of equality of every finite-ring ground-state correlation sequence. The imported correlation conclusion refers to the infinite-chain state obtained with BM's boundary prescription. No boundary-independent finite-size spectrum theorem is claimed.

[proved, scale conversion] Multiplying a Hamiltonian by s rescales time by s. Therefore the physical cell velocity and static susceptibility from (11) are

    v_cell=s vs,
    chi_cell=K/(pi s vs).                             (13)

The density response is (11) with vs replaced by s vs and the prefactor K/(pi s vs). Exponents do not change. For our v, its leading K correction simplifies algebraically to

    K=1-(2 sin k*/pi)g+O(g^2).

The field (12) holds the singularities fixed; we do not identify an interacting density from the free dispersion without a separate particle-counting argument.

[imported, optional response strengthening] Benfatto and Mastropietro, [“Drude weight in non solvable quantum spin chains,” Theorem 1.1, arXiv:1012.0403](https://arxiv.org/html/1012.0403), treats the same instantaneous Hamiltonian. For sufficiently small g its Euclidean current response has the form

    Dhat(p0,p)=(vs K/pi) p0^2/(p0^2+vs^2 p^2)+H(p0,p),
    |H(p0,p)|<=C|(p0,p)|^theta,  0<theta<1.

Here Dhat=-<Delta_x>-G_JJ, J_x=(a†_(x+1)a_x-a†_x a_(x+1))/(2i), and Delta_x=-(a†_x a_(x+1)+h.c.)/2. Thus D=vs K/pi and D/chi=vs^2. After multiplication by s, D_cell=s vs K/pi. Intersecting the two positive smallness neighborhoods verifies this additional application. This is a zero-temperature response result, not a finite-temperature conductivity theorem.

## 5. What follows about central charge

[proved, logical limitation] Equations (11)-(13) establish correlation exponents, density susceptibility and, with the second theorem, Drude weight. They do not state an interacting low-temperature free-energy coefficient, a Casimir coefficient or a complete stress-tensor scaling limit. The word “Luttinger” in their interpretation is not substituted for one of those statements. Normalized density correlations alone cannot exclude an additional neutral critical sector in a general microscopic theory. Here (9) excludes additional gapless syndrome sectors, but the imported response theorem still does not itself construct the complete stress tensor of the interacting logical spin chain.

[proved, sufficient remaining input for g nonzero] For H_exact, it would suffice to establish for its logical Hamiltonian

    f_BM(T)=e_BM-pi T^2/(6vs)+o(T^2),                 (14)

or a complete c=1 scaling-limit theorem. Equation (10) would then transfer the coefficient without further massive matching. The usual conformal interpretation of (14) is an additional identification, not a conclusion obtained just by knowing chi and D. This report does not independently prove (14) at nonzero g or import an additional Bethe-ansatz thermodynamics theorem.

[proved, free member has a complete thermal answer] At g=0, h=1/3, the quadratic logical chain gives its pressure directly. In normalized units,

    f0(T)=h/2-(T/(2pi)) integral_(-pi)^pi
                             log(1+exp[-epsilon(k)/T]) dk.

Subtract the zero-temperature integral of min(epsilon,0). The remaining integrand is

    -T log(1+exp[-|epsilon|/T]).

Outside fixed neighborhoods of the two simple roots its integral is exponentially small. Inside, change variable y=epsilon(k)/T. At each root dk=(T/v*)[1+O(T|y|)]dy; exponential tails justify dominated convergence and the error is o(T^2). Since

    integral_(-infinity)^infinity log(1+exp[-|y|]) dy=pi^2/6,

each root contributes -pi T^2/(12v*), and the pair gives

    f0(T)=e0-pi T^2/(6v*)+o(T^2),
    e0=(-sin k*-h k*)/pi+h/2.

Rescaling by s yields v_cell=s v* and the coefficient -pi/(6s v*). Equation (10) proves the same coefficient for the full frozen H_exact. Thus the thermal central-charge estimator -6v_cell(f-e0)/(pi T^2) tends to 1 for this explicitly specified member. This is an analytic statement, not a fit.

[proved, free continuum interpretation] The low-energy Fourier modes near +k* and -k* are two chiral components of one complex fermion. Linearizing the quadratic generator gives velocities +v_cell and -v_cell; the higher derivatives of the smooth cosine dispersion vanish under the low-energy scaling. There are no other roots and (9) excludes massive-cell contributions. The continuum CAR and filled-sea contractions therefore give the free Dirac theory. To check its central normalization, in chiral units let j=:psi†psi: with contraction <psi(z)psi†(w)>=1/(z-w). Wick contraction gives j(z)j(w)~1/(z-w)^2. For T=(1/2):jj:, the two double contractions in T(z)T(w) give (1/4)*2/(z-w)^4=(1/2)/(z-w)^4. Comparing with c/[2(z-w)^4] gives c_R=1 and similarly c_L=1. In the nonchiral-chain convention the thermal c is 1, not c_R+c_L=2. This identifies the free frozen member with c=1. It says nothing about unrestricted type-sector selection or the wrapper's second mode.

## 6. Comparison with the collaborator

[read, collaborator result] The collaborator's demonstrator note keeps nearest-neighbor physical interactions and obtains

    |Delta_XXZ|=3 |Jp/Jz| kappa^2/[4(4D^2-kappa^2)].

Fixing D,kappa,Jz and taking a nonzero small Jp puts the second-order chain close to XX. Its hopping is proportional to Jp. It correctly keeps the retarded-action and central-charge steps open. Its result is consistent with C47 and with (1)-(6).

[proved, comparison] The three routes have distinct merits and costs.

| Route | Alphabet / maximum physical support | Exact statement obtained |
| :-- | :-- | :-- |
| Collaborator's Jp tuning | 10 / two sites | Arbitrarily weak XXZ reference; matched-action RG still needed |
| C53 diagonal counterterm and field | 10 / four sites | Freely chosen reference g at fixed hopping, quartic detuning; same matching gap |
| C53 direct logical bond replacing PP | 10 / four sites | Exact frozen-sector BM chain, uniform massive threshold and thermal transfer |

The last line achieves a published-theorem application by changing the architecture. It inserts the logical hopping explicitly instead of generating it through virtual port excursions. It is a control model for the infrared claims, not a completion of the collaborator's nearest-neighbor proof. Blocking the four-site model would increase the on-site dimension and change the one-site formulation; we do not use blocking to claim an in-scope nearest-neighbor counterexample.

[proved, wrapper boundary] None of these arguments handles mobile cells, the unbounded length space or the additional length fluid. No c=2 claim for n=1 and no unbounded-c claim follows. The collaborator's selector audit also remains relevant: comparing lower bounds against lower bounds does not establish a positive dressed energy per dilute defect.

## 7. Numerical verification and reproducibility

[numerical] Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py`. Only NumPy, SciPy and the Python standard library are used. The output is [checks.json](checks.json). The checks rebuild the half-cell and full-cell matrices independently of C47's Python imports.

The second-order matrix error is 1.55e-15. Exact Fraction arithmetic verifies q(g) and the target anisotropy for g=0,+1/50,-1/50. The small-lambda diagonalizations verify a cubic residual for H_tune rather than silently dropping it.

For H_exact an open two-cell matrix of dimension 1024 at s=.005625, g=.02, h=1/3 satisfies H G=G H_log to 2.52e-16. The lowest four eigenvalues agree with the two-spin logical matrix to 1.32e-14. The next eigenvalue lies .401148117 above the ground, consistent with the analytic .4 threshold. Partition-function ratios at beta=5,10,20 satisfy (10).

Logical ring ED at M=8,10,12 uses every particle-number sector, with maximum block dimension 924. At h=1/3 the ground sectors are N=5,6,7 for both g=0 and g=.02, rather than N=M/2. The free sector energies agree with the parity-correct free-fermion spectra to 7.11e-15. The free infinite-volume density is .6081734479693927; the short rings display the expected integer rounding. These numbers neither determine the interacting thermodynamic filling nor certify .02 inside BM's unknown smallness radius.

The free thermal integral gives estimators 1.00294432, 1.00071790, 1.00017845 and 1.00004455 at T=.04,.02,.01,.005 in normalized logical units. These support the separate analytic calculation, not an interacting thermal claim.

[numerical, corrected test issue] The first ED run used the uniform vector as an ARPACK starting vector. In an exactly invariant one-particle subspace this caused `ARPACK error -9: Starting vector is zero`. A deterministic seeded generic vector avoids the degenerate Krylov start; the rerun passes and checks eigenvector residuals. No failed run is counted as a successful spectral check.

[proved, final boundary] H_tune only has a theorem-covered second-order reference. H_exact has an exact theorem-covered logical chain in a frozen sector at the same d=10, with four-site interactions. Its interacting response and exponent conclusions are imported rigorous results under the stated smallness and field schedule. Its free member has a proved thermal coefficient and the free-Dirac c=1 interpretation. Interacting thermal c, unrestricted type selection, nearest-neighbor realization and the wrapper remain separate unresolved tasks.
