# C39. Two folded buses and a geometric readout

[proved, finite result] Enlarging the initialized core from 28 to 46 letters gives an all-n construction with exact U(1), a complete frustration-free 2^n-dimensional code, and a second-order quadratic reference with 2n simple Fermi points, distinct sea centers and strictly separated unordered pair sums. The construction uses two independent folded buses with different stub hoppings and different port phases. A biased vacuum walker supplies a second geometric profile for the density readout. No independent physical letter for each flavor is required.

[proved, constructive scope] The parameters are specified by a terminating rational/algebraic search for each n. The proof of termination is analytic. It does not extrapolate the original loop numerics, assume genericity on a small parameter family, or infer momentum separation from hopping eigenvalue separation. The complete derivation and local rules are in [PROOFS.md](PROOFS.md).

[assumed inherited wrapper; proved alphabet count] Applying the C8/C32 initialized-core wrapper gives d0=96*46+48=4464. The frozen version fixes the marker positions and selects one length. The mobile version selects lengths n+3 and n+4, both carrying the same n-mode payload. This uses the inherited wrapper and exclusion theorem; no independent full-space enumeration of that wrapper was performed here. The number 46 counts the initialized core, not its spontaneous selector.

[assumed/heuristic infrared conclusion] Local massive-sector elimination and Gaussian stability remain the thermodynamic inputs. Given them, sufficiently weak fixed nonzero couplings give the usual frozen c=n and mobile c=n+1 conclusions. No central charge, every-ring uniqueness, or thermodynamic elimination theorem was established by these computations.

## 1. The mechanism

[proved] Let Gamma(w) be the positive Gram matrix of the folded-bus return filters,

    f_j(t;w)=int_0^t [exp(-s A(w))]_(0,j)
                      [exp(-(t-s) A(w))]_(j,0) ds,
    A(w)=4I-w adjacency(P_n),
    Gamma_ij(w)=int_0^infinity f_i(t;w)f_j(t;w) dt.

The two excited buses have stub hoppings w and wv. Their port phases pi/2 and pi give leading logical hopping

    T=-J[Gamma(w)+i eta Gamma(wv)],
    J=2 lambda^2 g^2 [q^(2n+2)/S^2] D^(-2n-4),
    S=sum_(ell=0)^(n+1) q^(2ell).

The factor i is a relative physical phase between two virtual channels. It cannot be removed by conjugating the logical data basis.

[proved] Put D_w=diag(w^(2j)) and C(w)=D_w^(-1)Gamma(w)D_w^(-1). This matrix has a positive limit

    C(0)_ij=(2i+2j+2)! /
      [(2i+1)!(2j+1)!(8)^(2i+2j+3)].

The generalized eigenvalues of D_v C(wv) D_v against C(w) satisfy

    r_j(v)=kappa_j v^(4j)(1+o(1)),  kappa_j>0,
    kappa_j=c_(j+1)/c_j,
    c_k=det C(0)[0:k] det C(w)[k:n]/det C(w).

The determinant expansion has a unique lowest-power minor for each coefficient. Its Newton polygon gives these distinct powers. The corresponding eigenvectors have explicitly computable limits supported on successive tails, given by columns of the inverse tail matrices of C(w).

[proved] Replace the vacuum kinetic edges by

    (q|j>-|j+1>)(q<j|-<j+1|).

Their unique zero vector has amplitude q^j. A repeated onsite operator which reads the data bit only when the vacuum token sits there compresses to q^(2j+2)/S times that data occupation. Choose q=(wu)^2. A retained weak readout then produces the field J mu D_w E_u D_w, E_u=diag(u^(4j)). This changes the code that underlay C28's four-profile obstruction. It does not contradict that obstruction.

[proved] A positive congruence removes D_w from the zero-energy equation, leaving

    h(k)=-2A cos k+2 eta B(v) sin k+mu E_u,
    A=C(w), B(v)=D_v C(wv)D_v.

At mu=0 its roots are theta_j +/- pi/2, with theta_j=-arctan(eta r_j). Thus v makes the centers lacunary. The derivatives under the field are opposite at the two roots,

    dp_(j,+)/dmu=-s_j, dp_(j,-)/dmu=+s_j,
    s_j=(x_j^T E_u x_j)/
        [2 sqrt(1+eta^2 r_j^2) x_j^T A x_j].

Choose u small first, then v small. The limiting eigenvectors show s_j has a positive leading coefficient times u^(4j), so the slopes are also lacunary.

[proved, all pair comparisons] Take consecutive center and slope ratios below 1/32. Pair-sum comparisons with different flavor multisets have a nonzero lacunary center combination. Those using the same flavor multiset split linearly in mu by a nonzero lacunary slope combination. Odd reciprocal half-periods remain separated by a fixed amount. This includes pair transfer, unequal-width backscattering, umklapp and repeated-point pairs. For a sufficiently small positive mu,

    delta4 >= min(|theta_(n-1)|/4, mu s_(n-1)/4)>0,
    delta2 >= delta4.

Every zero is simple, and the physical minimum velocity is positive. Distinct centers survive the small field.

[proved, finite parameters] Fix w=1/4. The determinant ratios and limiting eigenvectors specify terminating searches for eta,u,v,mu. Choose finite D next to approximate the Gram limit, then nonzero g to control all g^4 terms, and finally lambda to control higher port orders and readout leakage. All choices depend on n, not on the eventual ring length. The complete finite-g second-order quadratic projection retains the strict band conditions. Its non-Gaussian part can be made arbitrarily small relative to the retained scale. This is a finite-matrix statement, not the thermodynamic elimination assumption.

## 2. Exact code and microscopic definition

[proved] The alphabet has ten shared vacuum/data/marker letters and eighteen excited letters for each bus. The complete legal-cell dimension is (13n+14)2^n. Positive local kinetic terms, onsite excited mass and the grammar penalties give

    gap_core >= min(K,1-|g|,1+q^2-2q cos(pi/(n+2))).

The code is the weighted vacuum path tensor the unrestricted n data qubits. Its dimension is 2^n. Charge is data occupation plus the excited token charge 0,+1,-1. Every local term and port exchange conserves that charge.

[proved] The density readout is a retained order-lambda^2 perturbation. It is not included in the frustration-free zero-kernel statement. It need not preserve the code exactly; its leakage is included in the higher-order budget.

[proved] For each bus let H_r be the complete excited-cell matrix, and B_r the three pair-source maps with their port phases and right neutral parity string. The exact coefficient of lambda^2 is

    K2= -q^(2n+2)/S^2 sum_r b_r^2
       B_r^dagger (H_r tensor I+I tensor H_r)^(-1) B_r,
    b_1=1, b_2=sqrt(eta).

There are no cross terms between buses and no three-cell term at this order. The full three-cell coefficient is K2_12+K2_23. At finite g this bond is generally interacting; replacing it by the displayed Gaussian limit would omit real terms.

[numerical, independent local audit] [literal.py](literal.py) constructs the physical 2116 by 2116 bond and its boundary operators. Its restrictions at n=1,2,3, dimensions 54,160,424, agree with the independently built gauge matrices to maximum row-sum error below 5.9e-16. The charge commutator vanishes entrywise. Positivity and the complete kernel are proved in PROOFS, independently of these floating-point checks.

## 3. Certified small-n reference points

[proved, rational certificates] [certify.py](certify.py) uses exact Fraction arithmetic, polynomial interpolation, Sturm root isolation and bounded arctangent series. It certifies every root and every pair-sum comparison, including repeated indices. The points use w=1/4, v=1/1024, u=1/64, eta=1/100 and mu=10^(-(8n+10)). Conservative rounded-down bounds are

| n | Minimum root separation | Minimum pair-sum separation | Centers, decimal display |
| :-- | :-- | :-- | :-- |
| 1 | >3.14 | >5.11e-16 | -0.00999966669 |
| 2 | >0.0166 | >1.23e-28 | -0.0166086805, -5.39749e-15 |
| 3 | >1.42e-14 | >2.36e-41 | -0.0217407823, -1.42597e-14, -2.32720e-27 |

[proved, scope] These certify the limiting reference, not the finite-D diagnostic below. Their strict margins can be transferred to finite D and g by the proved approximation argument. The fixed numerical choices in this table are not asserted to work for arbitrary n. The all-n definition is the terminating search. The certificate runs took 2.5,6.5,8.2 seconds.

## 4. Complete effective bonds and two-/three-cell ED

[numerical, full source ED] [model.py](model.py) and [checks.py](checks.py) compute the full-in-g K2 using every excited-cell eigenstate. The per-bus source dimensions are 24,72,192 at n=1,2,3. The two-cell code matrices have dimensions 4,16,64, and the saved three-cell coefficients have dimensions 8,64,512. They are in [bond_n1.npz](bond_n1.npz), [bond_n2.npz](bond_n2.npz) and [bond_n3.npz](bond_n3.npz). No charge truncation enters these effective bonds.

[numerical, diagnostic parameters] The finite-coupling checks use a=3.5, D=.2, w=.8, v=.65, g=.4, q=.8 and eta=.6, plus .03 lambda^2 times the local conditional density readout. These values resolve the second-order matrix numerically. They are not claimed to satisfy the conservative all-n momentum hierarchy. Their scalar K2 entries are -0.06584624723, -0.03313431606 and -0.01843338048. Hopping entries and the complete matrices are saved, including the non-Gaussian remainder.

[numerical, actual cluster ED] [ed.py](ed.py) retains every legal state in exact Q=0,1 sectors with even parity for each excited bus. It diagonalizes the physical cluster before any code projection. Polar overlap with the unperturbed code defines the canonical comparison. The following Q=1 results show the fourth-order remainder after subtracting the complete second-order bond and readout.

| n | Cells | Sector dimension | Code dimension | Error/lambda^4 at .1 | At .05 |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | 2 | 338 | 2 | .030461 | .030374 |
| 1 | 3 | 5265 | 3 | .044852 | .044782 |
| 2 | 2 | 2224 | 4 | .020624 | .020539 |
| 2 | 3 | 60864 | 6 | .033706 | .033604 |
| 3 | 2 | 9110 | 6 | .013993 | .013922 |
| 3 | 3 | 356325 | 9 | .023737 | .023639 |

[numerical] Q=0 passes separately at all six sizes. First-order compression is zero, and the largest eigenvector residual is below 2.8e-12. An initial n=3 three-cell ARPACK run hit the 180-second cap. A block eigensolver initialized with the code completed that check in 118 seconds. It did not compute an additional excited level, so that run has no reported next-level gap. Smaller ARPACK runs do include that level.

[numerical, scope] The finite-coupling cluster ED covers Q=0,1 only. All charges are included in the full source-bond calculation. No full mobile-wrapper ED, wrong-marker-sector ED or central-charge estimator was run. The full physical 4464-letter chain was not diagonalized.

[numerical, normalization audit] Direct g^2 pair inverses at D=10,100,1000 approach the claimed Gram coefficient. Their relative errors at D=1000 are .0223,.0302,.0380 for n=1,2,3. Exact rational scalar inverses agree with independent Lyapunov Gram calculations to below 3.6e-18. These are consistency checks of the normalization and asymptotic power.

## 5. Mobile transfer and what remains assumed

[proved conditional on the inherited wrapper] C32's noncommuting source-resolvent estimate applies to the two-bus source. It equalizes the complete second-order bond between the four length pairs, including its scalar part, to arbitrary prescribed accuracy at finite nonzero wrapper ratios. The biased core stays frustration free and has the same payload at both selected lengths. The readout travels with the occupied letters and has the same compressed sum in both lengths. A finite exclusion of length densities gives the mixed C20 momentum margin once the payload points are fixed.

[assumed] This application inherits C8's initialized wrapper, occupation filter, global exclusion and identity-transport gluing. It does not constitute a new independent audit of all those finite constructions. Their alphabet count is explicit in PROOFS Section 7.

[assumed/heuristic] The remaining physical input is the irreducible low-frequency massive-sector elimination and weak-coupling RG matching described in C40 and ir_control. The narrow velocities and momentum margins enter that budget. Conversion g must suppress the normalized g^4 interaction before lambda suppresses higher port orders. Finite-cluster perturbation theory is not substituted for this thermodynamic statement.

[proved, relation to adjacent work] C33/C37 certify original-alphabet targets through n=23. C38 shows why a bandwidth hierarchy alone does not control the Fermi zeros. This construction instead changes the normalized matrix pencil, and uses a separate geometric field to remove its residual half-filling coincidences. It does not duplicate those original-loop certificates. C40's thermodynamic remainder remains an assumption here.

## 6. Reproduction and handoff

[proved, file record] All research files are in this folder. No network, subagents, cluster computation or edits outside the authorized folder were used. The project STATE.md was left untouched under the write restriction; [STATE.md](STATE.md) is the local handoff.

[numerical, commands] Use Python 3 with numpy and scipy. The exact certificate uses only the standard library.

    OPENBLAS_NUM_THREADS=1 python3 -B literal.py
    OPENBLAS_NUM_THREADS=1 python3 -B checks.py
    python3 -B certify.py 1
    python3 -B certify.py 2
    python3 -B certify.py 3
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B ed.py 3 3 1

[numerical, remaining ED commands] Replace the last three arguments by n in {1,2,3}, cells in {2,3}, and Q in {0,1}. Each invocation has a 180-second cap. The main proofs have not been independently refereed.
