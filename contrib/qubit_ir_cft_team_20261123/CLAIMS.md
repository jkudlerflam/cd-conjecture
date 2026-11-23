# Claim ledger

Status vocabulary: **proved**, **numerical**, **conditional**, **refuted**,
**unresolved**. Conditional results may be exact algebra within a specified
effective theory; their hypotheses are not proved for the target.
“Audited” means independent AI checking, not human review. Detailed reasons
are in [audit/AUDIT.md](audit/AUDIT.md).

| ID | Claim | Status | Evidence and scope | Audit |
|---|---|---|---|---|
| C01 | Fixed literal d=2 model, a=19/32, b=13/32, h=1783/1000, full fields including OBC ends. | proved | Definition; `audit/check_ring.py` independently compares frames and degree-one edges. | Accepted |
| C02 | Uniform rotation gives Delta=1, V=4 and imaginary hopping magnitude 3/16. | proved | Pauli identity; finite-chain frame checks <6e-15. | Accepted |
| C03 | This target satisfies the equal-hop tuning condition of delta=1/8. | refuted | Formal t2=1/8 gives direct/t2=3/2. Rotation does not control the expansion. | Accepted rejection |
| C04 | Saved entropy coefficients are near two on accessible scales. | numerical | `evidence/read_tables.py`: iMPS 1.97265 +/- .00986, OBC 1.95401 +/- .07500; rounded regressions. | Independently reproduced |
| C05 | L=512 is a converged target entropy datum. | refuted | Last-step filter fails; L=384 also retains a nonzero correction. | Accepted rejection |
| C06 | iMPS and OBC use independent microscopic implementations. | refuted | Same TeNPy model. Separate small-system literal-spin checks exist. | Accepted rejection |
| C07 | Raw target complex correlators, states and pair momenta are recoverable here. | refuted (access inventory) | Including ignored files: no HDF5 states or target raw arrays. Off-machine availability unknown. | Accepted |
| C08 | Exporter check zeros are measured zero errors. | refuted | `meas_obc.py` writes placeholders; its separate dense contraction test is different. | Accepted rejection |
| C09 | Fixed-target OBC odd-even data alternate in sign or include L=192. | refuted | Target has only L=64,96,128, all positive; other points supply the sign changes/L=192. | Accepted rejection |
| C10 | Microscopic N generates an exact U(1). | refuted | Pairing changes N by two; [H,N]!=0. No separate exact valley charge is established. | Accepted rejection |
| C11 | Staggered chirality equals smooth valley imbalance in the unbroken dictionary. | refuted | Translation-odd versus translation-even; no statement about every possible broken phase. | Independently derived |
| C12 | Staggered chirality and range-two bond can overlap the two pi-pair quadratures. | conditional | Exact symmetry characters match; leading nonzero overlap unmeasured. | Accepted with hypothesis |
| C13 | The microscopic compactification is the proposed integer bosonic vertex lattice. | unresolved | Explicit conditional ansatz; zero modes/current algebra not identified. | Accepted limitation |
| C14 | Generic allowed vertices in that ansatz are cos[2 sqrt(2)(r theta_c-s phi_s)]. | conditional | Parity, momentum and symmetry-orbit projection; balance/generic incommensurability required; special resonances separate. | Independently derived |
| C15 | A genuine generic pi Gaussian primary of exponent eta_pi has an allowed doubled vertex with x=2 eta_pi. | conditional | Charge/momentum doubling and homogeneous Gaussian dimensions; not an unidentified fitted function or descendant. | Independently derived |
| C16 | With the matched two primaries the primitive quartet is the only relevant generic harmonic even with mixing. | conditional | p=.9193, d_theta=.5083; x_(r,s)=2[p(r-t s)^2+s^2/d_theta]. x4=1.8386; s!=0 implies x>=2/d_theta≈3.9347; special resonances excluded. | Independently derived |
| C17 | Known exact symmetries forbid the quartet. | refuted | It passes all listed momentum/symmetry projections; accidental or critical-manifold zero still possible. | Accepted rejection |
| C18 | Quartet target coupling or locking length is determined. | unresolved | No thermodynamic matching at Delta/V=1/4; dilute/finite-level calculations cannot be transferred. | Accepted limitation |
| C19 | Mixing alone reconciles K_density=.675, eta_pair=.9193 with unit total-density normalization. | refuted conditionally | Canonical K_density*eta_pair>=1 versus .62053; relative phase-gradient admixture cannot fix it. Matching/carrier/finite-scale errors remain possible. | Independently derived |
| C20 | The actual pseudospin radius is Ks=1.967 in an unmixed dictionary. | unresolved | Equality needs an unmixed matched dictionary; mixed backscattering lower bound survives. | Accepted limitation |
| C21 | Small-ring pi pair, bond and chirality reach a shared pole. | numerical | L=8,10,12 ED; independent dense L=8 reproduction <1e-14. | Accepted |
| C22 | Opposite antiunitary characters or shared pi probes count distinct gapless modes. | refuted | L=8 imaginary bond/chirality cross residue has unit normalized squared coherence at one pole. | Accepted counterexample |
| C23 | Pair3 static maximum is always at pi on the studied rings. | refuted | Moves to adjacent bins on L=10,12; pi weight remains. No large-L extrapolation. | Accepted numerical check |
| C24 | Local response has a loose causal ceiling v<=153e/16 in physical units. | proved bound | Shifted bond norm51/32; onsite interaction picture; overlapping-bond commutator bound. No speed estimate or positive lower bound. | Independently derived |
| C25 | An unbroken symmetric two-boson Gaussian theory has speed pairs +/-v1,+/-v2, possibly unequal. | conditional | Symmetry-block Hessian; squared speeds are eigenvalues of H_u H_v. | Accepted algebra |
| C26 | The actual IR is a specified two-mode torus/orbifold/rational theory with known velocities and radii/mixing. | unresolved | Dynamic/current-algebra/compactification/longer-scale data missing. | Accepted limitation |
| C27 | The target is asymptotically critical with compact unitary Virasoro c=2. | unresolved | Near-two entropy slopes are finite-scale evidence; crossover to fewer modes or a gap remains possible. | Accepted final limitation |
| C28 | Proposed measurements discriminate the stated alternatives without assuming their identification. | conditional (measurement protocol) | Exact operators, anomalous pair data, matched current algebra, parity seam, explicit resource gates; no production authorized. | Accepted after corrections |
| C29 | The two scalar exponents alone force a relevant lock. | refuted | Abstract same-exponent example p=1.197775, d_theta=.229825, t=0 assigns the pair fit to neutral double transfer; every generic lock has x>=2.39555. No target overlap claimed. | Independently checked |
| C30 | For a matched pure pi annihilation primary, the normal/anomalous ratio scales as g4 r^(2-x4) in first order. | conditional | Gaussian charge-neutral three-point integral; 1/4<x_O<1/2, finite nonzero speeds, cutoff/charge matching and small perturbation required. Raw P can mix O and O-dagger. | Independently derived by auditor |

Dimensions inferred from the printed exponents require their identification
as asymptotic bulk values. Their estimator errors are not certified intervals
and do not include the unknown systematics. Compactification and a stable
microscopic CFT remain open.
