# C32. A mobile positive folded bus, with finite crossing bounds

[Proved result; incomplete overall task] A modified folded bus can be made
mobile at fixed onsite dimension d0=2736. It has a 28-letter positive core,
exact U(1) and parity, an exactly degenerate 2^n unperturbed code, and two
adjacent cell lengths carrying the same payload. Its second-order bonds
can be made arbitrarily close, including their scalar parts. A finite
parameter schedule proves that all n quadratic bands cross a common field.

[Proved scope] The crossing theorem concerns the common reference payload
bond. The finite mobile wrapper introduces the separately bounded length
dependence; it does not have an exact free Bloch Hamiltonian at nonzero
wrapper ratios.

[Open] I did not prove distinct centers and Sidon-compatible margins for
every n. Consequently the complete stable-phase hierarchy is conditional
on that missing lemma. The two- and three-cell ED below tests the controlled
fast-wrapper model, with active lengths and all legal core excitations.
It is not full microscopic ED of the d0=2736 chain. These are the two
parts of the requested task that remain unfinished.

[Proved record] Detailed definitions and derivations are in
[PROOFS.md](PROOFS.md). All work and generated files stayed in this folder.
No network, external packages beyond numpy/scipy, cluster computation,
subagents, commit, or push was used.

## 1. Mobility and the alphabet

[Proved construction] Direct padding of the supplied 23-letter cell is
not justified by C8's parent theorem. The original vacuum is a sine-wave
ground state, while that theorem requires a positive initialized core.
Also, C25 shows that a helix equalizer does not equalize a second-order
resolvent exactly at finite wrapper scales.

[Proved construction] Add five left-marker letters A_empty,A_v,A_0,A_+,A_-
to the original alphabet. Replace each hopping matrix element by its
positive two-state edge matrix, retaining the original off-diagonal
amplitude. Use positive onsite excited masses and the same number-
conserving conversions. Both markers and the express rail have an extra
energy D. The core now has n data sites between two markers, k=n+2 sites
in total. Its vacuum ground state is uniform on the k token positions.

[Proved] Positivity forces every zero vector into the vacuum sector,
and the connected vacuum path forces equal amplitudes. There is one code
vector for each n-bit data string. An explicit gap bound is

    gamma_core >= min(K,m-|g|,2[1-cos(pi/(n+2))]) > 0.

[Proved] The charge of a data letter is its bit plus its token charge;
marker charges are the token charges 0,0,+1,-1. Every local move preserves
total Q. The code realizes logical fermion number and parity exp(i pi Q).
The port operators change vacuum to excited species, so their code
compression is zero.

[Assumed inherited finite theorem] C8's complete initialized wrapper,
occupation filter, selector reset, and helix construction apply to this
positive core. Their hypotheses are checked in PROOFS Sections 1 and 2;
their unrestricted-ring exclusion proofs are inherited from the named
C8 files rather than re-enumerated here.

[Proved count] The full alphabet is

    selector(4) tensor
      [(endpoint/bulk tag(3) tensor filter(2) tensor core(28))
       direct-sum typed holes(3)]
      tensor endpoint helix tracks(2*2).

[Proved] Its dimension is 4(3*2*28+3)*4=2736, independent of n. The selected
physical cell lengths are n+3 and n+4, both containing n logical modes.
This count is not optimized. The added left marker also costs two extra
powers of D in the logical exchange, compared with the original design.

[Proved] At hole numbers p=1,2 the code factors into the positive charge
ground vector, the same payload code, and neutral filter/helix vectors.
A local marker move transfers a hole. After projection it is a positive
scalar times the identity on both payloads. Thus the leading length
Hamiltonian is the free hard-core hopping chain with

    J_l=u A_b^2 beta_v cos(1/(2b))^(2b),  b=n+3,
    beta_v=[sum_j chi_1({j})chi_2({j,k+1})]^2>0.

[Proved qualification] The physical ring keeps C8/C13's cyclic momentum
gluing. At finite wrapper ratios, the second-order scalar bond generates
length interactions, so the full effective length Hamiltonian is not
exactly free. Free length motion is the leading term of the controlled
hierarchy.

## 2. Equalizing the complete second-order bond

[Proved] Let r_p be the physical endpoint occupation weight. A neutral
helix compensator with expectation 1/r_p makes the projected source map
equal at both lengths. The projected core operator is nevertheless

    B_p=sum_j d_(p,j) B_j+a_pay+b_pay,
    d_(p,j)=<chi_p|D_j|chi_p>.

[Proved] The residual adjacency dependence is real. For small nonzero
w/U it obeys

    1-d_(p,j) <= [w(k+p-1)/(U-w(k+p-1))]^2.

[Proved] PROOFS Section 3 extends C25's bound to noncommuting logical
source matrices. With delta the actual fast charge/helix gap, s the core
scale, and t=s/delta, it proves

    s K2_(p,q)=K2_bare+E_(p,q),
    max_(p,q)||E_(p,q)|| <= C_n(w/U)^2+E_fast(t),
    E_fast(t)=O(t),

with an explicit operator-norm expression for E_fast, retaining the
compensator norms. First choose w/U, then s/delta. The conclusion is
arbitrarily accurate equality at finite nonzero parameters, not exact
equality at finite w and s.

[Proved] The same bound includes the scalar part. The scalar combination
c_22-c_21-c_12+c_11 is a genuine length-density interaction and is included
in the scale budget. It is not discarded as a payload energy shift.

[Proved] The slow charge tunneling gap is of order
U(w/U)^(2n+4), with a positive leading coefficient. It is often far smaller
than U. The equalization hierarchy uses this small gap.

[Proved] Two port insertions on different cuts leave an excited outer
cell, so there is no second-order three-cell term. One port insertion
cannot be undone by marker moves, so the mixed u lambda term vanishes.
Pure marker corrections act as identity on the payload. On open cell
chains all odd powers of lambda vanish by alternate-cell excitation parity.

## 3. Finite parameters for all-n crossings

[Proved] Set the neutral express phase to pi/(2n), initially with zero
loop coupling. At order g^2, the hopping coefficient has the expansion

    T2(D)=-c_n D^(-p_n) Gamma_n+O(D^(-p_n-1)),
    p_n=2n+4,  c_n=2/(n+2)^2,
    Gamma_ij=integral_0^infinity f_i(t)f_j(t) dt.

[Proved] The f_j are return kernels with one conversion at stub site j.
Their first nonzero terms are t^(2j+1)/(2j+1)!. They are linearly
independent, so Gamma_n is positive definite. A 2n^2-dimensional Lyapunov
equation defines Gamma_n by finite rational linear algebra at rational
mass m. No large-n numerical inference is used to prove positivity.

[Proved] Define

    rho_n=det(Gamma_n)/(tr Gamma_n)^(n-1)>0,
    L=2m+12,
    r=min(1,m)/(16L^2),
    B0=288/[(n+2)^2 m^3].

[Proved] A complex resolvent disk and Cauchy coefficient estimate give

    ||D^p_n T2(D)+c_n Gamma_n|| <= 2n B0/[D r^(p_n+1)].

[Proved] It is therefore sufficient to take

    D_n >= max(2/r,8n B0/[c_n rho_n r^(p_n+1)]).

[Proved] Set a_n=c_n rho_n D_n^(-p_n). The further choices

    0<g<=m/4,   g^2<=a_n m^5/(10^5 n),
    eta_n=g^2 a_n,   epsilon=eta_n/4,
    0<t_d<=eta_n m^2/(10^5 n)

keep both the higher-conversion terms and the nonzero loop perturbation
smaller than the crossing margin. At two opposite momenta all n
eigenvalues have opposite signs. Every ordered band crosses on both
arcs. This establishes the finite-parameter crossing part of task 2.

[Proved limitation] It does not establish simple Fermi zeros, distinct
centers, or Sidon-compatible pair sums. A definite diagonal hopping
matrix with a common field is an explicit counterexample to that inference:
all bands cross but every opposite-branch pair has the same momentum sum.

[Open] For the loop family, an all-n proof must exclude every relevant
pair-sum identity on its actual parameter space. C20's independent-
momentum volume argument does not prove that this small parameter family
meets its nonresonant region. Breaking an apparent symmetry and observing
nonzero splittings at small n are insufficient. I did not close this step.

## 4. Scales as functions of n

[Proved] The unconditional hierarchy has the following positive scales.

| Quantity | Bound or definition |
| :-- | :-- |
| Occupied core length | k=n+2 |
| Physical cell lengths | b=n+3, b+1=n+4 |
| Core gap | min(K,m-|g|,2[1-cos(pi/(n+2))]) |
| Selector margin | order t_sel b^(-4), with C8's exact finite formula |
| Helix gap | at least J_h/[(b+1)(2b+1)] |
| Slow charge gap | C_(n+2,2) U(w/U)^(2n+4) to leading order, C>0 |
| Express attenuation in T | D_n^(-2n-4) |
| Gram lower bound | rho_n=det Gamma_n/(tr Gamma_n)^(n-1) |
| Dimensionless crossing scale | eta_n=g^2 c_n rho_n D_n^(-2n-4) |
| Physical narrow scale | at least (lambda^2/s) eta_n/2 |
| Length scale | J_l=u beta_n |

[Proved] For a desired tolerance tau_n, first make the missing-bond error
smaller than tau_n/4 using w/U. Next make the complete fast-wrapper error
smaller than tau_n/4 using s/delta. Compute the positive exclusion gap
Delta_n from the explicit wrapper/filter/selector bounds. Finally set
u=(lambda^2/s)j_l/beta_n and choose a fixed nonzero lambda so that finite
local perturbation errors are below (lambda^2/s)tau_n/4. Every choice is
made after the constants on which it depends. There is no dependence on
the eventual ring length.

[Proved] This procedure closes for every fixed n with tau_n=eta_n/100.
It preserves the crossing signs and leading mobile decomposition.

[Assumed for the stronger hierarchy] To preserve a stable n-mode payload,
one instead needs a simple-root velocity v_*>0 and a payload momentum
margin delta_*>0, then takes

    tau_n=min(eta_n/100,v_* delta_*/100,v_*/100).

[Proved conditional statement] Once those finite witnesses exist, the
same ordered choices close. C13's finite exclusion argument also chooses
an interior length density with a positive mixed margin. No additional
scale obstruction appears at fixed n.

[Open/assumed] The missing all-n payload margin prevents an unconditional
claim that the full stable-phase hierarchy closes for every n. Local
massive-sector elimination and the stated bosonization stability input
are also needed for a thermodynamic c=n+1 conclusion. Finite-cluster ED
and a global finite-volume code-band expansion do not prove that conclusion.

[Numerical, Lyapunov dimensions 2,8,18,32,50] The Gram spectrum at m=1
illustrates why a fixed exponential extrapolation is unsafe.

| n | smallest/largest Gram eigenvalue | determinant lower bound rho_n |
| --: | --: | --: |
| 1 | 1 | 9.2593e-3 |
| 2 | 2.2403e-2 | 1.0602e-4 |
| 3 | 1.9234e-4 | 7.4561e-9 |
| 4 | 1.4413e-6 | 3.3999e-15 |
| 5 | 1.0195e-8 | 1.0885e-23 |

[Numerical illustration, not optimized parameters] Substituting these
values in the conservative schedule at m=1 gives the following base-ten
logarithms. The conversion amplitude saturates the displayed upper bound.

| n | log10 D_n | log10 g | log10 eta_n |
| --: | --: | --: | --: |
| 1 | 29.57 | -92.55 | -365.21 |
| 2 | 38.80 | -160.31 | -635.93 |
| 3 | 50.13 | -257.98 | -1026.45 |
| 4 | 63.58 | -392.17 | -1563.09 |
| 5 | 79.17 | -569.21 | -2271.14 |

[Proved interpretation] These are sufficient bounds, not estimates of
the best achievable scales. They show finiteness at fixed n. They give
no useful uniform energy scale as n increases.

## 5. Exact source calculations and mobility checks

[Numerical setup] mobile.py builds the positive local core pieces and
the exact ordered-coordinate vacancy matrices. check_mobile.py uses
g=.45, m=.5, D=.2, t_d=.08, neutral express phase .4, and loop phase .9.
This is a diagnostic point with a measurable logical interaction, not
a claimed weakly interacting infrared witness.

[Numerical, core dimensions 30 and 88] At n=1,2 the code residual and
charge commutator are zero in the assembled matrices. Code orthogonality
errors are at most 2.3e-16. First-order port compressions vanish. Measured
core gaps are .5004961 and .4196287.

[Numerical, exact source spaces] Charge dimensions are 4,10 for n=1 and
5,15 for n=2. Excited-core dimensions are 24 and 72. The reachable helix
sector contains its ground state and all one-flip modes. Including it
gives source dimensions 480,1440 and 2160,7560. Every source eigenstate
is retained. Pair inverse contractions are evaluated spectrally, without
allocating the largest 57,153,600-dimensional pair matrix.

[Numerical estimator] The following table uses s/delta_charge=.01 and
reports the maximum over all four length pairs of the traceless
difference from K2_11, divided by ||tl K2_bare||. The limiting column
uses the exact averaged operator (9), not an extrapolation. The last
column is the full, uncentered length mismatch in s K2 units.

| n | w/U | slow charge gap | finite mismatch | averaged-limit mismatch | full mismatch |
| --: | --: | --: | --: | --: | --: |
| 1 | .16 | 8.6604e-6 | .0059382 | .0059066 | .0069370 |
| 1 | .10 | 5.0628e-7 | .0021268 | .0021146 | .0076757 |
| 2 | .16 | 1.1229e-7 | .0055363 | .0055028 | .0040924 |
| 2 | .10 | 2.5441e-9 | .0021248 | .0021092 | .0043729 |

[Numerical] The averaged mismatch decreases with w as predicted by the
quadratic missing-bond bound. The full mismatch also depends on s/delta
and need not decrease when w alone changes. At n=2,w=.16 the scaled scalar
length-density interaction is -9.2072e-4. It is retained in the output.
These parameter choices are not small enough to meet the full infrared
error budget.

[Numerical, independent audits] Direct enumeration of occupation patterns
before and after moving the cut agrees with beta_v to 5.6e-17. Direct
products of all physical helix overlaps agree with the transport factor
to 4.4e-15. A raw 240-dimensional one-cell source inversion agrees with
the stabilized spectral result to 1.1e-10 absolute error. Charge changes
of all six endpoint source operators are checked separately and have
zero residual. The core's positive pieces have smallest numerical
eigenvalue no lower than -8.5e-16.

[Numerical limitation] An earlier small-signal calculation, with a
traceless bond of order 1e-10, lost relative accuracy after scalar
subtraction. It is retained in small_signal_diagnostic.json and is not
used as precision evidence. Even the improved n=2,w=.10 traceless
fast-limit comparison has a visible numerical floor. The full operator
bounds in PROOFS, rather than those subtractions, establish convergence.

## 6. Mobile two- and three-cell ED

[Definition, exact scope] The diagonalized model is the fast-wrapper
limit at fixed w=.16. A length configuration p_j in {1,2} has one long
cell and M-1 short cells. Its Hamiltonian is

    sum_j B_(p_j) + lambda sum_j V_(j,j+1)
       -.2 lambda^2 sum_(allowed adjacent length swaps) Swap_length.

[Proved] B_p is the averaged core operator (9), and compensated projected
ports are the bare ports. The length swap acts as the identity on the
entire ordered core state. These operators follow from projection onto
the fast charge/helix/selector ground sectors; they are not obtained by
inserting the already-computed second-order logical bond.

[Numerical setup] ED retains every legal core state in the exact invariant
even-total-excitation sector, at total U(1) charges Q=0 and Q=1. There is
no cutoff on carrier excitations within these sectors. All nM logical
one-particle states are included. Other total charges, and therefore the
full many-payload interaction matrix, are not checked by this ED.

[Numerical estimator] The low-band eigenvectors are mapped to the code
by polar decomposition. The operator-norm difference is taken against
lambda^2 times the independently computed sum of length-dependent K2
bonds and the projected length hopping. The Q=1 results are

| n | cells | sector dimension | code dimension | error/lambda^4 at .12,.06,.03 | smallest observed gap |
| --: | --: | --: | --: | :-- | --: |
| 1 | 2 | 356 | 4 | .06047, .06023, .06017 | .96858 |
| 1 | 3 | 8019 | 9 | .07365, .07363, .07362 | .96085 |
| 2 | 2 | 2288 | 8 | .04523, .04494, .04486 | .57391 |
| 2 | 3 | 91872 | 18 | .06685, .06660, .06654 | .56972 |

[Numerical] The Q=0 sector dimensions are 274,5265,1256,39072 and give
the same fourth-order convergence. Final eigenpair residuals are below
2.4e-12 across all saved runs. The corrected n=2 three-cell run has maximum
residual below 9e-13. The predicted first-order block and code residuals
vanish. This is evidence for the mobile second-order reduction and the
absence of a three-cell term in the tested sectors.

[Numerical failure and repair] The initial n=2 two-cell ARPACK run at
lambda=.03 failed after 10001 iterations, with 8 of 9 vectors converged.
The initial three-cell run also became slow on clustered levels and was
interrupted. Enlarging the Krylov space resolved both. Each subsequent
diagonalization has a 180-second timeout. The completed three-cell n=2
rerun took about 235 seconds for its three Q=1 diagonalizations together.
No failed result is reported as a successful check.

[Limitation] The full mobile microscopic Hilbert space includes vacancy,
filter, helix, selector and wrong-length excitations. It was not
diagonalized. The finite source calculations in Section 5 include the
relevant second-order fast excitations, but they do not turn this ED
into a full microscopic test. Thus task 4 is only partially satisfied.

## 7. Small-n band checks

[Numerical setup] band_checks.py uses g=.05,m=.2,D=1,t_d=.001, express
phase .4 and loop phase .9, then a common field half the smallest
phase-rotated Hermitian hopping eigenvalue. It extracts the quadratic
part of the exact second-order bond, locates roots by bracketing and
Brent solves, and computes velocities by the Hellmann-Feynman formula.

[Numerical, core dimensions 30 and 88] The results are

| n | hopping scale | centers modulo pi | delta2 | delta4 | minimum velocity / hopping scale |
| --: | --: | :-- | --: | --: | --: |
| 1 | 5.7549e-8 | 0 within numerical error | 2.63623 | 1.01072 | 1.93649 |
| 2 | 2.9856e-9 | -.00105711, +.00105711 | 1.16025 | .00422846 | .202721 |

[Definition] delta2 is the minimum separation of distinct Fermi points
modulo 2pi. delta4 compares every two distinct unordered pairs, including
repeated indices, by their momentum sums modulo 2pi. This is a
conservative number-conserving quartic test. Centers are quoted modulo pi
because pair transfer depends on twice the center.

[Numerical limitation] These are floating-point witnesses, not interval
certificates or interacting central-charge measurements. The non-Gaussian
norm divided by hopping scale is .000707 and .002781 at these diagnostic
points. No assertion that those displayed points already meet every
infrared smallness condition is made.

[Numerical, direct leading coefficient] check_asymptotic.py gives relative
Gram errors at D=10,100,1000 of .7065,.1154,.01196 for n=1 and
.8560,.1754,.01829 for n=2. The direct g^2 extraction confirms the
normalization and convergence mechanism used in the finite bound.

## 8. Files and reproduction

[Proved record] The implementation and checks are

| File | Content |
| :-- | :-- |
| mobile.py | Positive core, vacancy matrices, full source measures, helix and bond contractions |
| check_mobile.py | Four length pairs and finite-wrapper equalization scans |
| audit_checks.py | Raw inverse, charge, positivity and direct cut-overlap checks |
| ed_mobile.py | Active-length ED in the fast-wrapper limit |
| band_checks.py | Lyapunov Gram matrices and small-n band witnesses |
| check_asymptotic.py | Direct g^2 resolvent coefficient extraction |
| PROOFS.md | Local construction, all bounds, selection rules and open lemma |

[Numerical reproduction] Run each command inside this folder with
OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.

```text
python3 -B check_mobile.py
python3 -B audit_checks.py
python3 -B band_checks.py
python3 -B check_asymptotic.py
python3 -B ed_mobile.py 1 2
python3 -B ed_mobile.py 1 3
python3 -B ed_mobile.py 2 2
python3 -B ed_mobile.py 2 3
```

[Proved record] ED resumes complete saved sectors and lambda values.
To force a fresh ED, move its corresponding ed_n*_M*.json file aside.
The other checks overwrite only their own outputs in this folder.

[Open next step] The decisive remaining finite-dimensional question is
the all-n nonresonance lemma in PROOFS Section 7. The mobility and finite
crossing estimates do not resolve it. A full microscopic mobile ED would
be a separate numerical extension, not a substitute for that lemma.
