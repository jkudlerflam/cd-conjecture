# C15 local state

Completed 2026-10-27. All requested analysis is written in RESULT.md, with an explicit residual problem rather than a claimed global proof. This is local working memory because the task forbids edits outside this folder.

## Verdict

[Proved / assumed standard RG and exact-solution inputs.] c<=1 holds on the quadratic, uniform U(1), and zero-field XYZ loci; in their regular controlled continuum basins; and in the regular strong-repulsive constrained-liquid regime. Large-field and large-|V| gapped regions are controlled. Singular zero-velocity and degenerate neighborhoods do not inherit a relativistic c-theorem bound. A finite-coupling Ising times LL point is not excluded by a universal relevant-operator count.

## Results to retain

- Convention: H_spin=sum[Jx XX+Jy YY+Jz ZZ+D(XY-YX)-hZ]. Fermions have tau=-(Jx+Jy)+2iD, delta=Jx-Jy, V=4Jz, mu=V-2h. Do not confuse h with C6's opposite field convention.
- Free s^2=4D^2-(Jx-Jy)^2. Chiral c=1 interior: s^2>0, |h|<sqrt((Jx+Jy)^2+s^2). For s^2<0, h=+/- (Jx+Jy) gives Ising, or two Majoranas when both are zero. Thresholds can be quadratic, cubic or flat.
- K>=1/2 is false globally. Strong-repulsive no-adjacent gas has K=(1-rho)^2, approaching 1/4 as rho->1/2. At rho=1/3, K=4/9; pairing and the leading density umklapp are irrelevant.
- Atom-pair conversion cos(2theta1-theta2) has dimension 1/K1+1/(4K2) only for independent fluids. General phase covariance M gives M11-M12+M22/4. Momentum mismatch and locality are necessary tests.
- On delta=0, the pair band is E=V+4|tau|^2 cos^2(P/2+arg(tau))/V-2mu, with |2|tau|cos(...)|<|V|. Its single minimum is shifted with the atomic center. Low-energy pair binding starts at V=-2|tau|, where phase separation also starts; longer clusters preempt isolated pairs.
- At tau=0, V>>|delta|, mu~delta^2/V, second-order virtual pairing gives two-site hopping -1/2 and diagonal -sum_00 1/(1+n_left+n_right), in units delta^2/V. It equals -L+2N+graph_Laplacian+(2/3)N3. All one-sublattice uniform states attain E=-L+2N. This is a ferromagnetic degeneracy, not two LLs.
- Exact line: tau=0, 0<2|delta|<V, mu_f=V-sqrt(V^2-4delta^2), or h_f=sqrt(V^2-4delta^2)/2. Each shifted bond is a positive rank-one projector. An invertible staggered filter proves L+1 ground states on even rings. Exact one-defect branch E_-(k)=V-sqrt(V^2-4delta^2 sin^2(k/2)) is quadratic, k being two-site-cell momentum. The filter does not imply spectral similarity.
- Total N in that second-order theory is accidental. A third-order source for particles at distance three has coefficient 2delta^3/(3V^2). Exact staggered charge survives, total charge does not.

## Numerical evidence

- checks.py, full parity-sector ED L=10,12,14,16. Ising and XX calibrate entropy and independent Casimir estimators. Attractive probe (-.55,-.45,-.55,.1,0) has c_chord 1.12275 -> 1.02716 and c_half fit 1.00162. No c>1 evidence.
- Free collision (-.75,-.25,0,.25,1) has entropy near 1/2 but exact cubic dispersion, verified through q=2pi/256. Nearby untuned interacting probe has no c=3/2 evidence.
- Fixed-N t-V checks L=24,N<=4 support clustering below V=-2. Repulsive constrained gas at (L,N)=(12,4),(15,5),(18,6) gives K_structure .45325,.44929,.44738 and c_energy 1.01829,1.01164,1.00805.
- singular_corner.py, L=24,N<=6 resolved by staggered charge, confirms second-order bound and equality sectors.
- factorization.py, dense L=4,6,8, confirms ground degeneracies 5,7,9 and exact defect branch to 3.2e-15.
- operator_checks.py verifies the entire second-order operator at L=6 (zero discrepancy) and the third-order source at L=8 (2/3).
- corner_probe.py at (.49,-.51,2,.01,sqrt(15)), L=10..16, c_chord falls .86784 -> .26271 and half entropy fails logarithmic growth. Not a tuned critical-point exclusion.

## Remaining work

1. Perturb the exact singular line and determine whether more than one finite velocity can appear. Need higher-order effective theory or larger numerics, not a c-theorem through the z=2 point.
2. Find and tune an intermediate-coupling Ising/LL candidate, specify symmetry/locality of its operators and compute the microscopic coupling-map rank. Four parameters do not exclude the possible continuum codimensions.
3. Verify literature phase-diagram scope against original sources when network use is permitted. No source download, DMRG or global scan was done here.

No computation is running. Sibling C16/C17 progress was read at completion; their results were not independently rerun. No file outside this folder was modified.
