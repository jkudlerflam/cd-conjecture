# C102 progress

2026-10-29

[Proved conditional on the standard XXZ Bethe/CFT formula] The fixed-transverse-coupling XXZ line violates convexity. In spin units, A(-1/2)=3 sqrt(3)/8, A(0)=1, and A(1/2)=3 sqrt(3)/4. The midpoint defect is 1-9 sqrt(3)/16>0. All three bonds are critical and have gauge-optimized M=1/2.

[Proved conditional] Completed the full curvature proof. A is strictly concave in XXZ anisotropy and superadditive on its critical cone, with equality only on a common ray. The failure of subadditivity persists on the antiferromagnetic subcone. A increases with both transverse and longitudinal coupling separately inside the cone. Its Heisenberg supporting plane proves A<=pi M for this family.

[Proved] Completed the gauge and extreme-bond audit. The gauge envelope equals A. Two exposed Ising extreme bonds have A=0 while their midpoint has A=2 at M=1. Convexity survives on cones with linear bulk energy and on the free hopping disk, where A=2M. A commensurate complex-hopping segment has the opposite Jensen sign, ruling out global concavity as well.

[Numerical] Completed 136 recorded eigenpairs for XXZ, XYZ, spin-1 bilinear-biquadratic, seeded critical XXZ/Ising, complex hopping, and classical extreme bonds. Largest matrix dimension 12870, largest residual 2.29e-11. Exact XX energies, AKLT energies, full/sector comparisons, gauge/unitary identities, and commensurate hopping identities passed. Spin-1 critical-segment midpoint defects are positive in two estimators, but the thermodynamic sign is unresolved by these small sizes.

[Completed] RESULT.md gives the conclusions, PROOF.md supplies the derivations, NUMERICS.md gives the tables and limitations, and checks.py plus summarize.py reproduce the numerical record. The Bethe component uses exact velocity formulas and bulk quadrature, not a finite-volume root solver. No proof for arbitrary two-level or projector bonds is claimed.

[Scope] Read the project memory, the full background including its round-one update, REPORT sections 9.68 through 9.71, and the C101 results. No network commands or cluster jobs will run. Only this task folder will be modified, so project-level STATE.md will remain untouched. A local STATE.md will record the outcome.
