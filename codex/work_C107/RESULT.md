# C107. Compactness and semicontinuity

2026-10-29. Investigation completed. Full derivations are in [PROOF.md](PROOF.md), and the small checks are in [CHECKS.md](CHECKS.md).

[Open] I did not prove a finite universal F(d), prove upper semicontinuity of physical A, or construct normalized fixed-d bonds with A tending to infinity. The compactness route remains possible, but its local boundedness input is not a consequence of the elementary properties of ground energies.

## What was established

[Proved] The literal set M(h)=1 is unbounded because scalar shifts and coboundaries are free. The unit sphere in the quotient by those gauges is compact. Thus the geometric premise can be repaired without changing A/M.

[Proved] For Q_L=(6/pi)L(Le-E_L),

    |Q_L(h)-Q_L(g)| <= (12/pi)L^2 M([h]-[g]),
    |E_L-Le| <= 2M.

Each Q_L is continuous. The continuity estimate grows with L, and the boundary estimate only gives |Q_L|<=12ML/pi. Neither controls the limit.

[Proved] A raw full-length limit need not equal the physical coefficient. For the ferromagnetic XX chain with hopping t and M=t,

    Q_L -> 2t on even lengths,
    Q_L -> t/2 on odd lengths,
    physical A = 2t.

Thus the ordinary limit does not exist and liminf gives A/4. In the gapped classical antiferromagnetic Ising chain, liminf Q_L is minus infinity because odd rings contain a frustrated bond. A vacuum-sector prescription is necessary. These examples do not invalidate the prescribed conformal vacuum expansion.

[Proved] Uniform two-sided convergence fails near an ordinary Ising critical point. Gapped bonds arbitrarily close to it have A=0, while the critical bond has A>0. With the field detuning chosen as L^-3, finite-L continuity forces Q_L to approach the critical value even though each detuned bond has zero asymptotic coefficient. This is compatible with upper semicontinuity, so it is not a counterexample to that property.

[Proved] The raw Q_L is not universally monotone. For the critical XY chain

    H_gamma = -sum [(1+gamma)XX/2+(1-gamma)YY/2+Z],
    M_gamma = sqrt(1+gamma^2),

Euler-Maclaurin gives, at fixed gamma>0 on even rings,

    Q_L = gamma + 7pi^2(4gamma^2-3)/(240gamma L^2)
                    + O(L^-4).

The coefficient changes sign at gamma=sqrt(3)/2. It tends to minus infinity as gamma tends to zero while M stays bounded. Consequently no finite bond-independent C makes Q_L+C M/L^2 eventually decreasing for every bond. This also holds for a d-dependent C, since the entire family has d=2. More general corrected quantities are not excluded.

[Proved logical obstruction] The compact family of abstract energies

    t in [0,1],
    e(t) = -3t^2,
    E_L(t) = L e(t)-t/(1+Lt^2)

has analytic finite-L energies, concavity in t, order-L Lipschitz control, uniformly vanishing boundary errors, length concavity, positive doubling defects, and open-chain almost additivity after adding a boundary constant. Yet

    lim_L Q_L(t) = 6/(pi t) for t>0,
    lim_L Q_L(0) = 0.

The coefficient is finite everywhere and unbounded on a compact parameter set. This is not a local spin-chain construction. It proves that the listed analytic and subadditivity properties do not suffice, even together. Positive degree-one homogeneity can also be imposed by taking the perspective of this energy family.

## The remaining hypothesis

[Proved] Let C be the normalized bonds with a well-defined finite physical A. Define

    A^*(h) = lim_(r down to 0)
             sup {A(g): g in C, M(g-h)<r},

with empty supremum zero. This is a nonnegative upper semicontinuous majorant, possibly infinite. The following are equivalent:

1. A is globally bounded on C.
2. A is locally bounded near every point of the compact quotient sphere.
3. A has a finite upper semicontinuous majorant on that sphere.

An exact extension additionally requires agreement with A on C. For the bound itself a majorant suffices. Assigning an arbitrary finite value at a nonconformal point does not establish local boundedness. Defining its value by A^* makes semicontinuity automatic but leaves finiteness as the original missing theorem.

[Proved conditional implication] If A diverges along normalized bonds, compactness supplies an accumulation point where A^*=infinity. A Lieb-Robinson velocity ceiling would force unbounded total central charge, but not many modes with velocities bounded below. Abstractly, k unit-charge sectors with v=k^-1/2 already have A=sqrt(k) although every velocity vanishes. No fixed-d realization is claimed.

[Proved conditional implication] A quantitative finite cover by neighborhoods with bounds B_j gives F(d)=max_j B_j. Another sufficient input is a common finite length L_0 and a certified correction r_d(L_0) with A<=Q_(L_0)+r_d(L_0). Then

    F(d) <= 12L_0/pi+r_d(L_0).

No such cover or correction was obtained. Compactness alone does not make the bound explicit.

## Checks and scope

[Numerical] The XY fermion sum agrees with nine full spin diagonalizations at L=4,6,8 to maximum error 1.95e-14. Analytic finite sums and quadrature at L=64,128,256 reproduce the sign-changing correction. The XX even and odd sequences and the abstract model's inequalities also pass their stated checks. Run `python3 checks.py`.

[Proved scope] Only work_C107 was modified. No heavy numerics, cluster jobs, network shell commands, or shared-project edits were used. The local STATE.md is the handoff required by the folder confinement rule. Primary literature was consulted through the browsing tool and is linked at the relevant statements in PROOF.md.

[Open next step] The compactness route needs a new lattice estimate establishing local boundedness near every singular bond, or a finite-scale upper observable with a controlled tail. Repeating continuity, open-chain subadditivity, or a C M/L^2 correction cannot supply that estimate.
