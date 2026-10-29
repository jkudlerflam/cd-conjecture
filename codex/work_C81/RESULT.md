# C81 result

[Heuristic, conditional conclusion] The explicit one-vacancy lift of C75
has d=16 and predicts c=15+1=16>15. It remains a conditional candidate.
The finite bond and small ED checks are complete. Uniform finite-J
matching, mobile defect stability and the microscopic conformal
identification are not established.

[Definition] The alphabet is vac plus
L_a,R_a,D_(A,a),D_(B,a),X_(A,a),X_(B,a), a=0,1, and M_0,M_A,M_B.
With n the occupied projector and iota the occupied inclusion, the bond is

    h16=-t sum_a(|vac,a><a,vac|+h.c.)
        +J(iota tensor iota)h15(iota tensor iota)^dagger
        -mu(n tensor I+I tensor n)/2.

The complete h15 move table, charge assignments and parameter schedule are
in [PROOF.md](PROOF.md). [bond16.py](bond16.py) is an independent literal
implementation. Every half-site term of h15 stays inside the occupied-pair
block. "One vacancy" means one added local state, with a finite vacancy
density, not one missing particle in the whole ring.

[Proved] Occupied number N, marker number M and the integer internal charge
Q are conserved. On the payload code Q=2N_log. The proposed locking is
N=16M, an energetic condition rather than an identity on the full space.

[Definition; inherited certified starting member] Set n=15, b=16,
D=1000, Delta=4-2cos(pi/17), phi_A=pi/30, phi_B=0, psi_A=pi,
psi_B=pi/2, with port phases -pi/4,+pi/4. Set

    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_16,
    a0=sqrt(2/17)sin(pi/17), t_d=beta 2^-84,
    g=2^-256, K_pen=2^1152, lambda=2^-640,
    f=lambda^2 a0^4 beta^2 g^2 2^-140.

Use C75's b=16 midpoint marker selector. A nested family allowing further
smallness is g=2^(-256-r), K_pen=2^(1152+4r), lambda=2^(-640-u),
J/t=2^-v, with the field adjusted by the same formula. Choose r, then u,
then v, and tune mu to half occupied density. All parameters are fixed
before L tends to infinity. The last choices are conditional on matching
constants, not a certified numerical interval. In particular K must grow
if g is reduced indefinitely. The finite-K normalized error contains
1/(K g^2).

The transfer audit is as follows.

| Part of the d24 analysis | C81 result |
| :-- | :-- |
| Ordered-coordinate squeezing and cyclic charge-sea projection | Proved verbatim for arbitrary h15. The coefficient is alpha_(L,N)=N/L-sin^2(pi N/L)/(N L sin^2(pi/L)). Both half-flux branches and the fermion parity twist are retained. |
| Phase-model vacancy fluid | Same free reference. K_gap=4 at half filling, K_particle=1, and alpha=1/2-2/pi^2. The Gaussian interpretation is a stated continuum input. |
| C75 finite band certificate | Inherited computer-assisted result. Thirty simple roots and all 465 unordered pairs give delta_2>0.0145540, delta_4>1.1444e-7. It was read, not rerun. |
| Fractional marker harmonics | Same derivation with b=16. The marker correlation exponent is 1/128. Every new uniform charge harmonic has reference dimension at least four. The potentially dangerous physical momentum margin exceeds 3.57625e-9. |
| Nonzero-port frozen sector selection | Requires a new proof because bad-pair projectors do not commute. The full-cell proof below supplies it using C75's all-length lower bound. |
| Frozen defect pressure | Transfers after replacing legal/illegal sectors by full finite-K cells. A new thermal bound is given in PROOF.md, Section 4. |
| Compatible occupied counts on arbitrary physical rings | Same kinematics. Rounding to N in 16Z costs O(t*16^2/L). Actual selection remains conditional on a positive mobile dressed-defect cost. |
| Local ultraviolet normal form | Applicable at each fixed finite K, with K-dependent constants. It gives a spatially weighted o(J) remainder, not the missing low-frequency energy comparison. |
| Finite-cell massive matching | Must use the dressed domain-wall code and the full cell space, including illegal words. The d24 charge-based fermionic grading is not directly available on C75's odd-Q auxiliaries. Its strong-norm application is retained as an assumption. |
| Mobile defect comparison and final RG/CFT step | Remain assumed. Scalar contact susceptibility, finite-ring projection and the ultraviolet normal form do not prove them. |

[Proved here, using C75's all-length bound] The frozen sector argument can
be repaired without conserved bad-pair labels. Use complete marker-bounded
cells of dimension 3*12^(l-1). Only the port couples adjacent cells, with
norm B=sqrt(2)|lambda|. Put

    gamma=2^-16-1/(K_pen-4)>2^-17.

For a compatible material ring, a configuration with k wrong cells of
total length D_wrong obeys

    E_wrong-E_selected
       >=gamma D_wrong-3kB-2B
       >=(gamma-5B)D_wrong>0.

Marker-free rings are also excluded. This is a nonzero-port statement for
the full frozen model, including illegal configurations within its cells.
It strengthens C75's zero-port conclusion, subject to its inherited
all-length inequality. The derivation has not been independently refereed.
It does not establish the same comparison after coupling to vacancies.

[Numerical] All checks use resolved parameters, not the hierarchy above.

| Check | Size and outcome |
| :-- | :-- |
| Independent bond reconstruction | 256 by 256 lift, 437 nonzero entries. Zero entry error against C75, zero Hermiticity error and zero N,Q,M commutators. |
| Exact cyclic projection | Full L=4,N=2 space, dimension 1350; retained seas dimension 330. Projection error below 1.83e-14. |
| Finite-J ED against first order | All 120 lowest levels retained at J=.0008,.0004,.0002. Error/J^2 approaches 61.196; eigenpair residuals below 1.20e-14. |
| Literal two-cell vacancy rings | N=4,M=2,Q=2 at L=5,6, dimensions 2320 and 6960, with illegal strings retained. Six low eigenpairs at three J values; residuals below 1.49e-11. |
| Tagged-contact phase fluid | L=8,10,12, dimensions 35,126,462. Free energies and contact coefficients agree with exact formulas; defect relaxation has the expected J^2 coefficient. |
| Full-cell frozen selector control | Lengths (2,2) versus (1,3), dimension 1296 each. Positive energy differences at lambda=0,.005,.01. |

Reproduce with `python3 -B checks.py` and `python3 -B auxiliary_checks.py`
from this folder. Parameters, spectra, residuals and scope are recorded in
[checks.json](checks.json), [auxiliary_checks.json](auxiliary_checks.json)
and PROOF.md, Section 10.

[Assumed for c=16] The remaining inputs are finite-K massive matching in
the required coefficient norms, low-frequency vacancy matching, a uniform
common-bulk-subtracted mobile defect inequality, stable half-density
tuning, multichannel RG for the retarded endpoints, and sufficient
stress/current or scalar thermal control to identify sixteen conformal
modes. PROOF.md, Section 9 tags each dependency. The result is an explicit
smaller candidate with new frozen-sector control, not a completed
counterexample or a microscopic central-charge measurement.

All new files are inside work_C81. No external files were changed and no
network commands, cluster jobs or subagents were used. Local STATE.md
provides the handoff because the project STATE is outside the allowed scope.
