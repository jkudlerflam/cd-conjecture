# Sector-resolved finite rings and the earlier DMRG audit

All new results here are **[numerical]**, on periodic L=8,12,16,20
with fixed g. No central charge is fitted from the small-ring gaps.
No finite-size trend is promoted to a phase-boundary theorem.

## Method and convergence

The exact real-basis blocks in EXACT.md retain Pi and F, each +/-1.
We solve (+,+),(+,-),(-,+); physical translation supplies the identical
(-,-) spectrum. The whole physical spectrum is reproduced by all four
blocks at L=4,6,8. An initial 12-level benchmark of all three distinct
blocks at g=1 took 0.005,0.020,0.201 seconds at L=8,12,16, respectively.
This justified stopping the larger calculation at L=20, where each
block has 262,144 states, and requesting 16 eigenvectors per block.

All reported ground states are in (+,+). Let Delta_N be the first
excitation in that same block, Delta_F the (+,-) ground energy minus
E0, and Delta_P the (-,+) ground energy minus E0. These are different
measurements, not interchangeable definitions of a mass gap. The
first neutral excitation has momentum pi relative to the ground state.
The Pi=-1 spectral duplication follows from symmetry even in an
unbroken phase.

Eigenvector residuals in the principal and neutral-weight scans are
below 6.3e-11 in the energy units of H. Translation is diagonalized
inside each returned energy cluster, with its leakage outside that
cluster recorded. The initial 16-vector windows omit partners of some
degenerate *excited* levels at g=1/8,L=12 and 16. Their momenta must not
be used. The follow-up `weak_resolution_check.json` requests 32 vectors
and checks these windows; it is the reference for their resolved
multiplets. The ground and first neutral levels are resolved already
in the initial run. A final returned cluster can be cut by the window;
neither an unresolved momentum nor an unreached level counts as absent.

## Gaps distinguish the rational targets from order

At g=1/2, L Delta_N at L=8,12,16,20 is
`8.1364, 8.2142, 8.2930, 8.3682`. At g=1 it is
`16.5732, 16.6473, 16.6565, 16.6528`.

| g | Delta_N, L=20 | Delta_F, L=20 | Delta_P, L=20 |
|---:|---:|---:|---:|
| 1/8 | 0.102487 | 1.228322 | 0.350043 |
| 1/2 | 0.418408 | 1.143766 | 0.489495 |
| 1 | 0.83264 | 0.88754 | 0.77261 |
| 2 | 2.330200 | 0.175901 | 2.227431 |
| 4 | 8.99241 | 0.00173 | 8.90852 |

The target points have spectra compatible with finite-size critical
scaling. At g=4 the F partner splitting falls as
`0.245810,0.043849,0.008517,0.00173`, while the other excitations stay
large. Calling that first splitting a critical gap would misidentify
the ordered control. At g=2 the partner also shrinks faster than 1/L
over this range, and L Delta_N rises; together with the order data
below this favors an ordered crossover, though it is not an extrapolated
proof or a located transition. Equal coefficients therefore do not
currently improve the critical candidate.

## Density and neutral dimer fluctuations

Use `M_Z(pi)=sum_j (-1)^j Z_j`. The table gives
`<M_Z(pi)^2>/L²`, which can detect order even though the mean vanishes
in the symmetric finite-size state.

| g | L=8 | L=12 | L=16 | L=20 |
|---:|---:|---:|---:|---:|
| 1/8 | .267658 | .179975 | .135553 | .108713 |
| 1/2 | .328422 | .228828 | .175754 | .142730 |
| 1 | .428782 | .323234 | .261800 | .221129 |
| 2 | .657233 | .608129 | .584944 | .572865 |
| 4 | .895736 | .895117 | .892543 | .889950 |

For bond types A=XY, B=YX, D=ZZ, and h=2A+gD, the stored connected
dimer quantity is `Var(sum_j (-1)^j O_j)/L²`. At g=1/2 the A/B values
fall from `.249926/.254748` to `.102268/.113523` between L=8 and 20;
at g=1 they fall from `.246537/.252909` to `.103111/.116009`.
ZZ and energy dimers, uniform density, and bond
fluctuations at 2pi/L are also saved. These measurements reveal no
finite-size dimer plateau at the target points. They cannot exclude
weak order beyond L=20. The g=4 density plateau is qualitatively
different and cannot be explained by a small velocity alone.

## Slow neutral mode: measured weight versus assumed velocity

The follow-up projects bond operators onto the returned neutral
eigenvectors. It reports the fraction of the **full** equal-time
spectral weight carried by each retained energy cluster; the denominator
does not discard unreached levels. At L=20 the first pi neutral pole has:

| g | v_N=L Delta_N/(2pi) | YX dimer weight fraction | XY dimer weight fraction |
|---:|---:|---:|---:|
| 1/8 | .326227 | .427611 | .005290 |
| 1/2 | 1.331835 | .536466 | .085528 |
| 1 | 2.650378 | .610585 | .245841 |

At g=1/8 this is predominantly the slow B bond, and v_N is close to
the first-order value 1/pi=.318310. This is concrete evidence for an
induced slow dispersion, not an extra gapless sector. At moderate g
we have not measured the scaling dimension x_N, so a CFT relation
would read `v_N=v x_N`; it would not by itself measure v. Nor does
the lowest 16-level window reach all fast-bond states. Several
2pi/L pole lists are empty for this reason; their full static norms
are nonzero. The saved operator-weighted mean excitation energies
are spectral averages and are not used as gaps or velocities.

## Audit of C144 and C157

The exact map `H(g)=2 H_C144(0,g/2)` makes the existing
[C144 alpha=0,beta=0.5 record](../../codex/work_C144/typhon_out/d128_b0.5_r0.jsonl)
directly relevant to g=1. It is an open L=128 calculation with
chi sequence `[16,32,64,96,128,128]`, last discarded weight
`5.73e-11`, and last energy change `-1.76e-11` in C144 units.
Its midpoint entropy is .914500. Separate even/odd cut fits give
`1.06392` and `.79866`; the earlier smoothed summary fit is `.920`.
These estimators are different, so their arithmetic average must
not be substituted for the smoothed value.

This supports the c=1 interpretation, but a single size with strong
cut-parity dependence does not certify c or a parameter interval.
The small sweep change is also not a size extrapolation. Earlier
alpha=0,beta=.15,.3 records correspond to g=.3,.6 and have smoothed
fits .951,.947. They are nearby evidence, not exactly g=1/2.
No new DMRG or claimed bond-dimension extrapolation is hidden in this audit.

[C157](../../codex/work_C157/RESULT.md) constrains the second-order
kernel near a *tuned* slow-chain cancellation and states finite-coupling
and higher-order caveats. Its attraction and dimerization calculation
is useful there. It supplies neither an all-order c bound on g=1/2,1
nor an exact solution of the whole two-term line. Pending remote jobs
mentioned in old reports are not treated as completed evidence.
