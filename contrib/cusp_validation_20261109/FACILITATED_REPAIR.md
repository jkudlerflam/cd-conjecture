# A microscopic repulsion knob and its unavoidable companion

AI assistance: GPT-6 Codex, 2026-11-09. Proposed facilitated perturbation by the parallel cusp-scattering agent; independent reduction and matrix verification here. This is a candidate modification, not a successful counterexample or a demonstrated removal of binding.

## 1. A genuine nearest-neighbour perturbation

[proved, algebra] Add `lambda F` to the microscopic cusp Hamiltonian, where

```
F = sum_j (n_j X_(j+1) - X_j n_(j+1))
  = (1/2) sum_j (X_j Z_(j+1) - Z_j X_(j+1))
  = sum_j (n_(j-1)-n_(j+1)) X_j .                         (1)
```

The second equality uses `n=(1-Z)/2`; the single-site terms telescope on a ring. It is a one-site-translation-invariant, Hermitian nearest-neighbour qubit bond, a transverse component of the Dzyaloshinskii–Moriya interaction. Its compression `P F P` vanishes exactly. An occupied site can facilitate creation only beside an occupied site, leaving the blockaded space. Annihilation would require two adjacent particles in the starting state.

Let `lambda=r Delta`. The leading perturbation out of the blockaded space is then `Delta(T+r F)`, with `T` as in [PROOF.md](PROOF.md). The leading second-order correction is

```
-(Delta^2/V) P (T+rF) R (T+rF) P .                       (2)
```

All cross terms in (2) are of the same perturbative order.

## 2. The complete additional second-order terms

[proved, algebra] For periodic rings with L >= 7,

```
-P F R F P
 = P [ sum_j (b†_(j+1)b_j + h.c.)
       -2 sum_j n_j +2 sum_j n_j n_(j+2) ] P,            (3)

-P (T R F + F R T) P
 = P sum_j (n_(j+2)-n_(j-2)) X_j P .                    (4)
```

Equation (3) supplies the intended distance-two repulsion `2 lambda^2/V`, together with a real nearest-neighbour hopping `lambda^2/V` and a chemical-potential shift. Equation (4) is an unavoidable density-assisted creation/annihilation term with coefficient `Delta lambda/V`. It changes N by one. Treating (3) as an independently tunable density interaction while omitting (4) is not the microscopic reduction of (1).

Proof of (3). In a blockaded configuration an empty site between two occupied sites has zero flip amplitude: its left and right facilitations cancel. Every other allowed virtual flip has amplitude plus or minus one, creates exactly one adjacent bond, and has denominator 1. The squared amplitudes sum to

```
sum_j (n_(j-1)-n_(j+1))^2 = 2N-2 sum_j n_j n_(j+2).
```

The minus sign of second-order perturbation gives the diagonal part of (3). A particle hop has path `10 -> 11 -> 01`; its two amplitudes have opposite signs, so the second-order hopping has coefficient plus one. These are the only off-diagonal paths returning to the blockaded space.

Proof of (4). Consider creation at an initially empty site j. The right adjacent-pair path creates `(j,j+1)` and removes `j+1`. Its facilitated removal amplitude is `1-n_(j+2)`; if `n_(j+2)=1` that amplitude cancels. The left path creates `(j-1,j)` and removes `j-1`, with amplitude `n_(j-2)-1`. Each nonzero path has denominator 1. Including the minus sign of the reduction leaves `n_(j+2)-n_(j-2)`. Reverse paths give annihilation. Creation from the vacuum cancels, but creation next to an existing particle at distance two survives.

[numerical] `check_facilitated.py` checks (1) projected into P and both complete matrices (3) and (4) against the direct full-spin reduction for every ring L=7 through 12. All entries agree exactly in the numerical calculation, with maximum absolute difference 0.0. Results are in `facilitated_results.json`. Runtime was below one second.

## 3. The tuned leading model

[proved, algebra] The unwanted real hopping and chemical-potential shift in (3) can be compensated by changing the allowed bare microscopic parameters. With `Delta=1` for readability, choose

```
w_bare = -lambda^2/V + i/(2V),
mu_bare = (3-2lambda^2)/V + mu_eff/(2V).
```

Then, in units `|t2|=Delta^2/(2V)`, the full leading constrained Hamiltonian is

```
H_r = P sum_j [ (i b†_(j+1)b_j + h.c.)
                -(b†_(j+2)b_j + h.c.)
                +(-2+4r^2) n_j n_(j+2)
                -(2/3) n_j n_(j+3) -mu_eff n_j
                +2r (n_(j+2)-n_(j-2)) X_j ] P.          (5)
```

Thus the induced repulsion is `v=4r^2`, but the simultaneous number-changing coefficient is `2r`. At the `v=1/4` threshold suggested by the fixed-two-body determinant, `r=1/4` and the extra coefficient is `1/2` in the same energy units. It cannot be omitted on grounds of powers of `Delta/V`.

Equation (5) still has the vacuum as an exact eigenstate, since the number-changing term contains a density difference. A one-particle state is not an invariant sector: the term creates a particle at distance two on either side with opposite signs. For a one-particle Bloch wave of momentum K, the amplitudes to a fixed distance-two pair are proportional to `exp(2iK)-1`, hence to `sin K`. They do not vanish at the incommensurate valleys.

[scope] At generic valley momenta, continuum charge-one vertices are off resonance. That does not make their virtual effects negligible compared with `v`: such effects begin at order `r^2`, the same order as the desired repulsion, and their denominators and form factors must be computed. The fixed-N two-body scattering problem for `H_0+v sum nn` is not the scattering problem for (5). It is therefore premature to claim that the microscopic perturbation crosses the binding threshold.

## 4. Reproduction and next calculation

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/cusp_validation_20261109/check_facilitated.py --max-L 12 --output contrib/cusp_validation_20261109/facilitated_results.json
```

A valid next calculation should use (5), including the density-assisted term, and determine its dressed one-particle dispersion and two-particle threshold, allowing particle-number mixing. If a stable two-mode region survives, it can then be mapped to the actual nearest-neighbour bond with the explicit bare tuning above. This audit establishes the exact leading model to test; it does not establish that such a region exists.

A bounded diagnostic, registered before execution, will compare complete-model ground energies and particle-number mixing on rings L=12, 16, 18 at `r=0, 0.2, 0.25, 0.3, 0.4`, `mu_eff=-2.2`. It retains every blockaded particle-number sector. The criterion is only to establish whether the omitted mixed term materially changes the finite-system state; these sizes and observables cannot establish or exclude two conformal modes.

After that first diagnostic, the coordinator requested an extension to L=20, `r=0.5`, and `mu_eff=-2.27`, just below the old model's dilute threshold. These additions test the same finite-system quantities. Six comparisons against the direct microscopic second-order matrix are required before trusting this implementation of (5).

## 5. Finite-ring diagnostic result

[numerical] `repair_ed.py` retains every nonvacuum blockaded basis state (15,126 states at L=20), and compares its lowest levels to the exactly decoupled vacuum with energy zero. Its complete leading Hamiltonian agrees with six direct second-order microscopic matrices at L=7, 8, 10 to `6.7e-15` or better. The 48 parameter/size runs are saved in `repair_ed_results.json`; runtime was about 7 seconds.

At `mu_eff=-2.27`, L=20:

| r | full ground energy | mean N | variance of N |
|---:|---:|---:|---:|
| 0 | 0, vacuum | 0 | 0 |
| 0.2 | -0.063206 | 2.5979 | 0.4669 |
| 0.25 | -0.120771 | 2.6097 | 0.4832 |
| 0.3 | -0.170413 | 2.5977 | 0.4823 |
| 0.4 | -0.400494 | 4.3910 | 0.9199 |
| 0.5 | -0.601993 | 4.2284 | 0.8730 |

The model with only the proposed positive repulsion `4r^2 sum nn` would remain in the vacuum on this ring: adding a positive semidefinite term cannot lower its nonvacuum levels below the original vacuum. The complete model instead has a negative-energy state already at `r=0.2`, below the proposed `r=0.25` binding threshold. Thus the mixed term changes the onset materially and cannot be assessed merely by moving `V2` through the fixed-N resonance determinant.

At `mu_eff=-2.2`, L=20, the original model has N=2 and energy -0.110013. The complete `r=0.3` model has mean N=4.7167, variance 0.9644 and energy -0.465763. The altered density is another reason to retune the chemical potential before comparing putative phases.

These are exact-diagonalization diagnostics of small rings, not a central-charge estimate, a thermodynamic compressibility calculation, or a proof of phase separation. The ring ground state changes particle number at small r even in the unmodified family, so interpolation of these sparse sizes would be misleading. The retained conclusion is the finite-model mismatch between a repulsion-only repair and the actual microscopic deformation.

Reproduce the complete diagnostic with

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/cusp_validation_20261109/repair_ed.py --max-L 20 --output contrib/cusp_validation_20261109/repair_ed_results.json
```

The coordinator subsequently requested a separate small-r diagnostic at `r=0.05,0.10,0.15`, `mu_eff=-2.30,-2.27`, with and without an additional real hopping `q r^2`, `q=-0.4746`. That tuning comes from an independent dilute self-energy calculation, and is only an approximate perturbative valley retuning. The new output will retain the first three levels' translation phases, particle numbers and variances, without interpreting them as conformal towers.

## 6. Small-r valley retuning: promising finite spectra, unresolved limit

[numerical; perturbative tuning supplied by the coordinator] To compensate the different atom self-energies, add an *extra effective* real NN hopping `q r^2` to (5). The bare microscopic real hopping is then

```
Re(w_bare) = -lambda^2/V + (Delta^2/(2V)) q r^2.
```

The coordinator's leading self-energy calculation gave `q=-0.4746` at `mu_eff=-2.27`, and `q=-0.36413` at `mu_eff=-2.30`. These are approximate dilute tunings, not a guarantee of exact degeneracy at finite r, density or ring length. `repair_ed.py` accepts the additional hopping explicitly.

The refined scans keep all particle numbers and compute the first three nonvacuum eigenstates. Their mean N and translation expectation identify two atom-like branches and one momentum-pi pair-like branch where all three occur among those levels. The quantity in the next table is their finite-ring energy combination `E_pi-E_a-E_b`. It is a diagnostic, not the infinite-line binding energy: atom momenta are quantized, the shallow pair wraps the ring, and the dressed atom minima may move.

At `mu_eff=-2.30`, using the preliminary `q=-0.4746` tuning:

| r | L=20 | L=24 |
|---:|---:|---:|
| 0.05 | -0.029195 | -0.035532 |
| 0.10 | -0.005284 | -0.019327 |
| 0.15 | +0.028537 | +0.002628 |

The sign progresses in the direction expected for weakening binding, but the positive r=0.15 margin shrinks by more than a factor of ten from L=20 to 24. It does **not** establish infinite-volume unbinding. These data are in `repair_ed_refined_results.json` and `repair_ed_refined_L24_results.json`. The L=24 scan took about 51 seconds, with 103,681 nonvacuum blockaded states per matrix.

The corrected `q=-0.36413` tuning at `mu_eff=-2.30`, L=24 gives:

| r | E_a | E_b | E_pi | E_pi-E_a-E_b |
|---:|---:|---:|---:|---:|
| 0 | 0.050311102 | 0.050311102 | 0.059136461 | -0.041485743 |
| 0.10 | 0.035377015 | 0.038369719 | 0.054307661 | -0.019439073 |
| 0.15 | 0.017353788 | 0.026228573 | 0.045793700 | +0.002211339 |

For nonzero r, the atom-like translation phases in the script's convention are `pi/12` and `11pi/12`, with translation magnitude 1; the pair-like phase is pi. At r=0 the numerical eigensolver mixes the degenerate atom states, so individual translation expectation phases are not momenta and have magnitude below 1. Their equal energies and N=1 remain well defined. The comparison is explicitly to the two finite-ring branches, not automatically to twice the globally lowest dressed-atom energy.

The corrected-tuning implementation passes 12 direct microscopic second-order matrix checks. Its output is `repair_ed_retuned_L24_results.json`, and the targeted run took about 11 seconds. Reproduce it with

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/cusp_validation_20261109/repair_ed.py --refined --max-L 24 --only-max-L --q-coefficient -0.36413 --r-values 0 .1 .15 --mu-values -2.30 --output contrib/cusp_validation_20261109/repair_ed_retuned_L24_results.json
```

The `--refined` default reproduces the L<=20 preliminary tuning scan. Add `--max-L 24` to reproduce its L=24 extension. The result supports computing the dressed few-body scattering problem for this *actual* microscopic deformation. Finite-ring sign changes, weak-coupling self-energy extrapolation and the availability of a nearest-neighbour bond are insufficient to claim a counterexample.
