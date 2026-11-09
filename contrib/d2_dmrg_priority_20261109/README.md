# Microscopic qubit chains: independent Davidson DMRG checks

AI assistance: GPT-6 Codex, 2026-11-09. This work targets the literal nearest-neighbour qubit Hamiltonian, not a projected constrained model.

## Registered protocol

Use C183's ED-validated two-site Davidson DMRG implementation, unchanged, with the physical uniform Pauli field on every site. In open chains the missing half field is added at each end. No transverse edge compensation is added: the model is the physical Pauli bond, so boundary telescoping terms are not silently removed.

First target is C183's retuned V/Delta=16, v=0.3 transverse-DM bond:

```
Jx=0.4979983, Jy=-0.5020017, Dz=0.0152122, Jz=4,
Gamma_xz=-Gamma_zx=sqrt(0.3)/4,
hz around -7.9463.
```

1. Check small L against independently assembled sparse exact diagonalization, including energy, density and entropy.
2. Locate a gas window using an L=24 field ladder; report raw density and boundary profiles rather than treating dressed-vacuum density as gas density.
3. Continue useful points to L=32,48 and, if feasible,64 with increasing bond dimension and repeated sweeps. Record local eigensolver residuals, energy convergence, truncation and entropy profiles.
4. Calibrate with the XX chain and two independent XX chains. A finite-size entropy reading is a diagnostic only; a credible counterexample needs larger-size convergence, an independent estimate, and a check for a long crossover.

No c>1 pass is registered for these short chains. The screen is favourable only if profiles and two-size entropy slopes remain above 1.2, controls distinguish c=1 from c=2, final local residuals are below 1e-8, and last-sweep energy movement and discarded weight are reported. Any failures remain visible.

The first L=12 result is independently ED-validated. It has substantial uniform Y response and low entropy; X is zero as required by the surviving antiunitary symmetry. A registered L=24 ladder now uses hz=-7.9500,-7.9480,-7.9470,-7.9463,-7.9455,-7.9450,-7.9440,-7.9420,-7.9400 with chi=16,32,48,48,64,64. This deliberately brackets the narrow proposed gas window and possible packed-phase jump. It is a density search, not a c estimate.

After the ladder output through hz=-7.9440, two follow-ups were selected and recorded before their L=48 output: hz=-7.9450 (raw n_mid≈0.17 at L=24) and -7.9440 (first profile fit above 1.2). They use chi=16,32,48,64,64,96,96. Matching the same fields at 24 and48 tests whether the profile excess survives a size increase. Same-code XX and two-XX-leg controls use L=24 and48; exact correlations predict central-cut two-size slopes 1.06374 and2.12030 respectively. The 48-to64 low-density control slopes have severe particle-number aliasing (1.756 and4.239), so that size pair will not be used as an uncalibrated c estimate.

The coordinator requested an actual-density valley-tuning scan, since C183's real-hopping tuning came from the far-field band. Registered before its output: at L=24, hz=-7.9450, fixed Jx-Jy=1 and all other couplings fixed, scan Jx+Jy=-0.0065,-0.0060,-0.0055,-0.00525,-0.0050,-0.0045,-0.0040034,-0.0035. Record density changes along with entropy; a larger entropy at a different density is not automatically a better c>1 candidate. This runs after the matched L=48 tests.

Amendment before any real-hopping-scan output: fix hz=-7.9463 instead of -7.9450. The matched L=48 run at -7.9450 reached raw central density0.2397, whereas the coordinator's independent L=64 run at -7.9463 reached0.1716. This uses the longer-chain density information to target the intended region and avoids drifting toward quarter filling. The eight Jx+Jy values and convergence schedule remain as registered.

After the real-hopping scan, the coordinator authorized one density-matched follow-up around its zero hopping response: fixed Jx+Jy=-0.0058, L=24, hz=-7.9460,-7.9455,-7.9450, same convergence schedule. The wrapper now measures <XX+YY> directly and saves a local MPS checkpoint in /private/tmp/cd-conjecture-mps. At most one new L=48 point will be selected if the density-matched data warrant it; no unlimited parameter grid is planned.

One final matched-size point selected after the three-field follow-up: Jx+Jy=-0.0058, hz=-7.9450. Its L=24 raw density0.1679 is close to the original matched point's0.1728, with Sbar higher by0.0255. Run L=48 at chi=16,32,48,64,64,96,96, local residual target1e-8 (the registered physics tolerance), and report all residuals and last-sweep energy changes. This is a bounded check of actual-gas real-hopping retuning, not selection by an uncalibrated c fit. No further grid is scheduled here.

## New priority: rational parity-class candidate

After the transverse-DM tests, the coordinator reported substantially stronger c≈2 finite-size readings in the simpler partially polarized C181 family and requested an independent solver and rational-coefficient simplification. This is a new, explicitly registered stage. Use the literal qubit Hamiltonian

```
H = sum_j [10 Z_j Z_(j+1) + (X_j X_(j+1)-Y_j Y_(j+1))/2
           +(X_j X_(j+1)+Y_j Y_(j+1))/4000
           +(X_j Y_(j+1)-Y_j X_(j+1))/160
           -(199759/10000) Z_j].
```

It has no transverse DM term and conserves parity. First validate both parity selectors against sparse ED at L=10. Then run L=32 and64 in each parity using the independent Davidson solver, chi96 and128 respectively, local residual target1e-8, repeated final sweeps. Measure the selected parity directly and remove the selector energy using its measured expectation. Compare the lower-energy parity at each length; never form a ground-state size slope by silently keeping the same parity when the minimum switches. The exact coefficient simplification and the new sizes are set before this stage's output. A two-size reading alone remains insufficient for a thermodynamic CFT claim.

Optional continuation uses `full_chain.py --resume CHECKPOINT` and the local `dmrg_restart.py` copy. Only optional MPS initialization and max-iteration-hit diagnostics are added to the C183 solver; the numerical updates are unchanged. The restart path is separately checked against exact L=10 parity-sector results. Checkpoints live outside Git in `/private/tmp/cd-conjecture-mps`, while observables and logs are committed.

Before any rational L=64 entropy output (while its final two chi128 sweeps run), register a secondary shell-smoothing diagnostic: integrate the piecewise-linear entropy profile over cut coordinates l in [3L/8,5L/8], divide by L/4, and form the same six-times-log-size slope. Endpoint half weights implement the trapezoid when boundaries are integer. A central-half [L/4,3L/4] average is also retained for window sensitivity. Preserve the primary central-three-cut estimator and do not choose the largest result. For an exact CFT profile the fixed-fraction geometric average contributes an L-independent constant, apart from quadrature error. No additional optimization is involved.

Final rational stage status: both parities at L=32,64 and the bounded L=32 continuations are complete. The lower computed branches give central-three-cut c_eff=2.12889, but the registered 1e-8 local residual target fails (maxima around1e-6 at32 and1e-5 at64). The nearby unrounded128 data drift toward c≈1. See SIMPLE_PARITY_CANDIDATE.md for both the positive intermediate-size signal and its limitations. No further optimization remains running.
