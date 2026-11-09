# Independent d=2 priority audit, 2026-11-09

AI-only review by OpenAI Codex (GPT-6 family). Repository snapshot for the
collaborator results: `524aaeda0821d152edfe0f5e6778d140a040f3f0`. The parent
agent handles Git synchronization and new microscopic DMRG. This folder
contains lightweight independent algebra and analysis, not new large DMRG.

**The full microscopic tests have not confirmed c>1.** The transverse-DM
points read near one. The simpler parity family initially read about two
at 32→64, but its L=128 follow-up falls to 0.738–0.931 across every parity
pairing and remains incompletely optimized. The positive effective-model
data below must not be substituted for that negative microscopic size
check. The most informative next calculation is now a fixed-Hamiltonian
solver/χ convergence study, specified in [NEXT_STEPS.md](NEXT_STEPS.md).
The latest full-chain evidence is reviewed [separately](../d2_tenpy_20261109/REVIEW.md).

The new analytic correction is in [UNIFORM_Y_SOURCE.md](UNIFORM_Y_SOURCE.md):
the surviving antiunitary permits a uniform Y source, while C183 checked
only X. Exact finite-ring canonical perturbation theory finds the first
source at order V^-5 in the cusp scaling. At the C183 V=16 target its leading
coefficient is only about 2.54e-7. This corrects an all-orders symmetry claim;
it does not by itself close the route or establish the need for a counterterm.

## Best next targets, without duplicating queued work

| Priority | Target | Why it is useful | What still decides it |
|---|---|---|---|
| 1 | Strict solver and χ escalation for the **same** microscopic C181153 L=128 bond, both parities | Resolves whether the failed larger-size check is physical crossover or optimizer bias | Converged fixed-H states, stable entropy across χ/starts, and the matched L=64 comparison |
| Supporting check | Existing rounded parity-family independent-solver work | Tests code independence and simplification | Its own local residuals and larger sizes; it cannot override the different unrounded bond's L=128 result |
| Conditional follow-up | Larger V with matched effective parameters and actual density | Tests finite-V effects after numerical conditioning is understood | Validated solver and physical-field matching; no automatic transfer of the effective-model c≈2 reading |
| Deprioritize for now | Further transverse-DM grids or the near-balanced φ=1.555 point | Existing transverse-DM full-chain readings are near one; φ=1.555 has stronger shell ambiguity | A specific new mechanism or evidence, rather than selecting a larger finite-size fit |
| Deprioritize | Symmetric parity cusp or far-vacuum V=8 fitted coefficients | Both can show c≈2 before their pairing scale is resolved; the V=8 “unbinding” fit used the wrong chemical-potential regime | They are controls for crossover, not decisive first targets |

For reference, the tested C183 Pauli parameters are

```
Jx = 0.4979983, Jy = -0.5020017, Dz = 0.0152122, Jz = 4
Gamma_xz = -Gamma_zx = 0.13693064, hz near -7.9463
```

The parent agent's calibrated microscopic runs at this field read close
to one, separately from C186's second-order unbinding prediction. The
original V=40 transverse-DM point has decreasing central entropy and
density with L; it is not a demonstrated finite-density critical target.

The leading V=40 parity bond at \(\phi=1.53\) has
\(J_x-J_y=1\), \(J_x+J_y=\cos(1.53)/80\),
\(D_z=\sin(1.53)/160\), \(J_z=10\), and no transverse DM. The independent
L=128 finite-difference chemical potential is -2.1518497755, giving the
**leading-order starting field** \(h_z=-19.9759490611\). Microscopic and
finite-size corrections still require density tuning. Larger V is motivated
by C181/C179's analysis of the number-changing two-pair vertex; this audit
does not certify that vertex's thermodynamic irrelevance. The actual
microscopic runs instead used hz=-19.97589125. **Keep that actual field
fixed for the convergence diagnostic**; the effective-model starting-field
estimate above is not a reason to change the Hamiltonian during that test.

Work already registered upstream should be collected before launching
duplicates: C186's effective transverse-DM DMRG jobs 17604628 and 17604527;
C183's full-chain parity/generic screen jobs 17604545 and 17604499; C187q's
local L=192 off-symmetric continuation and Typhon L=256 job 17604910.
These IDs describe the repository's recorded jobs, not a fresh live queue
query by this agent.

## Independent reading of the C181/C187 raw data

`audit_existing_dmrg.py` recomputes the estimators independently from the
committed JSONL records. `existing_dmrg_audit.json` records their hashes,
actual model settings, sweep completion, entropy sensitivity and gaps.
Fixed windows use cuts \(|l-L/2|\le W\), W=1,6,12, exactly as C181 registered.

| Effective model | c_two(64→128), W=1,6,12 | c_profile(128), standard fit | polarization 64→128 | L·parity-gap 64→128 |
|---|---|---:|---|---|
| c=1 control | 1.0609, 1.0088, 1.0454 | 0.9622 | — | — |
| c=2 control | 2.1538, 1.8562, 2.0817 | 1.9269 | — | — |
| φ=1.53 | 1.6546, 1.8336, 1.9077 | 2.1408 | 0.1163→0.1323 | 2.2385→1.5270 |
| φ=1.555 | 2.1804, 2.0533, 2.1804 | 2.1207 | 0.0087→0.0651 | -0.2232→1.5001 |

Both controls pass the stated calibration at 64→128. At L=128 the candidate
profile fits remain about 2.0–2.18 when varying smoothing (none, 8,12,16 cuts)
and fit interval (outer sixth removed through outer third removed). The
c>1 reading is not an artifact of a single profile window. It is still a
finite-size effective-model result.

The final candidate discarded weights are 8.13e-8 and 8.49e-8, just below
the registered 1e-7 threshold. Raising χ from 256 to 384 increases the W=6
entropy by about 0.0069, corresponding to about 0.060 in the 64→128 c
estimate if the smaller size is fixed. The final same-χ sweep changes it by
about 0.00052, corresponding to 0.0045 in c. These observed convergence
changes are small compared with the separation from c=1, but are not rigorous
error bars or a substitute for higher χ at larger L.

The density is homogeneous but not size-matched in the bulk: at φ=1.53 its
central-half mean is 0.21939 at L=64 and 0.20371 at L=128, although N/L is
fixed to 3/16. At L=128 the 16-site-smoothed central profile varies by less
than 0.8%. This supports a fluid rather than a droplet at these sizes, while
the boundary-induced density drift remains a systematic error in c_two.

Within these effective-model data, φ=1.53 is the cleaner comparison. At
φ=1.555 the negative L=64 parity gap
and near-zero polarization show a shell crossing. C181's literal registered
condition that L·gap grow by less than 50% has no meaningful positive
baseline at that point; the broad “first-stage pass” language in REPORT
should be read with this exception, not as a literal pass of every original
inequality. φ=1.53 does not have that ambiguity.

The C187q L=192 candidates in this snapshot are incomplete: only χ=64 has
been reached and discarded weight is about 1.7e-5. Their entropy estimates
and parity gaps must not be used for a phase verdict. The old L=192 c=1
control is also a different Hamiltonian: it explicitly set t2=0, while the
C181 reference uses the default t2=-1. Its apparent 128→192 c≈0.44 is
meaningless. Commit 524aaed records this error and introduces `ctl1fix`
configs restoring t2=-1 at L=192 and 256. Compare actual configs, not labels.

## Scope ledger for the newest collaborators

- **C186:** The general second-order nearest-neighbour reduction and the
  transverse-DM number-changing vertex agree with our earlier independent
  local-algebra check. At u=1.035, v=0.1 and 0.3 its finite-ring,
  particle-cutoff-converged calculation resolves no bound pair or cluster
  below the dressed atom onset. This is promising few-body evidence; it
  does not compute the finite-density pseudospin theory including the
  number-changing vertex or establish a microscopic c=2 phase.
- **C177/C183/C186:** The old V=8 parity-class “unbound” fitted Hamiltonian
  was extracted far below the gas onset. Virtual denominators depend on
  chemical potential. Near onset the microscopic parity chain again shows
  binding. An apparently unbound far-vacuum fit must not be used as a gas
  Hamiltonian without re-extraction at the physical field.
- **C181:** Strongest current effective-model finite-size evidence is the
  imbalanced |w|=1 fluid at L=128. The earlier negative 64→96 preview and
  stale RESULT.md verdict are superseded by the committed backfill data
  for this effective-model purpose. The microscopic chain has now been
  tested separately: its initial 32→64 excess does not persist in the
  unfinished L=128 optimization. Solver/χ convergence, the FFLO momentum
  shift, larger effective-model sizes, FES and finite-V vertices remain
  unresolved; there is no confirmed microscopic c≈2 phase.
- **C178/C187 symmetric line:** A 1/L parity-like gap does not discriminate
  when the paired control also passes. The L≤144 entropy and parity data
  do not resolve the weak pair gap. This explains why the symmetric line
  is a poor route to a quick decisive numerical answer.
- **C182:** Large entropy slopes in a droplet or phase-separated state are
  not central-charge evidence. Its |w|≥2 cluster regime is not the best
  next target. The |w|=1 off-symmetric profiles above pass its homogeneity
  check at L=128.
- **This contribution:** Corrects the uniform-Y symmetry inference and
  supplies exact finite-ring source coefficients plus reproducible raw-data
  audits. It supplies no thermodynamic phase theorem and no d=2
  counterexample claim.

Parent-agent microscopic results, including its independently validated
MPO and new finite/infinite DMRG, are separate contributions. No human
physics review is claimed for this folder.
