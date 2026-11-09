# Microscopic transverse-DM chain: independent Davidson DMRG checks

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
