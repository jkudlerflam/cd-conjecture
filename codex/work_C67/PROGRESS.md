# C67 progress

2026-10-28

Read the project memory, C63 RESULT.md in full, C6's analytic and numerical report, the gap-reduction README, and the grid's stratum, ring-ED, open-chain, and finite-entanglement analysis code. The grid has no README in the local snapshot. Its calculations are in progress.

[proved, inherited algebra] The uniform-unitary C2 row has five real couplings before scale. Its strict internal-Z2 portion requires Jx != Jy and D b != 0. At D=0 or b=0 there are additional uniform antiunitaries. Spatially combined antiunitaries still exist in the strict portion.

[assumed] Exact-model CFT identifications and controlled continuum matching are inherited inputs. The charged-sector Sugawara exclusion from C63 does not apply to a generic chain with only fermion parity.

The new calculation tested the interacting Ising/chiral-liquid junction and all four Z2 entries in the grid candidate registry. Entanglement shape, size growth, and parity/momentum-resolved spectra were retained. A vanishing velocity was not assigned a conformal central charge.

Only files in this directory will be written. No network, cluster, or large-system calculation is planned. The project-level STATE.md is outside the writable scope; local STATE.md will carry the handoff.

## Completed work

[proved] Wrote the canonical parameter reduction, the direct absence-of-additive-charge check in the strict stratum, the explicit form of the imported gap cone, and the free cubic-junction expansion. Added the factorized XYZ-field construction and proved cancellation of the DM sum on both homogeneous product states. Its gapped small-D neighborhood is conditional on the stated standard parent-gap and stability inputs.

[numerical] Completed 70 point-size ED cases with L<=20, every momentum and both parities. Largest block dimension 26272. Largest retained eigenvector residual 9.3e-10. Dense independent checks give spectrum error 6.3e-15 and intertwining error 9.5e-15. Ising and XX ground energies agree with exact values through L=16 to 4.7e-14.

[numerical] All four flagged grid points fail to supply consistent high-c evidence. The strongest chord signal is Z2_21_295, c_chord(20)=2.259132, but fixed-mod-4 size slopes are 0.960061 and 0.867243. The inherited one-site FES runs for the two original survey points fail their residual filter and contain unphysical energies. The better 11-site-cell run for Z2_1_358 gives a preliminary slope 1.020698 through chi=24.

[numerical] Followed the Ising mass condition in (-1,-0.5,0.25,D,-g). One descendant difference collapses as D increases, becoming negative at D=0.24,L=20. Fixed-coupling runs at g=1.9891 distinguish the Ising-like finite-size entropy from a momentum-changing branch at D=0.30. The latter's c_chord(20)=1.193648 disagrees with its L=16 to 20 entropy step 0.229217.

[heuristic] The useful junction target is D near 0.20 and g near 1.989, where a slow branch appears to lose velocity. No extra finite-speed Ising sector was identified. The data do not establish the interacting exponent or an exact critical coordinate.

[assumed if used globally] Exhaustion by one-boson or two-Majorana continuum matching remains missing. RESULT.md explicitly leaves strong-coupling multicritical and irreducible c>1 possibilities open. The requested throughout-the-stratum argument has not been established, and no supported counterexample is reported.

RESULT.md, STATE.md, all scripts, raw ED records, and processed checks are saved locally. No file outside C67 was modified. No network commands, cluster jobs, subagents, or new tensor-network calculations were used.
