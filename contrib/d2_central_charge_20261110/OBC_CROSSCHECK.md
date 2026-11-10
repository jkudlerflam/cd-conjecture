# Same-boundary implementation check — 2026-11-10 18:50 UTC

The open-chain entropy estimates came from the prior independent C190
solver, whereas this continuation uses TeNPy. Before attributing their
different scaling trends to physics, repeat the identical open L48 even
state at chi128 in TeNPy and compare the energy and fractional-window
entropy with the pinned previous output. This checks that both solvers
represent the same couplings, boundary convention and entropy normalization.
The only change to the current runner is bc_x=open and the open-chain
profile prefactor. Existing periodic and infinite records are unchanged.
