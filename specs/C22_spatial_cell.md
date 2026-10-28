# C22: A route-A cell without a history-state computation

Read `../../specs/_routeA_background.md` first. If `PROGRESS.md` exists in your folder, continue from it.

**Your part.** Find a different way to build route A's cells: a spatial, computation-free design in which a cell of length growing with n holds n logical fermion flavours directly in its spatial structure (for example n repeated sub-blocks, a "bus" or wire that carries virtual excitations between the ports and the flavour sites, or a frustration-free ground space with a degenerate code), so that second-order (or higher-order) exchange through small ports gives independent flavour hopping with distinct phases and chemical potentials, and exact logical parity. The point is to avoid a universal history machine and the controller problem entirely.

Tasks.
1. Propose concrete cell designs built from one repeated nearest-neighbour bond on a small alphabet, plus a cell-length selector (C7 or C18). For each, compute the code (ground space of a cell), the port compressions (they must vanish at first order for route A, or give first-order Majoranas for a route-B variant), and the second-order exchange.
2. Respect the pole-count bound (2n ≤ 4(k − 1) for k accessible energies): show where the n + 1 distinct energies come from (spatial position, distances to the ports, sub-block sizes) and how distinct hopping phases and chemical potentials arise.
3. Work out the scale hierarchy: cell gap, port visibility, logical bandwidth, selection margin, and the conditions for the logical chain to sit in the stable incommensurate region of `../../contrib/ir_control/` and `../work_C20/`.
4. Verify the best design by exact diagonalization of two and three cells (n = 1, 2), and give d0. If every spatial design fails, state the obstruction precisely and whether it is specific to second order.
