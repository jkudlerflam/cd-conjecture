# Complete momentum scan — 2026-11-10, before its output

The three-momentum ED calculation treats k=0, even parity as the ground
branch, confirmed by previous full-space ED through L20. Complete the
lowest-energy scan over *all* momenta and both parities at L16,20,24.
This checks that assumption at L24 and records all candidate soft momenta,
without inferring a conformal mode count from a few selected levels.
An eigenvalue residual is retained for every block. Degeneracy within a
single momentum/parity block is not counted by a one-eigenvalue solve.

The finite-DMRG engine termination tolerance for new jobs is aligned to
1e-9 in total energy; the unchanged final four-record acceptance gates
still require an energy span below 1e-7, entropy span below 1e-5, local
residual below 1e-8, norm error below 1e-5, and independent energy agreement.
This avoids redundant sweeps below the registered acceptance precision.
