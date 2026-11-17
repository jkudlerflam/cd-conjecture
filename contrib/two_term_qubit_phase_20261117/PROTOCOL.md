# Modest two-term phase diagnostic

Target: `H(g)=sum[2XY+gZZ]`, g=1/2 and 1, fixed independently of length.
No new DMRG is planned. Use one CPU worker and one BLAS thread.

1. Check the real-basis transformation, exact symmetry blocks, physical
   translation and their multiplicities against complete small physical
   spin matrices. The odd physical-parity blocks have equal spectra;
   compare them before using that identity to save work.
2. Benchmark the full low-spectrum calculation at L=8 and L=12. Increase
   only if time and memory remain modest. Use compatible rings L divisible
   by four for the principal scaling, to avoid folding a twist into g->-g.
3. Retain all four exact internal-symmetry sectors (one duplicated by
   translation). Record same-sector neutral excitations separately from
   parity changes and from nearly degenerate ground partners. Resolve
   translation within returned degenerate eigenspaces and check residuals.
4. Record density order at momenta 0 and pi; staggered kinetic, slow-bond,
   ZZ-bond and total-energy fluctuations. Normalize squared extensive
   order parameters by L². A finite-size symmetric ground state can have
   zero mean order and nonzero order fluctuations.
5. Compare g=1/2,1 with a weak-coupling point and a strong Ising point.
   Interpret gaps using their sector, momentum, size dependence and
   operator weight. Do not infer C from a small gap or from the number of
   microscopic Majoranas. Do not regard an unreached eigenvalue as absent.

The analytic audit separately checks faithful term-by-term free duality,
the standard regular difference-form Yang–Baxter charge condition, weak
coupling generated hopping, and the limitations of a large-|g| expansion.
Neither a failed integrability test nor small-ring evidence is a global
interacting central-charge theorem.
