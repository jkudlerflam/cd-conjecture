# Narrow audit of the reciprocal-spectrum loophole

The task is to assess resonant designs with physical d<=10, without
confusing failure of C84's sufficient certificate with a universal bound.
Keep all new files here. No large numerical campaign or new DMRG.

1. Audit C84/C86/C87 hypotheses and the full Bloch symbol of the literal
   nine-letter single rail: four L/R bits, two carriers, three markers.
   It has n bits in n+2 physical sites, including the unused boundary bit
   position and marker. No independent vacancy is present.
2. Derive coefficients from that Hermitian nearest-neighbour bond.
   Do not choose an effective quartic interaction independently.
3. If the bare two/three-flavour symbol has only one Fermi pair, record
   the obstruction before discussing nonexistent pair-transfer channels.
4. Test the smallest explicit way to obtain two Fermi pairs without new
   letters: physical background fields and adjacent background-bit swaps.
   Verify their projections, then extract the first generated quartic
   interaction from the actual pair resolvent, including its signs.
5. Audit all momentum-allowed channels. Any weak-coupling stability
   conclusion must distinguish a scaling-dimension calculation from a
   complete RG or volume-uniform microscopic CFT theorem.
6. Preserve frozen markers, finite-K errors, length selection, and the
   distinction between one finite two-flavour seed and an all-n family.
   A full real symmetric two-mode onsite matrix need not remain realizable
   with the same small set of literal counterterms for arbitrary n.

Use small algebra and one worker with OPENBLAS_NUM_THREADS=1. Benchmark
before increasing matrix sizes. Prefer an exact narrow obstruction or
an explicitly qualified surviving seed to unsupported central-charge claims.
