# Bounded field-diagram audit

2026-11-18. One worker, one BLAS thread. No DMRG, cluster jobs, detached
processes, resumable production, or additional agents. The fixed bond is
9/16 XY + 7/16 YX + ZZ and the complete field is -h Z at every site.

1. Validate literal Pauli matrices against a separately assembled real-frame
   matrix for open chains and even rings, including L=2 mod 4 seam twists.
   Check parity, translation, antiunitary symmetries and observable mappings.
2. Coarse survey L=12,16, h=0,.5,1,1.4,1.6,1.7,1.738,1.7466,1.755,1.8,
   1.9,2,2.2,2.5,3,4,6. At most 360 seconds for this batch, with an alarm
   enforcing the limit. Six eigenpairs per parity, residual below 1e-8.
   Measure unrestricted wavevector structure factors for Z, the literal
   bond, distance-two XX-YY, and distance-two chirality; retain full
   connected translation-averaged correlators. Resolve translation within
   each returned degenerate subspace before reporting ground observables.
3. Adaptive follow-up may use L<=20 and selected fields, with at most a
   further 600 worker-seconds total. Choose fields from the coarse results;
   log the choice before execution. Do not fit critical exponents from
   these rings if the inferred pitch exceeds the system length.
4. Refit saved target/control ladders at chi>=64,96,128,160 and on the last
   four points, requiring at least three points and increasing xi.
   Keep branches separate, compare xi0 and xi1, record cat-sector changes,
   and do not treat least-squares residuals as physical error bars.
   Convergence flags, transfer-eigenvalue phases and MPS correlators are
   unavailable in the committed tables and cannot be reconstructed.

No thermodynamic phase boundary is to be manufactured from a coarse grid.
Finite-ring crossover brackets and inherited FES estimates will be labeled
separately. Any extended liquid or ordered interval inferred by continuity
will be identified as an inference rather than a certified interval.

## Adaptive decisions

The completed coarse survey (about one worker-minute) finds a Neel doublet
at low h, severe size-dependent crossings near 1.7--1.76, and rapidly
size-convergent high-field spectra above 1.9. Follow-up 1 is L=12,16 at
h=1.64,1.66,1.68,1.76,1.77,1.78 (cap 90 seconds). Follow-up 2 is L=20
at h=0,1.66,1.68,1.7466,1.77,1.8 (cap 420 seconds). These choices test
the low-field onset, one proposed Ising point, and departure to the
high-field region; they are not a dense fine scan. All follow-ups together
remain subject to the 600-second cap registered above.

The L=20 batch hit its 420-second alarm after completing h=0,1.66,1.68,
1.7466. Its h=1.77 calculation was interrupted before a record completed;
h=1.8 was not started. No incomplete eigenpair was used and this batch
will not be resumed. The saved JSON has `complete: false`.

Follow-up 3: resolve momentum blocks at L=12,16 and h=1.738,1.7466,1.755,
using at most 90 seconds, one thread. Measure the K2 spectral weight at
momentum shifted by pi within the ground parity, and comparison channels
at shifts 0,+/-2pi/L. This explicitly targets the relevant operator instead
of counting generic low eigenvalues. Also retain parity-changing transverse
spin structure factors. Together with follow-up 1 (about 25 seconds) and
the interrupted follow-up 2, this stays below 600 follow-up worker-seconds.

The operator-block test completed in 4.9 seconds. At h=1.7466 its
chirality-carrying pi gaps are 0.08269 (L=12) and 0.03742 (L=16), with
strong drift in L times the gap. Follow-up 4 uses the much smaller
momentum blocks at **L=20, h=1.7466 only**, cap 90 seconds, to test this
drift at a third size. It reuses the known ground parity/momentum from
the completed full-ring record. No restart of the timed-out h=1.77
calculation is authorized or planned. Total worst-case follow-up time
is now below 550 seconds.
