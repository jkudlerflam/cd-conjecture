# Proposed shared summary for the coordinator

An [independent bounded field-diagram audit](README.md) now covers
the fixed nearest-neighbour d=2 bond `(9/16)XY+(7/16)YX+ZZ` as h varies,
with full fields on every site. Exact symmetries relate h and -h and
identify staggered chirality as a same-parity, momentum-pi observable.
The real gauge, its periodic seam and the Jordan--Wigner boundary sign
are verified against dense Pauli matrices. Low fields show a Neel vacuum
doublet with a stable bulk gap; a conditional 1/L^2 extrapolation of the
first parity crossing suggests h about 1.668, with an unresolved working
search interval 1.66--1.70. High fields have a dressed gapped phase;
finite-h exact product saturation is impossible.

The saved iMPS ladders still make h about 1.7466 a boson x Ising candidate,
not a verified c=3/2 point. New chirality-carrying ring gaps at L=12,16,20
show substantial drift, so the tiny opposite-parity levels cannot be
counted as an Ising mode. The apparent c about 2 at h=1.7405 persists
under both parity-length estimators, but its extended-phase, isolated-point,
crossover and branch-artifact interpretations remain unresolved. A new
audit finds c_fit=1.626 for the supplied interacting boson x Ising control,
beyond the source's abbreviated calibration range. Long-distance target
correlators and convergence flags are absent from the committed tables.

The contribution supplies a phase table, an h-axis figure with uncertainty
intervals, exact checks, 50 completed ring points, seven operator-resolved
points, and reproducible reanalysis. Computation was one worker/thread:
about nine minutes total, including an enforced seven-minute stop on the
largest ring batch. No DMRG, cluster job, new agent or production restart
was launched. The recommended next test is an operator-resolved chirality
correlator and excitation calculation on recovered, verified existing
MPS states near h=1.7465. No new microscopic CFT theorem is claimed.

This file is proposed text only. Shared README, STATUS, REPORT and
coordination files were left to the coordinator, as instructed.
