# Current qubit c equals two candidates

2026-11-23. OpenAI Codex, GPT-6 family. **[Repository status review;
numerical table checks, not independent many-body reproduction.]** Reviewed
GitHub main through c3ef5408, the latest issues and the sole open pull
request. The latter concerns a different bound program and contains no
new qubit phase claim. No new simulation or production job was launched.

The repository now has stronger numerical evidence for c approximately
two in a simple four-term qubit family, alongside the older large-blockade
family. Neither has a settled thermodynamic compact unitary Virasoro c=2
identification. The previously reported c=3/2 chirality point was withdrawn
on November 19; historical November 18 summaries must be read with that
correction.

## The simple chirality family is the most direct new numerical target

With physical d=2 and Pauli matrices,

    H(delta,h) = sum_j [((1+delta)/2) X_j Y_(j+1)
                       + ((1-delta)/2) Y_j X_(j+1)
                       + Z_j Z_(j+1)] - h sum_j Z_j.

All terms are nearest neighbour or on site, with one-site translation
invariance. The field is counted once per site. The
[updated contribution](../d2_xy_yx_chirality_imps_20261118/README.md)
and [REPORT sections 9.209b and 9.212](../../REPORT.md) give:

| Point | Infinite-chain entropy fit | Open-chain entropy fit |
| --- | --- | --- |
| delta=3/16, h=1.783 | 1.973 +/- 0.010 | reported 1.96 +/- 0.08 |
| delta=1/8, h=1.7405 | 1.919 +/- 0.009 | reported 1.90 +/- 0.12 |

The open-chain fits use L=96 to 384 with bond dimension up to 1024.
The two methods share the same TeNPy model code. They are complementary
methods, not an independent second-code reproduction. Model energies
were separately compared with small-system exact diagonalization.

The bounded table check here reconstructs the prescribed central-window
chord coordinate and applies the documented final-entropy-step filter
abs(delta Sbar)<=0.02. From the rounded committed tables it gives
1.9540 +/- 0.0750 and 1.8966 +/- 0.1236, consistent with the published
values. The nearby one-boson point delta=3/16,h=1.772 gives
0.9984 +/- 0.0104. These errors are regression errors; they do not include
the full finite-size, fit-window and bond-dimension uncertainty.

In particular, the target L=512 entries fail the entropy convergence
filter and must not be included as converged c measurements. At L=384,
the final entropy step is still 0.011 or 0.017. Removing the smallest
included length moves the published fits appreciably, and local slopes
oscillate by about 0.5. A much longer crossover is not excluded.

The numerical phase interpretation is a c approximately two wedge, with
an indicative window 1.780 less than or approximately h < 1.7875 at
delta=3/16. It narrows near delta=7/32 and is absent on the tested
delta=1/4 cut. These are numerical boundary estimates, not rigorous
intervals. Two nearby critical lines and finite-scale crossover remain
alternatives to an extended phase.

## New correlators support one stability test but do not settle all modes

[REPORT section 9.214](../../REPORT.md), added November 23, reports
K_s=1.77(12),1.87(21),1.97(18) at the delta=3/16 point for
L=192,256,384, using the earlier pair and one-body estimators. Within
the proposed two-valley theory, this makes pseudospin backscattering
irrelevant. The target's charge estimates disagree between pair and
structure-factor methods, and known boundary effects are not removed.
At delta=1/8 the pair estimate drifts too much to determine K_s.

Consequently these measurements support the two-valley interpretation;
they do not establish two independently gapless modes or exclude every
parity-allowed locking interaction. The suggested pair-binding lower
edge is also unsettled: the one-body correlator stays algebraic over the
observed range, and the open-chain odd-even energy difference changes
sign with length.

The latest below-tuning scans at delta=1/16 and 3/32 find gapped regions
through about h=1.725. A possible c approximately two sliver below
delta=1/8 remains unresolved. Failed infinite-chain runs that selected
the wrong high-energy branch are explicitly excluded by the source.

## The older large-blockade candidates remain viable and conditional

The [STATUS qubit section](../../STATUS.md) and REPORT section 9.203
retain the symmetric V/Delta=40 example and the asymmetric phase-1.53
p153c neighbour. Their open-chain readings from L=256 to 512 are
1.983 and 1.954 respectively at bond dimension 2048, compared with
1.95 on an exact c=2 control. Their last entropy changes slightly exceed
the registered convergence threshold. The near-two infinite-chain
results for the projected effective model are not measurements of the
full microscopic qubit Hamiltonian.

The symmetric model permits number-changing quartet locking. The
asymmetric deformation may suppress it through momentum mismatch, but
the [connected response calculation](../qubit_detuned_quartet_20261118/README.md)
does not yet supply the continuum coupling, and the
[dressed-medium audit](../qubit_finite_density_modes_20261118/README.md)
does not yet give the thermodynamic response and mismatch needed for a
stability verdict. No new November 19--23 result closes those gaps.

## The two parameter regimes are related by a simple rotation

REPORT section 9.213 identifies an on-site rotation mapping the chirality
family to

    H = sum [(XX-YY)/2 + (delta/2)(XY-YX) + ZZ] - h sum Z.

The check here verifies this bond identity directly to floating-point
precision, with Z and hence the field unchanged. It lies in the same
broader pairing-plus-DM family as the older symmetric candidate, at
V/Delta=4 rather than 40. This exact algebraic relation does not establish
a continuous path within one phase or identify the low-energy theory at
stronger coupling.

For the user's simplicity preference, the delta=3/16 point is especially
attractive: its displayed bond coefficients are 19/32,13/32,1, with field
1.783. It has four literal Pauli monomials and coefficient spread 4.389
in this basis and full-on-site-field convention, after fixing the common
scale by ZZ. The delta=1/8 point has spread 3.978. The field window is
still narrow; these observations do not demonstrate a broad robust phase.

Our current assessment is therefore: prioritize the simple delta=3/16,
h=1.783 point for further validation, while retaining the phase-1.53
candidate as a distinct parameter regime with its own stability tests.
The important remaining evidence is controlled longer-distance scaling,
an independent implementation, and identification of the two modes and
all allowed locking processes. No new universal bound on c is implied.

## Reproducing this review's bounded checks

Run `python3 contrib/qubit_c2_status_20261123/check_tables.py` from the
repository root. It uses only the standard library, reads the committed
rounded tables, checks the local rotation and writes [checks.json](checks.json).
This is a reanalysis of saved observables, not a rerun of DMRG or iDMRG.
