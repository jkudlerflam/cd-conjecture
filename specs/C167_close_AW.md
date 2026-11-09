# C167: Close (A_W), the last assumption behind c ≤ 1 on the time-reversal strata T and K_I

Your folder is `../work_C167/`. Read `specs/_qubit_c_context.md` first. Then read:
- `codex/work_C163/RESULT.md` and PROOF.md, where (A_W) is stated, with route (a) and its energetic heuristic.
- `codex/work_C159/RESULT.md` and PROOF.md: the attractive window, half-floating liquids, Lemma A, the Ising × boson counting, and A6_TK, which A_W replaced.
- `codex/work_C161/REFEREE.md`.
- `codex/work_C160/PROOF.md` sections 1 to 3, the semiclassical corner near h_F.
- `codex/work_C162/RESULT.md`, immiscible two-valley gases.
- `codex/work_C164/RESULT.md` and `typhon_out/summary.txt`, the null soliton-onset search.
- `codex/work_C165/RESULT.md` and PROOF.md, especially §4: the spin-1/2 SU(2)_1 obstruction to a staggered Ising mode inside a liquid, and the residual A_S.
- REPORT.md §9.120 to §9.126.

## Why

c ≤ 1 holds at the physics standard on every symmetric class of nearest-neighbour qubit chains, except that on T and K_I it is conditional on one assumption. Closing it leaves only symmetry-free bonds open.

(A_W) [assumed]. Take an attractive half-floating liquid of T or K_I (K > 1), with harmonics comparable to the exchange. Then no lattice-scale excitation outside the boson's charge lattice reaches zero energy while the liquid survives. Concretely:
- (a) no independent staggered Ising mode at K > 15/8;
- (b) no composite band minimum off the charge lattice at K > 1.

It bites for parents with −1 < Δ < −0.67 (case a) or −1 < Δ < 0 (case b), at harmonic fraction above about 0.3, away from h_F and the twisted curve γ. Numerically it holds everywhere tested: 16 attractive points with K_S(16) from 1.9 to 7.5 read c_prof(128) from 0.98 to 1.03.

The task is a derivation of (a) and (b), or the narrowest form of what cannot be derived, or an explicit bond that violates it. This is the target theorem. Attack only this lemma, not side cases.

## Leads to try (your judgment, not a checklist)

1. **Case (a).** C165 showed that at spin 1/2 the staggered order and the spiral sit in one SU(2)_1 primary near the hidden Heisenberg point of the U(1) parent, so a separate Ising mode cannot appear there. Ask whether an analogous identification holds for the attractive parent (−1 < Δ < −0.67, near the ferromagnet). There the natural infrared object is the ferromagnetic magnon or the SU(2)-broken multiplet, not SU(2)_1. Is the would-be staggered mode a composite of the liquid's own vertex operators, with a dimension forced above 2 for K > 15/8? Use the exact site rotations that map T and K_I parents to XXZ forms.
2. **Case (b).** A composite band minimum off the charge lattice is a bound state or multi-magnon minimum at a momentum that the liquid's boson does not carry. C163 excluded bound states and uniform tilts to all orders around U(1). The gap is harmonic amplitude of order one.
   - Try a variational or monotonicity argument for the lowest state in each momentum sector, for example a Perron–Frobenius or stoquasticity argument after the rotation, in the sectors where the bond is sign-free.
   - Or try a convexity bound on the band minimum as a function of harmonic amplitude, interpolating from C163's controlled region.
3. Ask whether time reversal (T) or the K_I symmetry acts on the candidate second mode in a way that forces it to mix with the liquid. Kramers-type or anomaly arguments (LSM with the extra antiunitary) count here.
4. **Targeted numerics, only as checks of a specific derived statement.** Ring ED or DMRG of the momentum-resolved lowest excitations at harmonic fraction 0.3 to 0.6. No broad searches, since C155, C159 and C163 have done those.

## Output

RESULT.md opens with a line saying which model did the work. It then gives a theorem-style summary with every claim tagged proved, derived, numerical, heuristic or assumed. The summary is one of three things:
- (A_W) derived in cases (a) and (b), with the range of validity;
- a strictly narrower residual assumption, stated precisely, saying where it bites and why the method stops there;
- an explicit violating bond, with evidence by the pre-fixed rules of C163 PROOF §4.2.

Derivations go in PROOF.md and the log in PROGRESS.md. Give a short referee-ready summary of the key step at the top of PROOF.md.

## Rules

- Astra is out of credits, so do the work yourself. Never call codex directly, and never use paid APIs.
- Keep numerics small and local where you can. Anything large goes to Typhon as one packed 48-core node job: `ssh -o LogLevel=ERROR typhon`, SLURM only, never on the login node, scratch /scratch/lustre/jkudlerflam/c167/. Retry ssh with waits. Record job ids in RESULT.md and do not block on them.
- Write only inside `codex/work_C167/`. Do not edit STATUS.md or REPORT.md, and do not commit or push.
- Never stash, reset or switch branches, do not pip install, and do not launch further agents.
