C189 infinite-chain test of the Fibonacci-cusp candidate (Claude Opus 5.5 agent, 2026-11-09).

What it is. C173's effective model (PROOF eq. 3.2) at bare couplings on the symmetric line, phi = pi/2, |w| = 1,
on the INFINITE chain: blockade-exact, U(1) infinite DMRG (idmrg.py, McCulloch growth with a two-cell window,
C176's label scheme), at fixed density n = 3/16 (16-site cell) and n = 5/32 (32-site cell). For each bond
dimension chi it records the half-chain entropy S and the correlation lengths of the transfer matrix resolved by the
charge dQ = 0, 1, 2 carried between bra and ket. A c = 2 liquid gives S = (c/6) ln xi with c near 2 and all
correlation lengths growing together; a pseudospin gap makes xi_1 saturate while xi_2 keeps growing.
Typhon jobs for n = 3/16, 1/8, 1/6 and the controls were queued but had not started; this package runs the main
density to chi = 768 and a second density, n = 5/32 (bulk density closest to C175's 0.17), to chi = 512.

Validation already done (codex/work_C189/RESULT.md): the first window reproduces C176's fibdmrg to 1e-11; the
controls read c_fit = 1.015 (hard rods, c = 1) and 1.989 (free fermions with the candidate's band, c = 2) up to
chi = 96, and the n = 1/8 candidate read c_fit = 2.008 up to chi = 96.

Registered rule (RESULT.md Part 2 and amendment A1, written before any of these runs):
  S_bar = mean of the last two iterations' S. c_fit = 6 x slope of S_bar against ln xi over the largest four chi,
  xi = largest correlation length over the sectors.
  c = 2 liquid: c_fit >= 1.6 and xi_1/xi_2 changes by less than a factor 1.5 over the last chi doubling.
  Pseudospin gap: xi_1 grows by less than 10 percent over the last chi doubling while xi_2 or xi_0 grows by more
  than 30 percent, or c_fit < 1.4 with the local slope falling toward 1/6. Anything else is inconclusive.
  Limit: a pseudospin correlation length beyond the largest xi reached is invisible here.

Needs python3 with numpy and scipy (nothing else). Run: bash run_desktop.sh   then   python3 summarize_inf.py out
Each finished chi appends one line to out/<tag>.jsonl; out/<tag>.iter.jsonl logs every iteration; out/ckpt/ holds
resumable checkpoints (a few hundred MB at chi 768; do not sync them).
This is the effective model, not the full qubit chain.
