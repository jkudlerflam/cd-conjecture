# Current d=2 numerical status: ground-sector arithmetic check

OpenAI Codex, GPT-6 family, 2026-11-10. Numerical data audit only; no new
simulation or human review. Sources were pulled at commit 56bbf259.

The C190g reruns supersede this chat's earlier unconverged L=128 downturn.
They support a c approximately two candidate through L=256, but every
L=256 candidate still fails the registered bond-dimension entropy test C4.
The committed C195 L=512 snapshots are running, not finished results.

Checking the raw completed files identifies a small correction to the
current summary: sym40_n316 has **odd** lower-energy states at both L=128
and 256. Their final-sweep Sbar6 values are 1.8310312793 and 2.0558758053,
giving c_two=1.946292. The value 2.074 reported in the summary instead
matches the even-sector comparison, which is not the ground branch.
The energy splittings selecting odd are about 1.80e-5 and 6.31e-6, well
above the observed last-sweep energy motion of order 1e-9.

Using the same C1–C4 selection hierarchy as the production analysis gives:

| Key | Ground parities, 128→256 | c_two at W=6 |
|---|---|---:|
| p153c | even→odd | 1.935896 |
| sym40_n316 | odd→odd | 1.946292 |
| rat_n316 | odd→odd | 2.032657 |
| sym40_n18 | even→odd | 1.978055 |
| p153_n316 | odd→even | 1.995620 |
| sym80_n316 | odd→even | 1.942079 |

Thus the current ground-branch range is **1.94–2.03**, rather than 1.94–2.07.
This does not alter the scientific verdict: a stronger candidate, with
remaining convergence and infrared checks. The canonicalized final
entropies differ slightly from the final-sweep values; both estimates
are retained in results.json and the correction holds for either.

Reproduce with `python3 contrib/d2_status_audit_20261110/check.py`.
The script reads the original data without modifying it, records hashes
of selected files, preserves the failed C4 status, and writes results.json.
