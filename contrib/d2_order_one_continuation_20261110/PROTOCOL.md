# Continued search: registration before new numerical output

OpenAI Codex, GPT-6 family, 2026-11-10. The user requested continued work
until a promising candidate is found. Baseline local commit: b8b6179e.
The previous contribution is immutable input. New records and checkpoint
copies live separately. Publication remains blocked by automatic approval
review; no push is attempted in this continuation.

The microscopic target and preference remain unchanged: one-site translation
invariance, physical qubits, nearest-neighbour interactions, a unitary CFT
with any c>1; prefer coefficient spread R<=3, label 3<R<=10 exploratory,
and keep wider bridges separate. Record invariant scales, all terms,
cancellations and actual tuning resolution.

## Root comparison: test the strongest existing finite-size lead

The unmodified physical bond is

H=sum[(XX-YY)/2+(XY-YX)/16+ZZ]-(7/4)sum Z.

It has R=28 and is not called comparable. The lower variational branches
read entropy slopes about2.7 at24/48/96, with strong density oscillations.
The96-site result lacks a converged bond-dimension comparison.

Before the new output, fix these tests:

1. Copy both96-site chi192 checkpoints and records to this contribution's
   `validation` directory and a new local checkpoint directory. Continue to
   chi256 using the unchanged C190 solver and prior wrapper, at least3 and
   at most6 sweeps, the existing1e-8 energy/1e-3 entropy step stop rule and
   local residual target1e-10. Compare each parity's entropy to its converged
   chi192 input and report energy variance. A step exceeding0.01 keeps the
   largest-size truncation test unresolved.
2. Run both192-site parity sectors atchi192, with intermediate chi16,32,64,
   128,192 stages rather than jumping from32 directly to192. Keep exact
   physical fields and the same random seed17. Use at least3 final sweeps
   and at most8. Preserve every sweep-boundary checkpoint. This is a size
   check at fixed final chi; it is not automatically a converged192-site
   answer. Increase chi only after inspecting convergence and a new entry.
3. Preserve both registered entropy windows, all parity branches, global
   variance, full density and entropy profiles, transverse/longitudinal
   correlations and the exact XX/twoXX finite-size estimator calibration.
   Analyze homogeneity and ordered-cat signatures before assigning modes.

An initial promising **numerical candidate**, not an established CFT, should
have reproducible excess above one beyond calibrated finite-size/chi errors,
two size pairs without a collapse toward one, no diagnosed cat/phase-jump
artifact, and a credible route to comparable couplings. The earlier1.2
screening flag is not an exclusion threshold: modest persistent excess is
retained for careful testing. A new preferred or exploratory point must be
tested at more than one size and both relevant parities before promotion.

Independent branches are registered in `bridge_scan/` (fine field and
moderate real-hopping continuation) and `attractive_scan/` (rational
attractive-junction route). `diagnostics/` reads previous checkpoints without
altering them. No projected or longer-range control is counted as a target.

## Independent periodic-ring comparison of the preferred attractive lead

After the existing C184 data identified phi=1.2, Delta=.8, V=-2 and before
new ring output: compare that `pref_ref` and its fixed rational neighbor
`pref_rat`=(Jx,Jy,D,Jz)=(14,-5,11,-12)/24 at L=14,16,20, both parity
sectors, at least four low Ritz vectors, with the previous translation-orbit
completion and residual checks. The exact field is absent. Preserve raw and
resolved entropies and gap semantics. Report size/window/sector drift without
using a small-ring gap or Casimir veto. Prior same-code XX/Ising/twoXX
controls at16/20 remain the estimator calibration. This is an independent
finite-ring check of a preferred-ratio lead, not a CFT certificate.

## Independent entropy-spectrum consistency check

Before inspecting new Renyi results, compute alpha=1,2,3 entropies from the
same saved Schmidt spectra for the two attractive points and R=28 bridge.
Use the identical central-three and fractional cut windows. Compare
`c_eff(alpha)=12*Delta S_alpha/[(1+1/alpha)*log(L2/L1)]` across sizes and
against the saved XX and two-XX controls. Alpha=1 reproduces the registered
von Neumann estimator. Retain all values, including disagreement; finite-size
Renyi corrections can be large, so agreement is supporting evidence and
disagreement triggers investigation rather than a universal rejection rule.
This is postprocessing of the same state, not an independent ground-state
method. It may help identify changing Schmidt weights or cats that distort
the von Neumann slope. Never read a checkpoint while its owner is updating
it: require completed record status and unchanged record/checkpoint hashes
after contraction, or discard that row. Do not modify original checkpoints.

## Independent periodic check of the new R=8.75 bridge

The bridge scan found `(Jx,Jy,Jz,D,hz)=(1/2,-1/2,5/8,1/10,-7/8)`
with both24→48 entropy-window slopes above1.68, small density modulation,
and decreasing transverse correlations. Before new ring output, register
the same full-Hilbert-space PBC L14,16,20 comparison used for the attractive
points: both parity sectors, k=4, tol=1e-11, translation-orbit completion,
and both raw/resolved data. Keep the size and momentum selection rules
unchanged. No parameter adjustment based on these rings is allowed without
a new adaptive registration. New JSONL/config/checkpoints are separate.

Before its output, also register an independent random-seed23 run of this
R=8.75 point at L48, even parity, chi192 (stages16,32,64,128,192; at most
8 clean sweeps per stage; unchanged stopping gates). Compare with the bridge
agent's seed17 result, including entropy and global variance. The purpose is
to expose metastable selection; it does not certify the true ground state.

After the completed24/48/96 results and before its output, add L64 for this
same R=8.75 point, both parity sectors, seed17, chi192 with stages16,32,64,
128,192 and at most8 clean sweeps. This intermediate size tests whether the
entropy trend is smooth rather than supported only by a doubling sequence.
No couplings, fit windows or branch-selection rules are changed. This is
one sequential worker in `independent_sizes`, with separate checkpoints.

## Focused stopping rule after the192-site candidate result

At17:39 UTC, both R=8.75 parity sectors have completed192-site chi192
optimization, with the even-branch entropy size/profile estimates agreeing
near1.45 and declining order diagnostics. The registered selected-sector
chi256 check remains mandatory before the candidate assessment closes.
No further field grid or weak single-size promotion will be started in this
continuation. If the wider R=28 exploratory192-site jobs are still unfinished
when that selected candidate check passes, defer them at their most recent
saved sweep checkpoint. Preserve all partial output and label it interrupted,
never completed or converged. This allocation decision is made before
their final192-site output and does not depend on its numerical value.
The task is to identify a promising comparable-coupling candidate, not to
exhaust every registered optional follow-up or claim a verified CFT.
