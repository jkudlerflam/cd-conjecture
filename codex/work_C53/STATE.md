# C53 local state

Completed 2026-10-28. Read RESULT.md for the result and PROOF.md for derivations, theorem statement and hypotheses.

Two modifications must not be conflated.

1. H_tune retains the original port coupling and adds a four-site diagonal counterterm q(g)=9g/64-243/11648 plus field. Its second-order hopping is still -9/32, anisotropy g, with a proved O(lambda^3) spatial remainder. It is not covered by the published instantaneous-interaction theorem as a microscopic model.
2. H_exact deletes the PP port bond and inserts direct logical exchange on each BPPB block. It retains d=10 and original cell checks, but physical support is four sites. In frozen PBP, H=A+s H_BM acting on bare c qubits, s=9lambda^2/16. Exact code invariance, noncode ground-relative threshold 2/5 and pressure error <=T log(1+15 exp[-2/(5T)]) are proved. This is an architectural change, not a correction of order lambda^3.

BM 0907.2837 Theorem 1.1 and 1012.0403 Theorem 1.1 are verified from primary sources and apply to H_exact's exact logical sector. Use v(±1)=1/2, small g and h(g)=1/3-g-nu(g,k*) with k*=arccos(-1/3). Correlations and responses follow; an interacting thermal/stress theorem is still missing. Unknown smallness radii do not certify g=.02. For g=0, h=1/3, the full frozen sector has proved thermal coefficient with c=1 and the free Dirac continuum interpretation.

Unrestricted type selection, original nearest-neighbor scope, mobile-cell matching and c=2 for the wrapper remain open. The collaborator's Jp tuning remains useful inside the nearest-neighbor class; it does not eliminate the retarded-action issue.

checks.py runs small local ED and logical particle-number ED, and writes checks.json. The first ARPACK run failed from a uniform Krylov starting vector; deterministic generic starts and explicit residual checks fixed it. Final run passes. Largest full matrix is 1024, largest logical number block is 924. No cluster jobs or subagents. No writes outside work_C53. Parent STATE not editable under this task's scope. Primary-source browsing followed the higher-priority session rule despite the background's offline instruction.
