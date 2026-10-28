# C17 local state

Completed 2026-10-27. This is the task-local handoff because the task forbids
editing the project STATE.md outside work_C17.

The task was to audit physical mechanisms for c>1 in a one-site-uniform
nearest-neighbor qubit chain. RESULT.md is the complete tagged report with
derivations, limitations, tests and ranked open mechanisms. PROGRESS.md records
the inputs and completion state. checks.py reproduces the small numpy/scipy
checks; checks_output.txt records all passing results.

Closed mechanisms include a literal conserved one-site marker with a walker,
literal physical even/odd subchains, exact U(1) using the supplied XXZ phase
classification, exact SU(2) using the Heisenberg reduction, the existing free
and XYZ classes, and completed Abelian parton cosets under their stated IR
identification. K_FF(2)=4 closes only longest-zero-segment selection.

Still open are isolated neutral-sector multicritical points, emergent
higher-level currents without the required microscopic current embedding,
strongly interacting incommensurate multiple fluids, quantum-selected cells
with dressed payloads, unfamiliar constrained neutral sectors, and more
remote genuine noncompact possibilities. They overlap and are ranked in
RESULT.md section 9. No proof of c<=1 and no c>1 counterexample was obtained.

The fresh numerical evidence is calibration and finite algebra, not a search.
PBC ED at L=6,8,10,12,14 gives XX c_S(L=14)=1.027712 and c_Cas=0.999204;
Ising gives 0.502011 and 0.499976. Both controls have v=2. All checks passed.
The report cites C15/C16 progress for their complementary results and does
not claim to have independently verified those new calculations.

No work remains in progress. No files outside this folder were modified.
