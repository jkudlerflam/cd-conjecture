# C88 result

[Proved finite ingredients; conditional infrared conclusion] The explicit bond in [bond13.py](bond13.py) has eleven material letters and V_0,V_1, so d=13. On its invariant pin-free material sector, C85's n=13 certificate applies unchanged. The conditional count is

    c=13+1=14>12=d-1.

[Proved] The removed pin counts squeezed L_1-marker junctions on the physical ring. Every literal move conserves this count, including illegal transitions, colored hopping and the port. In its zero sector the full finite-K material Hamiltonian equals C85's. This is an equality of material operators, not of the finite-K two-color and three-color mobile bonds.

[Proved frozen selection; assumed mobile selection] The frozen ground state on compatible rings lies in the pin-free, length-15 sector. The extended F3 bound includes pinned length-15 cells, whose individual excess is at least 2cos(pi/16)-2^32/(K-2^17)>1 at the certified schedule. On the mobile ring, ground-state selection remains conditional on the common-bulk-subtracted dressed-defect comparison, explicitly extended to pin defects. PROOF Section 2 states this assumption as S_pin. Conservation alone does not select a ground state. The projected half-density reference pin cost is positive and of order J.

[Proved] Legal coloring is unique at nonzero occupancy, N_matbad<=N_bad, and every projected swap has amplitude -t. There is no extra color mode. The complete leading grammar correction is

    -(t^2 C+J^2 A+2J^2 lambda^2 Z)/(JK).

C is the color-ambiguity contact counter, with coefficient one only on R_0 and E followed by a vacancy. A is the vacuum-tail counter and Z the wrong-end port counter. The mixed hopping/material term vanishes. The first term produces a t^2/(JK) bit field and its coupling to clock and vacancy fluctuations, without a port factor. Higher terms are retained in the Feshbach operator and F4 remainder.

[Numerical] Nonzero-bit, nonzero-port tests with two vacancies give K^-2 suppression of illegal material weight. In the 1045-dimensional closure, the weights at K=100,1000,10000 are 1.33164e-7,1.52294e-9,1.54444e-11; the uncolored control approaches 2.64428e-5. A separate two-cell, two-vacancy port test has dimension 3078. The complete leading correction passes every legal ring of lengths 3,4,5, up to 19267 states, with discrepancy below 4.5e-16. Frozen F3 and local F4 have separate full-matrix checks. All eight scripts in [verification.json](verification.json) pass.

[Proved, regenerated integer/rational certificate] The unchanged C85 material pipeline gives 26 simple Fermi points, 351 separated unordered pair sums, delta_2>0.0156742281 and delta_4>3.2513114e-6. The finite-g, physical finite-K and isolated two-cell port errors sum to less than one eighth of the 2^-35 normalized allowance. The exact root, budget and selector certificate files agree with C85. The carrier still has 29 vertices, including its inactive final data endpoint, and its product resolvent has dimension 841.

[Proved, restricted reduction test] No further letter can be removed by the same superselection reasoning. Within one legal pin-free two-cell component, every remaining letter can be created from a word that lacks it. [sector_checks.json](sector_checks.json) gives explicit witnesses for all thirteen letters. This excludes deletion from the retained dynamics, not a redesigned d=12 mechanism.

[Assumed] A microscopic c=14 ground-state result still needs massive-sector and finite-J vacancy matching, S_pin and the full mobile dressed-defect comparison, interacting multichannel stability, and conformal-observable control. The mobile K bounds are implicit and include the t^2/(JK) contact field. No central charge was measured and no unconditional counterexample is claimed.

The complete construction, derivations and scope are in [PROOF.md](PROOF.md). [MEMBER.json](MEMBER.json) specifies the conditional family. Run `python3 -B run_all.py` from this folder to reproduce the checks. No files outside work_C88 were modified.
