# C88 progress

[Read] Project memory, simplest-background specification, C84, C85, C83 and the C84/C85 referee report H1 were read. The pinless automaton and the distinction between invariant-sector restriction and physical ground-state selection are the starting points.

[Completed] The literal d=13 bond is in bond13.py. Local copies of the ten C85 certificate sources are unchanged. All eight scripts in run_all.py pass, with exact certificates regenerated locally and compared to C85. No parent files are modified.

[Proved] Unique coloring and the mismatch inequality hold. The complete leading grammar correction has color coefficient one on R_0 and E only. The mixed hopping/material term vanishes. F3 now excludes frozen pin defects as well as wrong lengths. F4 uses the physical penalty JK. Every remaining letter has a creation witness inside the retained legal component, so none is removable by the same superselection argument.

[Numerical] Nonzero-bit, nonzero-port tests with two vacancies show K^-2 suppression, including a one-cell closure of dimension 1045 and a separate two-cell closure of dimension 3078. The uncolored controls retain illegal material weight. Exact legal correction enumeration reaches 19267 states. The regenerated material certificate has 26 roots and 351 separated pair sums, with C85's exact certificate outputs unchanged.

[Open, explicitly assumed] Ground-state selection of the pin-free mobile sector needs an energetic comparison. Conservation alone is insufficient. PROOF.md states S_pin as part of the inherited mobile dressed-defect assumption. Frozen selection is proved, and the mobile reference pin cost is positive. The other inherited matching, stability and conformal-identification inputs remain open. No full mobile ground-state theorem or measured central charge is reported.

[Delivered] RESULT.md, PROOF.md, MEMBER.json, STATE.md and reproducible local Python checks. The conditional candidate is d=13, c=14>12. No network commands, subagents, cluster jobs or writes outside this folder were used.
