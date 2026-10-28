# C57 handoff

2026-10-28. Partial result, no completed interacting thermal theorem.

[proved partial] PROOF.md Section 3 proves a four-history bootstrap from
uniform one-step mixed thermal/derivative bounds. It retains the reference
thermal drift through a Hessian term, and includes terminal scalars and
weighted history recurrences. The actual uniform GN norm and normalization
estimates for XXZ are not proved.

[proved coefficient] PROOF Section 6 computes the first interaction Taylor
coefficient of the XXZ free energy at fixed field and fixed density. At fixed
density it is T^2/(12t^2)+O(T^4), consistent with the first-order dressed
speed and thermal c=1. No fixed-nonzero-interaction asymptotic follows.

[read, important prior-work lead] Rutgers 127th Statistical Mechanics
Conference, B12, 17 December 2024, page 10 of
https://www.math.rutgers.edu/images/127_SMC_short_talks.pdf . Bruno Renzi,
Alessandro Giuliani, Vieri Mastropietro and Venanzio Alessio Santini announce
the leading thermal correction for weakly interacting quantum chains.
No manuscript or explicit XXZ hypothesis was located. Do not claim absence
of prior work. No author was contacted.

[read] SOURCE_AUDIT.md records primary equations. Closest scalar estimate
is Giuliani–Mastropietro 1211.3030 (3.90), with its continuation after
(3.91), plus (3.69), (3.84), (3.87). It proves the desired shape for Ising,
not for a marginally interacting Luttinger reference. BM 0907.2837
(81)–(82) and MP 2108.13708 (8.61) are beta-function estimates, not P.
BFM 1303.3681 Section 3.1 retains the free-energy scalar but only proves
the displayed limiting convergence. The physical free-energy theorem in
2305.06679, Theorem 1.4, explicitly assumes Conjectures 1.2 and 1.3.

[read] C55/C56 progress still leaves the retarded all-scale comparison open.
C58 is independently treating ensemble transfer. The collaborator thermal
notes already prove the single-scale winding and Matsubara localization
lemmas. Cite them rather than duplicating them. Their all-scale propagation
is not automatic.

[numerical] checks.py and checks.json contain passing small controls.
All changes are confined to this folder. The project STATE is outside
the permitted write scope. No subagents or shell network commands used.

[next work] Obtain the announced manuscript, or prove PROOF (5), (13) and
the terminal/reference-normalization requirements with one common finite-T
renormalized expansion. This is the remaining mathematical task, not a
routine consequence of bounded marginal couplings.
