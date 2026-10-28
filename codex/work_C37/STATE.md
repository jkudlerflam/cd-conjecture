# C37 local state

[Proved, completed finite task] RESULT.md and PROOFS.md certify the original frozen 23-letter candidate for every n=1,...,23. Parameters are D=1000, shifted stub mass 4, total neutral express phase pi/2, neutral loop phase pi, t_d=10^(-2000), normalized uniform field 10^(-1000). No alphabet change. Exact modular pair-product polynomials are squarefree modulo 65521 with i=24297. This excludes every nontrivial number-conserving quartic resonance. C33 positivity certificates and its certifier supply crossings at all 23 sizes, with an exact uniform normalized lower bound 10^(-250).

[Proved] Explicit denominator/discriminant momentum margins, velocity bounds, an analytic root-perturbation tolerance and a finite-g schedule are in PROOFS Sections 5-6. They are extremely conservative. Physical scales include alpha_n^4 g_n^2 lambda^2. At n=23 log10(delta_bound) is about -1.06948e15. Do not confuse this lower bound with a measured margin.

[Assumed/open] Local thermodynamic elimination and interaction-induced dressing remain the prior infrared inputs. The interacting schedule uses separate local-elimination and dressing constants. Reduce g for the retained interaction self-energy before choosing lambda; decreasing lambda alone does not suppress that normalized dressing. All n beyond 23 remains open. These certificates do not apply without recalculation to C32's different mobile core.

[Proved file map] Main exact arithmetic is modular.py, certify_pairs.py and pairs_n*.json. Positivity is inherited C33 code in positivity_C33.py, reused by complete_positivity.py and verified in verify_bounds.py. Independent exact checks are audit.py/audit_results.json. The successful literal numerical check is literal_coefficient.py/literal_coefficient_results.json. The earlier subtraction check in literal_results.json lost precision and is not certificate evidence. scales.py/scales.json give logarithmic displays. PROGRESS.md is the session record.

[Proved scope] All writes stayed in work_C37. The project-level STATE.md was not modified because the user forbade writes outside this folder.
