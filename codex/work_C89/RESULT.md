# C89 result

[Result] C88's S_pin is not discharged. The requested full finite-K mobile ground-state selection still needs a uniform low-frequency comparison of operator-valued defect energies after subtraction of the common bulk. The supplied vacancy normal form controls the ultraviolet interaction remainder, and does not establish that comparison. C88 also specifies no numerical J/t or complete numerical mobile K threshold.

[Proved] In the zero-port, K=infinity vacuum-clock model, removing any finite set of pins strictly lowers the energy. The variational bound is

    E_pinned-E_unpinned >= J sum_pins <A_terminal> > 0.

This holds at every finite size with J,t>0. It does not supply a thermodynamic lower bound on the contacts or include auxiliary and finite-K effects.

[Derived reference costs] At half occupancy the leading pin cost is 0.583287978365481 J. A length-14 or length-16 cell costs 0.000133769420252 J relative to the same material-length bulk. A compatible (14,16) pair costs twice that. Long-cell excess grows as 0.004825498950096 J a l, where a=0.29735763271532445.

[Proved conditional comparison] With W equal to the total material length of wrong cells plus the number of pinned length-15 cells, a sufficient bound is

    E_sigma-E_good >= J W
      [a(gamma-5sqrt(2)|lambda|)-epsilon_charge(J/t)-eta_K],
    gamma>2^-17.

Here epsilon_charge must come from the missing common-bulk comparison. The finite-K local budget contains 2/[(J/t)^2 K], not just 1/K. The lambda-squared port energy has scale J a lambda^2; its payload dispersion has scale J a lambda^2 E_unit, with log2(lambda^2 E_unit)=-2074.577410080.... All these terms are retained in [PROOF.md](PROOF.md). Their scale separation does not determine the missing charge coefficient.

[Proved, legal projected model] A one-marker ring has energy at least E_gap-4J. The selected-sector trial has energy E_gap+J a_L N e_15. At half occupancy every compatible N>=2865 therefore excludes the whole one-marker sector. This does not exclude an arbitrary long-cell island inside selected bulk.

[Proved, legal projected model] The marker-free sector has energy E_gap and is excluded on every compatible ring with positive mean contact by the same selected trial.

[Numerical] All 24 quantum-clock and 24 scalar-contact ED comparisons have positive costs. Maximum dimensions are 6750 and 6435, with eigenpair residuals below 8.89e-12. The quantum calculations retain distributed marker energies and contact-dependent clock hopping. They include length-15 pin and wrong-length defects at high occupancy, and b=3 defects at half occupancy. Separate half-density scalar tests use the b=15 reference costs. No half-density b=15 quantum calculation or scheduled-band simulation is claimed.

[Numerical] The independent operator-valued second-order calculation gives

    Delta E = J a_L Delta E_material
               -(J^2/t)(B_defect-B_good)+O_size(J^3/t^2).

The computed correction raises the pin cost in the tested quantum systems. It lowers the length-30 replacement cost, which remains positive. The sign cannot be assumed uniformly.

[Reproduction] Run python3 -B check_defects.py, then python3 -B summarize_checks.py and python3 -B second_order.py from this folder. Parameters and residuals are in ed_results.json, checks_summary.json, and second_order_results.json. All writes and Python execution stayed inside work_C89. No network commands, external edits, subagents, or cluster jobs were used.
