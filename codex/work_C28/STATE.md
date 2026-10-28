# C28 local state

[proved] The work is complete as a partial research result, documented in RESULT.md. A positive initialized 17-letter charged bus has an exact 2^n code at every n and exact U(1). Its vacuum walk is a Laplacian, with uniform history, so it meets the positive-core premise of the supplied vacancy wrapper. This is a change from C22's adjacency-walker code.

[proved] At n=2, endpoint densities give independent stored-data fields. A common field and flux pi/2+.005 produce four simple leading Fermi crossings, certified by exact rational inversion and elementary trigonometric bounds in certificate.py. The general repeated diagonal two-site readout has only four field profiles, so this readout method cannot specify arbitrary independent fields for n>=5. Scalar path phases do not repair C21; unequal charge conversions do not repair C22's leading phase structure.

[numerical and heuristic] The n=2 generated quartic interaction puts Cooper transfer on the irrelevant side of the leading one-loop BKT flow. The matching coordinates are x/g^2=-.09435238259 and y/g^2=.08124029522. All other distinct-pair momentum mismatches are nonzero at the reference. No central charge has been measured, and no finite-lambda phase certification is made.

[numerical] Complete two/three-cell legal-sector ED at n=1,2 and lambda=.08,.04 verifies the exact full-in-g second-order bond plus fields. The leading remainder is lambda^4 by an exact open-chain symmetry. Largest dimension 110592. See the ED JSON records and RESULT for errors and residuals.

[proved, file map] ff_bus.py is the new local positive core; bus_base.py copies the inherited C22 source calculation. explore.py supplies coefficient and scattering routines. verify.py, checks.py, certificate.py and handoff.py reproduce the final records. logical_chain_n2.npz is the ready-to-use 16 by 16 bond handoff, with leading, order-g^4 and full-in-g alternatives. PROGRESS.md records failed attempts and the corrected literal n=1 eigensolver issue.

[assumed/open] No all-n target, all-n RG sign theorem or universal no-go was obtained. Conditional alphabet counts are 68 frozen and 1680 mobile, not completed unbounded-c d0 values. The next work must address the all-n spectrum/field problem or change the encoding. The external project STATE was not modified, as required by the task's folder restriction.
