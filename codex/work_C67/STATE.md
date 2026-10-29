# C67 local state

2026-10-28. Task completed as a bounded analytic and numerical investigation, with the global claim still unresolved. No global proof or supported c > 1 point was found. Read RESULT.md for the status-tagged derivations and PROGRESS.md for the work record.

The canonical family is (Jx,Jy,u,D,b) multiplying XX,YY,ZZ,XY-YX,Z. The uniform-unitary C2 row has Jx != Jy and (D,b) != (0,0). Its strict internal-Z2 subset has D b != 0. It has no nonzero additive spin charge, since nonzero hopping requires neighboring charge coefficients to agree while nonzero pairing requires them to be opposite.

Controlled material includes the supplied free, XXZ and XYZ descriptions; the explicit gap cone u <= -sqrt((Jx+Jy)^2/4+D^2)-|Jx-Jy|/2 with b!=0; and the factorized family Jx=-A,Jy=-B,u=C,b=-2sqrt((A+C)(B+C)), A>B>0,C>-B. Its two homogeneous product states remain exact eigenstates for every D because the DM defects telescope. Their ground-state status at small D uses stated parent-gap and ordered-phase stability inputs. It is not asserted at arbitrary D.

All new ED used every momentum and both parities, four levels per block, L<=20. There are 70 point-size records and a maximum block dimension 26272. Maximum eigenvector residual 9.3e-10. Independent dense-spectrum/intertwining checks passed at L=6; Ising and XX energies match exact answers through L=16. See analytic_checks.json.

The four grid candidates are frozen in grid_candidates.json. Z2_21_295 has c_chord(20)=2.259132 but c_size=0.960061 for L=12,16,20 and 0.867243 for L=10,14,18. Z2_1_358 has incompatible mod-4 sequences; the inherited n=11 FES run through chi=24 gives 1.020698 and is only preliminary. The n=1 FES runs for these points fail convergence and contain unphysical energy outputs. The two patch candidates also fail entropy consistency. No grid point is assigned a thermodynamic c.

Junction family is (-1,-0.5,0.25,D,-g). The mass tuning sets the m=0 even/odd gap ratio to eight, so that ratio is not an independent Ising test. At D=0.20, g_L decreases to 1.99101535 at L=20; c_chord=0.502816 and the slow odd descendant difference is 0.131622. At D=0.24 that difference is -0.025741. Velocity extrapolations suggest a collapse near D=0.20 to 0.21, g near 1.989, with appreciable fit dependence. No interacting exponent or extra sector was established. At fixed g=1.9891 and D=0.30 a ground-momentum change produces a misleading chord fit near 1.19; the L=16 to 20 entropy step is 0.229217.

Scripts and output names are listed in RESULT.md. ed_core.py is copied from the corrected C6 kernel. Scripts write only locally. analyze.py reads selected inherited FES data as well as saved local runs, so summary.json's inherited-grid section could change if the upstream run is extended. Its current values are frozen in RESULT.md. No network, cluster, subagent, or tensor-network calculation was performed. Project-level STATE.md was not changed because it is outside the writable scope.
