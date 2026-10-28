# C26 local state

2026-10-28. Complete requested G5/G6 calculation. Read RESULT.md for every
derivation and PROGRESS.md for the handoff. No files outside work_C26 changed.

[proved] Sidon parameters preserve C19's exact code, positive weights,
source poles, parity and second-order hopping plus exchange-generated fields.
Bare logical parity conjugates rho to minus rho, while physical parity
induces the same code action. Thus each finite-order traceless canonical
effective coefficient is analytic in rho squared and divisible by rho squared.

[numerical, n=1,2] Two-cell K4 traceless norm/rho^2 = 1.2376052983,
4.3765621160. Connected adjacent-bond K4 values = 0.8622186591,
1.6116197398. Three-cell invariant dimensions = 55296,1048576.
The finite-coupling graph checks give fifth-order remainders. Full/reduced
operators and an independent three-cell n=1 calculation agree. Exact
rational Jacobi tests passed at n=1,2,4,8,12.

[assumed/heuristic] Local massive-sector elimination plus dressed parity-only
Sidon Gaussian stability gives a payload correction relative to
W=lambda^2 rho^2/Delta_src starting at lambda/Delta_src. No uniform-in-n
coefficient bound, compiler, full mobile transfer or measured central charge.
Payload-scalar terms can interact with cell lengths and must be controlled
separately by G3. C25 and C27 completed their separate calculations and
were read at final handoff.

Files: sidon_cell.py, checks_sidon.py, three_cell.py,
independent_checks.py; checks_results.json, three_cell_results.json,
independent_results.json; coefficients_n1.npz and coefficients_n2.npz;
three_cell_coefficients_n1.npz and three_cell_coefficients_n2.npz.
All generated coefficients act on the code and use the canonical polar gauge.
