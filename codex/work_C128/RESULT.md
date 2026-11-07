# C128 result

2026-11-06. Derivations are in [PROOF.md](PROOF.md), and reproducible tests and their limitations are in [CHECKS.md](CHECKS.md).

[Open] The interacting Fock domination remains unproved. No counterexample to finite C(d), or to every possible scalar Fock majorant for stoquastic positive-transfer chains, was found.

[Proved] The decorated-history construction is explicit for XX, Heisenberg, and XXZ at Delta>=0. The nonnegative checkerboard contraction is naturally a two-column matrix B_m=K_0 K_1, with

    Z_m(L,beta)=Tr B_m^(L/2).

A related positive operator U_m=K_0^(1/2) K_1 K_0^(1/2) gives T_m=sqrt(U_m) and Z_m=Tr T_m^L on even rings. The square roots need not preserve entrywise nonnegativity. A one-site positive operator must not be silently identified with a nonnegative history kernel.

[Proved] For every requested XXZ model, the number N_i of events incident to one column obeys the normalized generating-function bound

    sum_config weight(config) z^(N_i)/Z <= exp[beta(z-1)], z>=1.

This controls the event count at order beta M. It does not control temporal frequencies or give independent column clocks. Paired neighboring jumps occur at identical times. Their small-step weight is order t, while independent column events have weight order t^2. The continuous joint measure is singular relative to independent diffuse column clocks. Conditional bridge intensities also diverge as the reciprocal remaining time. These failures occur already in XX.

[Proved new XX bound] Using C125's positive free product and convex midpoint quadrature, for exchange velocity v=2M and every even L>=2,

    log Ztilde_L <= log 3+beta v[csc(pi/L)-L/pi],
    R_L <= 2 log 3+2 beta v[csc(pi/L)-L/pi].

Consequently

    R_L <= 2 log 3+(16sqrt(2)-64/pi) M beta/L, even L>=4.

The coefficient is 2.2555842822, smaller than C119's 10pi/3. Including L=2 gives coefficient 8-16/pi=2.9070418211. These also bound L(E_(2L)-2E_L)/M on the respective subsequences. The result adds across independent XX stacks. It does not improve the constants for C119's other free families without another argument.

[Proved] A literal linear-frequency Fock product majorizes the XX trace with

    S=2, V'=2pi v=4pi M,
    K=log 3-2 log(1-exp(-1))=2.0159625794... .

This is an aggregate trace bound, not a same-frequency spectral comparison or a positive-kernel comparison. It follows from the stronger bound above. The denominator is controlled without an extensive loss by Z_(2L)>=lambda_0^(2L).

[Proved lattice failures] Linear decay at every Matsubara frequency fails for the locally detected XX branch inherited from C125, whose gaps grow logarithmically at high frequency. The unchanged XX normalized-trace majorant also fails for Heisenberg. At L=4 its exact zero-temperature slope is 3-4 log 2=0.2274112778, above the proposed XX slope sqrt(2)-4/pi=0.1409740176 at the same M=1/2. Neither failure disproves a sufficiently enlarged scalar majorant.

[Assumed sufficient missing input] An interacting multi-event estimate must control the full normalized crossed partition function by finitely many bosonic species with additive decay costs, plus an explicit zero-frequency factor. PROOF equation (18) states one such estimate using the valid lattice rate alpha(omega)=omega/(4eM) at low frequency and log(omega/(4M)) at high frequency. Its logarithmic ultraviolet tail is summable and yields an explicit beta M/L bound, proved in equation (19). Thus high-frequency linear dispersion is unnecessary. The missing step is the many-state comparison, including zero-frequency sectors, which neither count control nor reflection positivity establishes.

[Numerical] There are 192 finite-Trotter contraction and positivity checks on XX, Heisenberg, Delta=0.5 and Delta=2, with m=1 through 4, L=2,4,6, and beta=0.2,2,6,8. Maximum relative contraction error is 5.58e-12. Thirty-five exact XX trace tests extend through L=1024 and beta=32768, with explicit product-tail bounds. The new inequalities pass after accounting for those tails. The generous linear Fock candidate also passes the sampled interacting finite-Trotter matrices, which is not a continuum or all-temperature proof. The unchanged XX bound fails in 24 interacting finite-Trotter samples.

[Work record] All files and Python execution stayed in work_C128. The local [STATE.md](STATE.md) serves as the handoff because outside writes were forbidden. No cluster connection, network shell command, parallel agent, or independent referee review was used. [run_typhon.slurm](run_typhon.slurm) is prepared for an extended transfer scan and was not submitted.
