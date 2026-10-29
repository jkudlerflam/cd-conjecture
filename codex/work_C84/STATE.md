# C84 handoff

Completed 2026-10-29. Eleven-letter frozen member, n=11 active flavors,
p=12 data sites, cell length 13. Read RESULT.md then PROOF.md.

Delete E_1 from C80 and add K on L_1 followed by any marker. The last
logical bit is fixed to zero. The folded graph retains 25 vertices and
its full 625-dimensional product propagator, but only D_0,...,D_10
convert. The data endpoint D_11 is inactive. Its partner X_11 and loop
edge remain. Legal dimension (5p+3)2^(p-1), vacuum code 2^(p-1).

Certificate uses 22 simple roots, 253 unordered pair sums, delta2 greater
than 0.0164202883, delta4 greater than 9.11916e-7. All decisions use
integer/rational tests. The extra normalized perturbation budget is
2^-35. Loop exponent 66 failed with the chosen brackets; exponent 55
passes. This rejected parameter is retained.

Final schedule is D=1000, Delta=4-2cos(pi/14), phi_A=pi/24, phi_B=0,
psi_A=pi, psi_B=pi/2, beta=1/d_13 with d_(j+1)=1004d_j-d_(j-1),
t_d=beta*2^-55, g=2^-216, K=2^2048, lambda=2^-480,
physical field=lambda^2*alpha^4*beta^2*g^2*2^-101.
The same length-13 midpoint selector as C80 applies. MEMBER.json defines
all couplings. Couplings are fixed for the whole Hamiltonian, not retuned
on cells of other lengths.

Finite-K zero-port splitting is included in physical units. Finite-g
error <3.380e-21, finite-K error <7.946e-59, local two-cell port error
<2.518e-29. Their sum is below one eighth of 2^-35. The uniform zero-port
selector follows by compression of C80 plus a positive boundary pin.

Independent legal matrices at n=1,2 agree with the C75 bit-carrying
parent with last bit fixed and last conversion removed. Full invariant
ED dimensions 2221,42424. All 4,16 levels exhibit fourth-order port
errors, residuals <1.6e-13. Exact grammar checks through p=7. All ten
charge-preserving single-letter merges audited. Exact additive-charge
space has dimension three, spanned by identity, marker and payload Q.

Lower bounds are scoped. Eleven is minimal for the faithful shifted-bit
two-rail product encoding. Charge and retained port/conversion structure
alone give seven; a dense frozen marker family alone gives only three.
No general ten-letter no-go was obtained. The fold was retained.

Infrared c=11>10 is conditional on nonzero-port dressed-defect comparison,
volume-uniform matching and multichannel stability. No measured central
charge, no full counterexample, no independent referee. Scripts use only
local copies of C80 sources, with changes documented in PROOF. All work
stays in work_C84; project STATE outside this folder was not modified.

The root certificate's finite_g_power=-199 is the inherited scalar
routine's threshold, not the final two-species schedule. Use -216 from
finite_budget.json and MEMBER.json. This distinction is inherited from
C80; all final acceptance bounds use the two-species prefactor.
