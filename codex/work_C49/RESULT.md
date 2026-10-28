# C49 result

[proved, verdict] This route gives an explicit quasi-local prethermal dressing and a sharper finite-volume spectral comparison. It does not prove Input S or c=n+1 for the microscopic chain at fixed nonzero couplings. The remaining estimate is stated precisely in PROOF.md Section 7. No counterexample to the proposed phase is found.

## Constructed dressing

[assumed, imported chain] C45 supplies the rooted-cell identity U0 H U0*=A+B+R, with A_j>=Delta Q_j, B locally block diagonal, and ||R||_a<=epsilon=C_SW lambda^3. The gap is the complete wrapper gap. Let d bound ||B||_a and put Lambda=d+2epsilon.

[proved conditional on that identity] There is an explicitly constructed finite product of local Schrieffer-Wolff transformations with

    m=floor[a Delta/(128e Lambda)],
    U H U*=A+D_m+O_m,
    [D_(m,X),P_X]=0,
    ||O_m||_(a/2)<=epsilon_m=2^(-m)epsilon,
    ||D_m||_(a/2)+||O_m||_(a/2)<=Lambda.

The retained Hamiltonian differs from C45's full second-order logical operator by at most 2epsilon in the weighted interaction norm. At fixed n and wrapper parameters, Lambda=O(lambda), so the leakage bound is exponentially small in 1/lambda. The proof spends locality weight at each step and never divides by a logical level spacing. It is a finite iteration, not a convergent all-order diagonalization.

[proved, locality and symmetry] The additional generators have total weighted norm at most 4epsilon/Delta and obey a volume-independent Lieb-Robinson bound. They preserve the exact charges and cyclic gluing. For the mobile representation, each generator on r cells has P_X on one side and conserves total length, so it vanishes unless those cells occupy between br and (b+1)r physical sites. This proves a physical support bound, rather than assuming that ordered-cell locality implies physical locality for arbitrarily long cells. The physical lift uses C44's fixed-outer-marker representation, as stated in PROOF Section 1.

## Finite-size spectral theorem

[proved] Let K=P(A+D_m)P, E=min spec K, and delta=Delta-2Lambda>0. If

    M epsilon_m^2 <= delta^2/8,

then every ordered eigenvalue k_i of K with k_i-E<=delta/8 has a corresponding ordered microscopic eigenvalue h_i satisfying

    k_i-8M epsilon_m^2/delta <= h_i <= k_i,
    |(h_i-h_0)-(k_i-E)| <= 8M epsilon_m^2/delta.

The eigenvalue count follows from the Schur complement, including multiplicities. The essential improvement is

    ||Q O_m P|| <= sqrt(M)epsilon_m,

because leakage events on disjoint supports are orthogonal on the code. Absorbing the complete local diagonal correction is essential.

[proved, scale] Relative to a reference spacing v/M, the resulting error budget is small when

    M << sqrt(delta v)/epsilon_m.

With v=lambda^2 c_n v_*, this certifies an exponentially growing but finite size window as lambda decreases. It does not prove that the interacting K has a Gaussian spectrum. At fixed nonzero lambda, these estimates do not give an o(1/M) error as M tends to infinity.

[proved, physical-ensemble boundary] The theorem is unconditional on the rooted tensor product once the imported local reduction is supplied. A conserved physical sector also satisfies it if its retained ground energy is within delta/16 of the unrestricted retained ground energy, with the low-energy window reduced to delta/16. The exact offset is exposed in PROOF Section 6. Comparing marker numbers and proving that offset for the desired dressed physical density remain open.

## What would suffice for c

[proved, conditional transfer theorem] PROOF Section 7 gives a weighted Schur-complement criterion. After matching a scalar shift s_M and a logical Hamiltonian K_M, it is sufficient to have a massive complementary compression and

    |F_M(z)-(K_M+s_M-z)|
        <= eta_M (K_M-E_M+v/M),     eta_M -> 0,

as two-sided quadratic-form inequalities in each fixed scaling-energy window. This implies o(1/M) errors in all fixed low-energy levels and excitation gaps. It does not require isolation or exact block diagonalization of the entire extensive code band. A low-low projected matrix element alone is insufficient because coupling to its complement can shift the levels.

[assumed, unproved application] No bound of this form is established for the chain. The estimates proved here instead give a possible relative error of order M^2 epsilon_m^2/(delta v). Proving the state-dependent improvement after the common bulk energy is removed is the missing dressing estimate. Logical RG and the CFT/current or thermal identification remain separate inputs. A freely chosen scalar s_M also cannot be used to identify a Casimir coefficient without controlling its 1/M term.

## Checks and limitations

[proved, obstruction to a shortcut] An exactly diagonalizable product of rotated auxiliary two-level systems has a positive physical gap, while its bare QHQ threshold above the true ground tends to zero exponentially with M. PROOF Section 8 gives the formula. Thus a bare complementary resolvent cannot be declared uniformly massive using the cell gap. The same example shows that vanishing vacuum overlap alone does not rule out a local unitary.

[proved, spectral-flow audit] Whole-band spectral flow lacks an established separating gap. Following only low logical levels gives gaps of order 1/M and no uniform locality estimate from the customary inverse-frequency filter. The product code also fails local indistinguishability. Neither observation rules out a special construction; both prevent importing a gapped stability argument from the cell gap alone.

[numerical, passed] The checks use matrices of dimension at most 256. They verify the new leakage, form, threshold, spectral and mixing bounds, the BCH identity to 1.24e-16, the product-rotation formula to 3.34e-16, and the iteration budget. They do not simulate the 4464-letter chain or measure c. Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py`.

[proved, files and scope] PROOF.md contains the full derivations, checks.py and checks.json the controls, and PROGRESS.md and STATE.md the handoff. No network, subagents, cluster computation or writes outside work_C49 were used. The primary Chen-Lucas-Yin theorem was not retrieved, and no chain-specific no-go is asserted from the local literature summary. The project STATE was left unchanged under the task's explicit write restriction.
