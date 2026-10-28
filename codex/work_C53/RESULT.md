# C53 result

[proved] A modification at d=10 puts the frozen demonstrator inside a published theorem, provided four-site physical interactions are allowed. Retain the original cell checks and replace the PP port interaction by a direct logical-spin interaction on every BPPB block. In the PBP sector the resulting Hamiltonian is exactly

    H_exact=A+s sum_j[-Sx_j Sx_j+1-Sy_j Sy_j+1
                                      +g Sz_j Sz_j+1-h Sz_j],
    s=(9/16)lambda^2,    A_j >= (2/5)Q_j.

The logical chain acts on the bare c qubits and preserves every cell code. No massive integration or Schrieffer-Wolff remainder is needed. The alphabet remains ten letters and the physical rules are one-site translation invariant. This changes the architecture and interaction range. It does not solve the original nearest-neighbor demonstrator.

[proved] Every frozen-sector ground state is logical, and every noncode sector starts at least 2/5 above the logical ground. Uniformly in volume,

    0 <= f_log(T)-f_exact(T)
                    <= T log[1+15 exp(-2/(5T))].

Consequently any independently established logical thermal power coefficient transfers to this microscopic frozen sector. [PROOF.md](PROOF.md), Sections 3-5, gives the operator and trace arguments.

[imported theorem; proved application] Benfatto-Mastropietro [Theorem 1.1, arXiv:0907.2837](https://arxiv.org/html/0907.2837), published DOI [10.1007/s10955-010-9920-0](https://doi.org/10.1007/s10955-010-9920-0), applies with v(±1)=1/2, all other v=0, sufficiently small g, and the theorem's field counterterm fixing pF=arccos(-1/3). Its full statement, boundary convention and hypothesis ledger are in PROOF Section 4. The field equals 1/3 at g=0 and is 1/3+O(g) in this schedule. Both Fermi velocities are nonzero and the half-filled quartic umklapp resonance is absent.

The conclusion includes density and pairing exponents K and 1/K and susceptibility K/(pi v_cell). The companion [Drude-weight theorem](https://arxiv.org/html/1012.0403) gives D=v_cell K/pi. The published smallness radii are existential. No nonzero numerical value of g is certified to lie inside them.

[proved, central-charge boundary] These interacting theorems do not themselves supply the thermal T^2 coefficient, a finite-size Casimir coefficient or a complete stress-tensor scaling limit. For nonzero g, that input remains missing in this report. At g=0, h=1/3, a direct free-fermion calculation gives

    v_cell=s(2sqrt(2)/3),
    f_exact(T)=e0-pi T^2/(6v_cell)+o(T^2).

The two simple Fermi points give one free Dirac field, with c_R=c_L=1 and thermal c=1 in the chain convention. The proof includes the root-integral calculation and the current stress-tensor normalization. This is a result for the full frozen sector of the redesigned model, not merely its compressed reference.

[proved, smaller modification] If the original PP port coupling is retained, adding the diagonal BPPB term

    lambda^2 (9g/64-243/11648) Z_c,j Z_c,j+1

and field -(9/32)lambda^2 h Z_c,j gives precisely the same logical reference through degree two, plus scalar -20505 lambda^2 M/11648. Its anisotropy is g and its hopping stays -9lambda^2/32. An explicit uniform O(lambda^3) remainder bound is derived in PROOF Section 2. This modification alone does not put the actual microscopic chain inside BM's instantaneous density-interaction hypotheses.

[proved, comparison] The collaborator's [demonstrator](../../contrib/decisive_attack_20261028/demonstrator/README.md) already makes the reference anisotropy arbitrarily small by tuning Jp/Jz while preserving nearest-neighbor physical interactions. Its hopping shrinks too, and its matched action still needs an RG extension. The diagonal counterterm keeps the hopping fixed but has the same matching gap. The exact redesign removes that gap by inserting the logical bond directly. None of these observations upgrades the original nearest-neighbor chain or the mobile wrapper to an unconditional CFT theorem.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. Results are in [checks.json](checks.json).

- Independent second-order reconstruction agrees to 1.56e-15. The tuned original chain retains a cubic residual at g=0,±.02.
- The exact redesign's full 1024-dimensional two-cell matrix intertwines with the logical matrix to 2.52e-16, matches its lowest four levels to 1.32e-14, and has its next level .401148117 above the ground, consistent with the .4 bound.
- Logical ED on 8, 10 and 12 sites at h=1/3 finds particle numbers 5, 6 and 7, respectively, for g=0 and .02. Parity-correct free spectra agree to 7.11e-15.
- The free thermal estimator approaches one, reaching 1.00004455 at dimensionless T=.005.

[proved, remaining scope] Unrestricted type-sector selection at nonzero coupling, interacting thermal c, a nearest-neighbor exact realization and the wrapper's c=2 remain unproved here. All writes are inside work_C53. Parent STATE.md was not edited because it is outside the allowed write scope; [STATE.md](STATE.md) is the local handoff.
