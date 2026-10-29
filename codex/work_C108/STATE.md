# C108 local state

[Proved, scope] Investigation completed on 2026-10-29. Read the project STATE.md, all background updates, REPORT Sections 9.68 through 9.77, C104, and the pairing README. All edits are local to work_C108. The parent STATE.md was not edited because the task forbids writes outside this folder.

[Proved, lattice result] A local replacement channel proves the chemical-potential interval [-4M/delta_-,4M/delta_+], with delta_+ and delta_- the charge-spectrum endpoint gaps. The distributional curvature has mass at most 4M(1/delta_-+1/delta_+), at most 8M for integer charges. This includes cusp atoms. It does not bound smooth curvature pointwise.

[Proved, conditional infrared result] The elementary Lieb-Robinson constant 6 exp(1) M bounds each Gaussian speed with nonzero local density response, assuming real-time distributional matching. Together with the non-drifting one-mode dictionary, D e''<=36 exp(2) M^2. Thus large curvature forces small D. No joint twist-density variational bound was found.

[Proved, conditional aggregate result] If Gaussian charge-coupling rows are independent, the response eigenvalues of kappa^(-1/2) D kappa^(-1/2) are v_a^2 and m<=r<=d-1. This proves A_U(1)<=6 exp(1) M(d-1) for that subclass. Nonzero rows alone do not imply independence. A one-charge continuum example has fixed D and kappa for arbitrarily many equal-speed modes.

[Proved, conventions] D=Z^T V Z/pi and kappa=Z^T V^(-1) Z/pi. In the drifting scalar model with current allowed to relax, D e''=v_R v_L, not the square of the mean speed. Locality still bounds both chiral speeds. A conventional Mott doping edge has e'' tending to zero, while the plateau contributes a cusp atom.

[Numerical] check.py and CHECKS.md contain passing small matrix and scalar checks. They do not establish lattice-to-CFT matching. No heavy numerics or cluster jobs were run.

[Proved, open tasks] A pointwise bound on smooth e''/M for a fixed integer charge remains unproved, with no microscopic counterexample supplied. More importantly for A, a microscopic aggregate bound on modes sharing the same charges is still needed. Do not replace nonzero mode weights by a rank count. Do not cite the locality result as following from Casimir asymptotics alone.
