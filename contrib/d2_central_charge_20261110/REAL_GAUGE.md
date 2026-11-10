# Exact real gauge, 2026-11-10

Use the on-site unitary `U=prod exp[-i (-1)^j pi Z_j/8]`. The identical
candidate becomes real with alternating couplings
`Jx(j)=1/2+(-1)^j/10`, `Jy(j)=-1/2+(-1)^j/10`, `Jz=5/8`, `hz=-7/8`.
There is no DM term in this computational frame. U is a product of on-site
rotations, so every spatial entropy and every energy is unchanged; the
physical model still has one-site translation symmetry and qubit sites.
This is the same exact gauge used in the prior C190 open-chain solver.

The periodic L8 even result must reproduce the original-frame energy and
entropy before larger real-gauge runs are used. Real-gauge infinite runs
start independently in an alternating-X product state and provide a
different initialization and arithmetic implementation. They are kept
separate in finite-entanglement fits. They can favor a real symmetry-broken
finite-chi branch; energy comparison with the unrestricted complex branch
is required. Reported X/Y observables are rotated back to the physical frame.

The local Lanczos stopping setting is tightened to an absolute Ritz residual
target of approximately 1e-10 by using min_gap=1 and P_tol=1e-20. Actual
untruncated residuals are still measured independently, with the unchanged
acceptance threshold 1e-8. This changes no Hamiltonian or fitting rule.
