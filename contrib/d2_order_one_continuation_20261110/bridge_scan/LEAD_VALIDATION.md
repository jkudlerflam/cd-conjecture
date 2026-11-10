# Validation of the r=5/2, D=1/10, h=-7/8 lead

Registered after complete L24/L48 output and before any L96 output.
Exact couplings are (Jx,Jy,Jz,Dz,hz)=(1/2,-1/2,5/8,1/10,-7/8),
q=0. Literal coefficient spread R=8.75. No coupling is retuned for size.

The original six-point promotion selected this point by the symmetric
closest-density rule. L24/L48 lower-sector entropy slopes are 1.703756
(fractional window) and 1.682601 (central three). L48 profile slopes are
2.022/2.081, with rms residuals 0.0065/0.0057. Long-distance XX/YY RMS
falls from 0.3619 to 0.2755. These data motivate a third size, not a CFT
identification. The two parity energies are both variational estimates.

Priority jobs are both exact parity sectors at L96, chi192, seed17, using
the unchanged microscopic C190 wrapper and convergence gates. L48 in both
parities then continues from chi128 to chi192. Original completed chi128
records are copied to `validation_snapshots` before continuation. This
establishes a direct finite-chi error estimate at L48. The L96 final chi
schedule alone does not establish bond-dimension convergence at L96.
The parent is independently running L48 seed23 at chi192 to test branch
reproducibility in its own destination; do not duplicate that job here.
This seed check does not change model selection.

All three sizes retain both entropy windows, all profile fits and density
oscillations, long-range XX/YY/ZZ diagnostics, sector energies, and sweep
errors. A promising lead needs both 48-to-96 entropy slopes above one
beyond numerical errors and no growing long-range-order signature; this
is still not a demonstration of conformality or an asymptotic c value.
The parent independently checks periodic finite-size spectra and variance.
The registered r3 scan and deferred controls remain pending while this
lead is validated.
