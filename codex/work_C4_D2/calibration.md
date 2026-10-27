# Calibration table

Raw estimates use the largest available L. Quoted uncertainties are conservative sensitivity scales, the larger of the window/model spread and the final-size-to-extrapolation correction. They are not confidence intervals or rigorous error bounds. The Fourier density splits the field equally between sites.

| Chain | expected c | L | entropy c_eff | entropy bulk | c_MV | c_E, raw v_T | extrapolated c_S | extrapolated c_MV | extrapolated c_E |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ising | 0.5 | 22 | 0.501281 | 0.500401 | 0.484874 | 0.502978 | 0.50000 ± 0.00128 | 0.50000 ± 0.01512 | 0.50000 ± 0.00298 |
| XX | 1.0 | 22 | 1.017161 | 1.003542 | 1.041772 | 1.023919 | 0.99864 ± 0.01852 | 1.00008 ± 0.04169 | 0.99996 ± 0.02396 |
| Heisenberg | 1.0 | 22 | 1.026860 | 1.010253 | 1.051319 | 1.029434 | 1.00412 ± 0.02274 | 1.00228 ± 0.04904 | 1.00345 ± 0.02599 |
| XXZ_0.5 | 1.0 | 22 | 1.020201 | 1.004938 | 1.046577 | 1.024942 | 0.99937 ± 0.02083 | 0.99986 ± 0.04672 | 1.00007 ± 0.02487 |

The quoted spreads do not include an unknown asymptotic correction law. In particular, logarithmic corrections in Heisenberg are not controlled by these polynomial fits. The benchmark value c=1 remains consistent with its observed drift.

All off-diagonal exchange entries vanish in these benchmarks. Couplings use Pauli matrices.

| chain | (Jxx, Jyy, Jzz) | field |
|---|---|---|
| Ising | (-1.0, 0.0, 0.0) | (0, 0, -1.0) |
| XX | (1.0, 1.0, 0.0) | (0, 0, 0) |
| Heisenberg | (1.0, 1.0, 1.0) | (0, 0, 0) |
| XXZ_0.5 | (1.0, 1.0, 0.5) | (0, 0, 0) |

The lowest-momentum velocity fails for Ising. Its q=1 estimate approaches 2.25, whereas the stress velocity approaches 2. Using the former as v would underestimate c.

The literal left-site field density also fails the largest-overlap Ising prescription. See methods.md and results/density_check.json.
