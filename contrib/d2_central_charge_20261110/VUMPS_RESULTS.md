# Independent VUMPS checks

All rows are retained; only full passes enter the separate slope table. Changes are relative to the original DMRG source, including for warm retries.

| Model | χ | Stage | Pass | ΔE per site | ΔS₁ | S₁ | ξ | Tangent norm |
|---|---:|---|---|---:|---:|---:|---:|---:|
| candidate | 64 | [Initial](vumps_results/candidate_infinite_L2_P0_s17_chi64_g_one_vumps.json) | False | -5.3992e-12 | 1.0292e-06 | 1.069193456 | 98.327790 | 1.0334e-07 |
| candidate | 64 | [Warm retry](vumps_results/candidate_infinite_L2_P0_s17_chi64_g_one_vumps_vumps.json) | True | -5.4583e-12 | 1.244e-06 | 1.069193670 | 98.313955 | 1.2248e-08 |
| candidate | 96 | [Initial](vumps_results/candidate_infinite_L2_P0_s17_chi96_g_one_vumps.json) | False | -2.0522e-11 | 6.1556e-07 | 1.084107786 | 119.498940 | 5.2165e-07 |
| candidate | 128 | [Initial](vumps_results/candidate_infinite_L2_P0_s17_chi128_g_one_vumps.json) | False | -1.0588e-12 | 5.8446e-06 | 1.089839855 | 134.810480 | 4.0842e-07 |
| candidate | 128 | [Warm retry](vumps_results/candidate_infinite_L2_P0_s17_chi128_g_one_vumps_vumps.json) | False | -8.0597e-12 | 3.2659e-05 | 1.089866669 | 134.935016 | 1.0737e-06 |
| candidate | 192 | [Initial](vumps_results/candidate_infinite_L2_P0_s17_chi192_g_vumps.json) | False | -4.6341e-10 | 3.1049e-05 | 1.093320725 | 152.184524 | 1.2884e-06 |
| candidate | 256 | [Initial](vumps_results/candidate_infinite_L2_P0_s17_chi256_g_vumps.json) | False | -1.2953e-10 | 5.2233e-06 | 1.094192219 | 162.390749 | 4.5422e-07 |
| candidate | 256 | [Warm retry](vumps_results/candidate_infinite_L2_P0_s17_chi256_g_vumps_vumps.json) | False | -1.3188e-10 | 2.6259e-06 | 1.094189622 | 162.317878 | 1.8289e-07 |
| ising_critical | 24 | [Initial](vumps_results/ising_critical_infinite_L2_P0_s23_chi24_a1_vumps.json) | True | -6.2096e-10 | -0.00030566 | 0.529167129 | 307.940062 | 9.0063e-08 |

No family has two fully accepted bond dimensions, so no VUMPS slope is assigned.

Source links, failed gates and checksums are in `vumps_summary.json`; raw iteration histories are in `vumps_results/`.
