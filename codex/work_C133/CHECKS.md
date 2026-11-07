# C133 checks

[numerical] Independent exact-energy and Gram integral controls passed.

{
  "route_A_charge_index": 3279.483870967742,
  "visibility_max_ratio": 1.0000000000000004,
  "gram": [
    {
      "n": 1,
      "trace": 0.00390625,
      "bound": 0.0058309037900874635,
      "ratio": 0.669921875,
      "integration_error": 3.265603120820532e-16
    },
    {
      "n": 2,
      "trace": 0.0039370226548424035,
      "bound": 0.0058309037900874635,
      "ratio": 0.6751993853054722,
      "integration_error": 2.938488583269658e-16
    },
    {
      "n": 4,
      "trace": 0.0039370688995750535,
      "bound": 0.0058309037900874635,
      "ratio": 0.6752073162771217,
      "integration_error": 2.940350432627623e-16
    },
    {
      "n": 8,
      "trace": 0.003937068899651663,
      "bound": 0.0058309037900874635,
      "ratio": 0.6752073162902602,
      "integration_error": 2.9403733039805796e-16
    },
    {
      "n": 12,
      "trace": 0.0039370688996516485,
      "bound": 0.0058309037900874635,
      "ratio": 0.6752073162902578,
      "integration_error": 2.9404225507697994e-16
    }
  ],
  "xx_energy_error": 3.1086244689504383e-15,
  "ising_energy_error": 1.7763568394002505e-14,
  "ashkin_teller_decoupled_error": 1.5987211554602254e-14,
  "c106_max_relative_error_last": 0.00020304376934143775,
  "matrixfree_controls": [
    {
      "family": "xxz",
      "L": 6,
      "energy_difference": 8.881784197001252e-16,
      "matrixfree_residual": 2.8256751426335126e-15
    },
    {
      "family": "ashkin_teller",
      "L": 4,
      "energy_difference": 5.329070518200751e-15,
      "matrixfree_residual": 9.198450688065657e-15
    },
    {
      "family": "spin1_blbq",
      "L": 4,
      "energy_difference": 3.552713678800501e-15,
      "matrixfree_residual": 3.372595760849947e-15
    }
  ],
  "rank_one_determinant_forbidden_coefficients": 6.28603239852682e-19
}

[numerical] C106 two-band free-chain test, t=1, a=2, k_F=pi/(2 ell), epsilon=2 sin(k_F)/sqrt(3). Both fermionic boundary sectors and their physical parity constraints were compared. At L a multiple of 4 ell the antiperiodic even-parity sector wins. This fixes shell contamination without tuning the Hamiltonian between L and 2L. The table uses L=64 ell.

| ell | epsilon | v/M | Q_L | A_eff/M | A_exact/M |
| --: | --: | --: | --: | --: | --: |
| 2 | 0.816496581 | 0.56694671 | 0.890738631 | 1.13412365 | 1.13389342 |
| 4 | 0.441884765 | 0.781261497 | 1.22733812 | 1.56269543 | 1.56252299 |
| 8 | 0.2252709 | 0.844047734 | 1.32594321 | 1.68824333 | 1.68809547 |
| 16 | 0.113180445 | 0.860478535 | 1.35174725 | 1.72109806 | 1.72095707 |
| 32 | 0.05665847 | 0.864635353 | 1.35827538 | 1.72940993 | 1.72927071 |
| 64 | 0.0283377698 | 0.865677682 | 1.35991232 | 1.73149414 | 1.73135536 |

[numerical] A second C106 test adjusts a with epsilon to keep v/M=1/2 at every parameter, while keeping the same commensurate Fermi points. Hence A/M=1 throughout. At L=64 ell the computed values are

| ell | epsilon | Q_L | A_eff/M |
| --: | --: | --: | --: |
| 2 | 1.05976678 | 0.785589192 | 1.00024322 |
| 4 | 1.0053723 | 0.785518239 | 1.00015289 |
| 8 | 0.618407244 | 0.785478332 | 1.00010207 |
| 16 | 0.331279452 | 0.785465691 | 1.00008598 |
| 32 | 0.168898687 | 0.785462263 | 1.00008161 |
| 64 | 0.0848772366 | 0.785461389 | 1.0000805 |
| 128 | 0.0424927992 | 0.785461174 | 1.00008023 |
