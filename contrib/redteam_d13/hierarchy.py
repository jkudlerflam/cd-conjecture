"""Attack (i), (ii), (iii): scale ledger, operator dimensions, compressibility.

All logarithms are base 2 and computed from exact integers where possible.
Inputs are C88's schedule (PROOF eq. (16)) and the C91/C93 colour weights.
Energies are in units of J unless stated.
"""
import json
import math
from fractions import Fraction
import numpy as np

out = {}
# ---- exact beta = 1/d_15 ----
d = [1, 1004]
for j in range(1, 15):
    d.append(1004 * d[-1] - d[-2])
d15 = d[15]
lg = lambda x: math.log2(x)
alpha_cell = math.sqrt(2 / 16) * math.sin(math.pi / 16)
log_beta = -math.log2(d15)
log_g, log_lam, log_K = -272, -608, 2560
log_Eunit = 4 * lg(alpha_cell) + 2 * log_beta + 2 * log_g
alpha_c = 0.5 - 2 / math.pi ** 2
log_Wref = lg(alpha_c) + 2 * log_lam + log_Eunit       # J alpha_c lambda^2 E_unit
log_tauW = 35.2 * math.log2(10)                          # referee_C88 Part III item 8
allow_rel_raw = -35 - lg(128 * 13) - log_tauW            # PROOF (30) tau / W_ref
out['scales_log2'] = dict(
    beta=log_beta, E_unit=log_Eunit, W_ref_over_J=log_Wref,
    lambda2=2 * log_lam, zero_port_split_2p84_over_K=84 - log_K,
    vacuum_tail_J_over_K=-log_K, grammar_port_2Jlam2_over_K=1 + 2 * log_lam - log_K,
    finite_g_relative_g2_over_beta2=2 * log_g - 2 * log_beta,
    material_allowance_relative_to_Wref=allow_rel_raw,
)
s = out['scales_log2']
out['ratios_to_Wref_log2'] = {
    'zero_port_split': s['zero_port_split_2p84_over_K'] - log_Wref,
    'vacuum_tail_A_term': s['vacuum_tail_J_over_K'] - log_Wref,
    'colour_field_t2_over_JK_at_t_eq_J': -log_K - log_Wref,
    'lambda2_scalar_part': 2 * log_lam - log_Wref,
}
# colour field: t^2/(JK) relative to W_ref is (t/J)^2 2^(-2560 - log_Wref)
rel_colour0 = -log_K - log_Wref
# largest t/J such that colour field (max weight 0.686) stays below the raw allowance
w_max = 0.6864370767
out['colour_field'] = dict(
    max_weight=w_max,
    log2_tJ_max_for_raw_allowance=(allow_rel_raw - rel_colour0 - lg(w_max)) / 2,
    log2_tJ_max_for_normalised_2m35=(-35 - rel_colour0 - lg(w_max)) / 2,
    C93_K_condition='K >= max(2^2560, 2^2243 (t/J)^2)')

# ---- (ii) operator dimensions at material density rho_m with b = 15 ----
# uniform lock-in of a cell harmonic Q_m = 2 pi ell/15 needs ell rho_m/15 integer;
# charge dimension (ell/15)^2 K_f (K_f = 1 for the free vacancy fluid).
lock = []
for ell in range(15, 61):
    for h in (1, 2):
        rho = Fraction(15 * h, ell)
        if 0 < rho <= 1:
            lock.append(dict(ell=ell, h=h, rho_m=str(rho), rho_m_float=float(rho),
                             vacancies_per_cell=float(15 / rho - 15),
                             dim_free=(ell / 15) ** 2))
lock.sort(key=lambda r: r['dim_free'])
out['cell_lockin'] = dict(
    note='uniform cell-harmonic terms; relevant (dim<2) only for rho_m >= 15/21',
    relevant=[r for r in lock if r['dim_free'] < 2],
    at_half_occupancy=[r for r in lock if r['rho_m'] == '1/2'])
# Luttinger parameter at small J/t: mean contact energy per bond acts as an XXZ
# nearest-neighbour V = J e_b; Delta_xxz = V/(2t).
mu_sel = 2 - 0.5 * sum((l + 1) * 2 * (1 - math.cos(math.pi / (l + 1)))
                       - l * 2 * (1 - math.cos(math.pi / (l + 2))) for l in (14, 15))
E15 = mu_sel - 2 * math.cos(math.pi / 16)
e_b = E15 / 15
out['luttinger'] = dict(mu_sel=mu_sel, E15=E15, mean_bond_energy=e_b,
                        K_f_formula='pi/(2(pi - arccos(Delta))), Delta = J e_b/(2t)')
for Jt in (1e-3, 5e-3, 1e-1, 1.0):
    De = Jt * e_b / 2
    Kf = math.pi / (2 * (math.pi - math.acos(De)))
    out['luttinger'][f'K_f(J/t={Jt})'] = Kf
    out['luttinger'][f'umklapp_dim_4Kf(J/t={Jt})'] = 4 * Kf

# ---- (iii) compressibility / phase separation, mean field in the contact energy ----
# f(rho) = -(2t/pi) sin(pi rho) + J e_b [rho^2 - sin^2(pi rho)/pi^2]
# f'' = 2 pi t sin(pi rho) + 4 J e_b sin^2(pi rho): convex iff J|e_b| sin(pi rho) < pi t/2
out['phase_separation'] = dict(J_over_t_instability_mean_field=math.pi / (2 * abs(e_b)),
                               attractive=e_b < 0)
# marker-vacancy binding scale (marker contact mu_sel/2 per adjacent bond)
out['marker_contact_half_mu'] = mu_sel / 2
json.dump(out, open('hierarchy.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
