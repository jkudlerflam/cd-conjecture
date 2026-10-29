"""The conditional d=12 member and its scale ledger (log2 units of J).
Material schedule: C84's certified n=11 member (p=12, b=13), unchanged.
Mobile parameters: t, w and the vacancy filling nu (vacancies per inter-cell gap).
Constraint values use the p=12 perturbative coefficients in perturbative_scales.json.
Every constraint below is a [heuristic] order-of-magnitude requirement."""
import os, json
from math import log2, cos, sin, pi, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
ps = json.load(open(os.path.join(HERE, 'perturbative_scales.json')))
seg = {(r['p'], r['w']): r for r in ps['segment']}
tun = {(r['p'], r['w']): r for r in ps['tunnelling']}
p = 12
X4 = seg[(p, 0.)]['X4']; b = seg[(p, -1e-3)]['X2'] / 1e-3
U2 = seg[(p, 0.)]['U2']; BL2 = seg[(p, 0.)]['BL2']; c12 = tun[(p, 0.)]['t_eff_coefficient']
# C84 schedule
d = [1, 1004]
for _ in range(12):
    d.append(1004 * d[-1] - d[-2])
lb_beta = -log2(d[13])
alpha = sqrt(2 / 14) * sin(pi / 14)
lg, ll = -216, -480
lE = 4 * log2(alpha) + 2 * lb_beta + 2 * lg
lW = 2 * ll + lE
entry_gap = 2 * cos(pi / 14) - 2 * cos(pi / 13)   # clock cost of the first vacancy entry


def ledger(lt, lw, tol_log2):
    X = log2(b) + lw + 2 * lt
    Xatt = log2(abs(X4)) + 4 * lt
    teff = log2(c12) + 13 * lt
    port_scalar = 2 * ll + 4 * log2(alpha)
    nnn = port_scalar + lw + 2 * log2(alpha) + 4        # walker distortion of the port amplitude (estimate)
    drel = log2(abs(2 ** lw) * alpha ** 2 + BL2 * 2 ** (2 * lt))
    return dict(log2_t=lt, log2_minus_w=lw, log2_entry_gap=log2(entry_gap),
                log2_X_repulsive=X, log2_X_attractive_t4=Xatt,
                log2_U=log2(U2) + 2 * lt, log2_t_eff=teff, log2_port_scalar_per_contact=port_scalar,
                log2_NNN_estimate=nnn, log2_W_payload=lW, log2_relative_clock_distortion=drel,
                tolerance_log2=tol_log2,
                ok_perturbative=lt < log2(entry_gap) - 3,
                ok_repulsive=X > Xatt + 3,
                ok_tolerance=drel < tol_log2,
                ok_hopping_beats_NNN=teff > nnn + 3,
                ok_hopping_beats_payload=teff > lW + 3,
                ok_hard_rods=X > teff + 3)


member = dict(
    tag='definition_conditional_mobile_family',
    physical_dimension=12, material_letters=11, vacancy_colors=1,
    vacancy_automaton='o(V)=0, I(V)={0}', vacancy_hop_letters=['L_0', 'L_1', 'E', 'M_A', 'M_B'],
    vacancy_context_switch='-J(|V X_s><V D_s| + h.c.)', contact_energy='w on the ordered pair (V, R_0)',
    onsite='genuine on-site material energies, o(V)=0',
    material_schedule='C84 MEMBER.json unchanged (p=12, b=13, n=11, K=2^2048, lambda=2^-480, g=2^-216)',
    active_bits_n=11, conditional_c=12, d_minus_1=11,
    pin='F=0 sector of the pinless grammar, where the material operator equals C84 (H_C84 = H_pinless + K F)',
    vacancy_filling='nu = vacancies per inter-cell gap, 0 < nu < 1/2',
    perturbative_coefficients=dict(X4=X4, X2_over_minus_w=b, U2=U2, BL2=BL2, t_eff_coefficient=c12),
    windows=dict(normalized_allowance=ledger(-30, -40, log2(2 ** -35 / 8)),
                 raw_tolerance=ledger(-88, -158, -163)),
    status='finite bond facts proved or numerical; mobile reduction heuristic; infrared assumed')
print(json.dumps(member['windows'], indent=1))
json.dump(member, open(os.path.join(HERE, 'MEMBER.json'), 'w'), indent=1)
