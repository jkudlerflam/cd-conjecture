"""Model checks for the scheduled n = 2 point: precision of the bond, symmetry,
hierarchy of scales. python3 -B checks.py -> checks.json"""
import json
import numpy as np
import mpmath as mp
from model2 import build, load, SCHED
from ferm import decompose

keys = ('11', '12', '21', '22')
args = [SCHED[k] for k in ('rho', 'g', 'w', 'a', 'eta', 'v', 'D')]
z = np.load('kbar_n2_rho%s_g%s_w%s_a%s_eta%s_v%s_D%s.npz' % tuple(args))
K = {k: z['K' + k] for k in keys}
h2, h1, info, ops, T = build('M1')
vp = info['vp_raw'] * 1e12                       # in the 1e-12 units of the stored bond
nb = np.array([bin(i).count('1') for i in range(16)])
sc, h, rest = decompose(K['11'], 2)
out = dict(
    params=dict(SCHED, q=0.8, K=10.0, U=1.0, n=2),
    vp_raw=info['vp_raw'],
    double_minus_mp_over_vp=max(float(np.abs(z['Kdouble' + k] - K[k]).max()) for k in keys) / vp,
    payload_number_violation_over_vp=max(float(np.abs(K[k][i, j])) for k in keys for i in range(16) for j in range(16) if nb[i] != nb[j]) / vp,
    onsite_one_body_over_vp=float(max(np.abs(h[:2, :2]).max(), np.abs(h[2:, 2:]).max())) / vp,
    fast_hopping_over_vp=float(np.abs(h[:2, 2:]).max()) / vp,
    hopping_singular_values_over_vp=[float(x) for x in np.linalg.svd(h[:2, 2:] / vp, compute_uv=False)],
    payload_interaction_max_over_vp=float(np.abs(rest).max()) / vp,
    fermi_points=info['fermi'],
    velocity_ratio_payload_max_min=max(abs(x[2]) for x in info['fermi']) / min(abs(x[2]) for x in info['fermi']),
    length_velocity_over_mean_payload=2 * info['jl'] * np.sin(np.pi * info['nu']),
    V_over_vp=info['V'], V_over_jl=info['V'] / info['jl'], jl=info['jl'],
    mixed_traceless_over_vp_M1=info['mixed_max'],
    mixed_traceless_over_vp_M2=build('M2')[2]['mixed_max'],
    gauge_phi=info['gauge_phi'],
    scalars={k: str(z['scalar' + k]) for k in keys},
)


def pockets(T, field, nk=40000):
    """Fermi points of h(k) = T e^{ik} + h.c. + field at half filling."""
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.linalg.eigvalsh(np.array([T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k) + field for k in ks]))
    ef = np.sort(E.ravel())[nk - 1:nk + 1].mean()
    res = []
    for b in range(2):
        e = E[:, b] - ef
        s = np.sign(e)
        for i in np.where(s != np.roll(s, -1))[0]:
            j = (i + 1) % nk
            res.append(float(ks[i] - e[i] * (2 * np.pi / nk) / (e[j] - e[i])))
    return res


# C39 density readout: field mu * q^(2j+2)/S on flavour j (S = sum_{l=0}^{n+1} q^(2l)), units v_p
qv = 0.8
Ssum = sum(qv ** (2 * l) for l in range(4))
r = np.array([qv ** 2 / Ssum, qv ** 4 / Ssum])
ro = {}
for mu in (0.0, 0.5, 2.0):
    kf = pockets(T, mu * np.diag(r))
    w0 = 2 * np.pi - abs(kf[1] - kf[0])
    w1 = abs(kf[3] - kf[2])
    ro[str(mu)] = dict(fermi=kf, pocket_widths=[w0, w1], width_difference=abs(w0 - w1))
out['readout_field'] = dict(r=list(r), by_mu=ro)
# core spectrum at the scheduled parameters: code dimension and gap
import sys
sys.path.insert(0, '../verify_C44')
from core46 import Core, Params
Hc = Core(Params(2, g=float(SCHED['g']), w=float(SCHED['w']), a=float(SCHED['a']), eta=float(SCHED['eta']),
                 v=float(SCHED['v']), D=float(SCHED['D']))).H()
ev = np.linalg.eigvalsh(Hc.toarray() if hasattr(Hc, 'toarray') else Hc)
out['core_lowest_levels'] = [float(x) for x in ev[:6]]
try:
    out['dps34_vs_dps50_max_abs'] = json.load(open('checks.json'))['dps34_vs_dps50_max_abs']
except Exception:
    pass
print(json.dumps(out, indent=1, default=float))
json.dump(out, open('checks.json', 'w'), indent=1, default=float)
