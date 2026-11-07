"""Monotone bracketing for PROOF Theorem 4.4. Phi(lam) <= lam*[2coth(lam) - piJ0_low(lam)] / (2 pi sin(pi/8) (1-2exp(-pi^2/lam))^2),
piJ0_low = max(ln2/lam, 1/2 + 2 sum_{j<=3} t^j/(1+t^j)), t = exp(-2 lam).
On a cell [l0,l1]: lam <= l1, coth(lam) <= coth(l0), piJ0_low(lam) >= second branch at l1 (decreasing in lam),
(1-2exp(-pi^2/lam)) >= value at l1. Also direct margin tanh(lam/2)-coth(lam)/2 (increasing in lam) at lam_d."""
import numpy as np, json
c = 2*np.pi*np.sin(np.pi/8)
def br2(l):
    t = np.exp(-2*l); return 0.5+2*sum(t**j/(1+t**j) for j in (1, 2, 3))
f_small = lambda l: (2*l/np.tanh(l)-np.log(2))/(c*(1-2*np.exp(-np.pi**2/l))**2)  # increasing in l
lam_s, lam_d = 1.0, 1.36
cells = np.linspace(lam_s, lam_d, 37)
worst = max(l1*(2/np.tanh(l0)-br2(l1))/(c*(1-2*np.exp(-np.pi**2/l1))**2) for l0, l1 in zip(cells[:-1], cells[1:]))
out = dict(f_small_at_lam_s=f_small(lam_s), worst_cell_bound_1_to_1p36=worst,
           direct_margin_at_lam_d=float(np.tanh(lam_d/2)-0.5/np.tanh(lam_d)),
           Phi_limit_lam_to_0=(2-np.log(2))/c)
print(out); json.dump(out, open('pos_rigorous.json', 'w'), indent=1)
