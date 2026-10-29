"""C93 independent small controls. Run only in this directory."""
from pathlib import Path
from itertools import product, combinations
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parent
assert Path.cwd().resolve() == ROOT

MATERIAL = ('L0', 'L1', 'R0', 'R1', 'DA', 'DB', 'XA', 'XB', 'E', 'MA', 'MB')
MARKERS = ('E', 'MA', 'MB')
DATA = MATERIAL[:8]


def charge(a):
    if a in ('L1', 'R1'):
        return 2
    if a in ('DA', 'XA', 'MA'):
        return -1
    if a in ('DB', 'XB', 'MB'):
        return 1
    return 0


def excitation(a):
    return int(a[0] in 'DXM')


def rules():
    ans = []
    def add(x, y):
        ans.extend([(x, y), (y, x)])
    for bit in '01':
        add(('R0', 'R'+bit), ('L'+bit, 'R0'))
        for s in 'AB':
            for h in 'DX':
                add((h+s, 'R'+bit), ('L'+bit, h+s))
    add(('R0', 'E'), ('L0', 'E'))
    for s in 'AB':
        for h in 'DX':
            add((h+s, 'E'), ('L0', 'M'+s))
        for m in MARKERS:
            add((m, 'D'+s), (m, 'X'+s))
    add(('DA', 'R1'), ('DB', 'R0'))
    add(('E', 'R0'), ('MA', 'DB'))
    add(('E', 'R0'), ('MB', 'DA'))
    return ans


def block_grading():
    checks = []
    for lengths in ((2, 2), (1, 3)):
        cell = [j for j, ell in enumerate(lengths) for _ in range(ell)]
        markers = {lengths[0]-1, sum(lengths)-1}
        choices = [MARKERS if i in markers else DATA for i in range(sum(lengths))]
        states = list(product(*choices))
        def values(w):
            q = sum(map(charge, w))
            s = sum((-1)**cell[i]*excitation(a) for i, a in enumerate(w))
            ns = [(charge(a)-(-1)**cell[i]*excitation(a))//2 for i, a in enumerate(w)]
            assert all((charge(a)-(-1)**cell[i]*excitation(a)) % 2 == 0
                       for i, a in enumerate(w))
            assert 2*sum(ns) == q-s
            return q, s, sum(ns)
        count = 0
        for w in states:
            old = values(w)
            for i in range(len(w)):
                j = (i+1) % len(w)
                for src, dst in rules():
                    if (w[i], w[j]) != src:
                        continue
                    z = list(w)
                    z[i], z[j] = dst
                    assert values(z) == old
                    count += 1
        checks.append(dict(lengths=lengths, states=len(states), directed_moves=count))
    # Code parity at n=13, independently of the vacuum clock position.
    for bits in product((0, 1), repeat=13):
        q = 2*sum(bits)
        assert (-1)**(q//2) == (-1)**sum(bits)
    # An odd-cell seam creates two carriers with the same staggered sign.
    assert ((charge('MA')-excitation('MA'))//2
            +(charge('DB')-excitation('DB'))//2) == -1
    return dict(sectors=checks, code_strings=2**13, odd_seam_changes_integer_charge=True)


def jw_control():
    eye = np.eye(2)
    z = np.diag([1., -1.])
    c = np.array([[0., 1.], [0., 0.]])
    def at(op, j):
        out = np.ones((1, 1))
        for k in range(4):
            out = np.kron(out, op if k == j else eye)
        return out
    p = [at(z, 0), at(z, 2)]
    pi = [at(z, 1), at(z, 3)]
    P = [p[j]@pi[j] for j in range(2)]
    S = [np.eye(16), P[0]]
    psi = [S[j]@at(c, 2*j) for j in range(2)]
    xi = [S[j]@p[j]@at(c, 2*j+1) for j in range(2)]
    errors = []
    for i in range(2):
        for j in range(2):
            errors.append(np.max(abs(psi[i]@psi[j]+psi[j]@psi[i])))
            errors.append(np.max(abs(psi[i]@psi[j].T+psi[j].T@psi[i]
                                     -(np.eye(16) if i == j else 0))))
            errors.append(np.max(abs(psi[i]@xi[j]+xi[j]@psi[i])))
    old_hop = at(c.T, 0)@p[0]@at(c, 2)
    ext_hop = psi[0].T@psi[1]
    expected = old_hop@pi[0]
    errors.append(np.max(abs(ext_hop-expected)))
    errors.append(abs(np.linalg.norm(ext_hop, 2)-np.linalg.norm(old_hop, 2)))
    assert max(errors) < 1e-13
    return dict(dimension=16, maximum_error=float(max(errors)))


def sea_controls():
    rho = .5
    def G(r):
        return rho if r == 0 else math.sin(math.pi*rho*r)/(math.pi*r)
    def density_product(sites):
        sites = sorted(set(sites))
        return float(np.linalg.det([[G(i-j) for j in sites] for i in sites]))
    p2 = density_product((0, 1))
    alpha = .5-2/math.pi**2
    assert abs(p2/rho-alpha) < 1e-15
    rows = []
    for ell in (8, 32, 128):
        cov = [density_product((0, 1, r, r+1))-p2*p2 for r in range(ell)]
        for r in range(2, ell):
            assert abs(cov[r]) <= 24/(math.pi**2*(r-1)**2)
        variance = (ell*cov[0]+2*sum((ell-r)*cov[r] for r in range(1, ell)))/ell**2/rho**2
        assert -1e-14 <= variance <= 44/ell
        rows.append(dict(block_length=ell, normalized_contact_variance=variance,
                         proved_upper_bound=44/ell))
    def contact_of_gap(u):
        r = 1/(1+u)
        return r-math.sin(math.pi*r)**2/(math.pi**2*r)
    h = 1e-5
    derivative = (contact_of_gap(1+h)-contact_of_gap(1-h))/(2*h)
    exact = -.25-1/math.pi**2
    assert abs(derivative-exact) < 1e-9
    v = np.sqrt(1/8)*np.sin(math.pi*np.arange(1, 16)/16)
    weights = (1-alpha)*np.cumsum(v*v)[:13]
    assert all(np.diff(weights) > 0)
    return dict(alpha=alpha, derivative_in_gap_density=exact,
                derivative_error=abs(derivative-exact), blocks=rows,
                color_min=float(weights[0]), color_max=float(weights[-1]))


def constants():
    L = 8*(1+math.log(32)+10)
    a13 = L+13*math.log(128)+1
    a14 = L+14*math.log(128)+1
    C0 = 3*8**14
    T = max(2*L, math.log(C0)+2)
    assert T == 2*L and a14 < 256
    # The actual fixed-cutoff constants are intentionally kept symbolic.
    # Their logarithms fit in floating point; the constants themselves do not.
    a0 = 256
    mu = 4*a0+4
    logv = math.log(64)+mu
    logB = math.log(4/mu)+logv  # Neglected +17 is reported by a separate bound.
    assert 17 < math.exp(10)  # B < 2*(4 v/mu), since v is vastly larger.
    return dict(L=L, a_payload=a13, a_hypothetical_fourteen_mode=a14,
                finite_material_dimension=C0, conditional_T=T,
                conditional_kappa_over_gap=1/32,
                conditional_log10_Adec=(22*math.log(2)+L+math.log(16386))/math.log(10),
                material_support_rate=a0, LR_weight=mu, log_LR_velocity=logv,
                log_B_lower=logB, log_B_upper=logB+math.log(2),
                color_schedule_power=2243)


if __name__ == '__main__':
    out = dict(status='all small controls passed; mobile massive matching remains open',
               block_grading=block_grading(), jordan_wigner=jw_control(),
               sea=sea_controls(), constants=constants())
    (ROOT/'checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))
