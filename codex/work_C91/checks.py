"""C91 finite controls. Run only from this folder. No parent imports."""
from pathlib import Path
from itertools import product
import json
import math
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
assert Path.cwd().resolve() == ROOT


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first,) + rest


def gap_matrix(n, total):
    states = list(compositions(total, n))
    lookup = {x: i for i, x in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for col, state in enumerate(states):
        for j in range(n):
            if state[j]:
                for k in ((j-1) % n, (j+1) % n):
                    out = list(state)
                    out[j] -= 1
                    out[k] += 1
                    h[lookup[tuple(out)], col] -= 1
    return states, h


def phase_filter_check():
    states, h = gap_matrix(4, 4)
    a = np.diag([int(x[0] == 0) for x in states])
    e, u = np.linalg.eigh(h)
    mat = u.T @ a @ u
    omega = e[:, None] - e[None, :]
    delta = 0.7
    low = np.exp(-(omega / delta)**2)
    f = np.zeros_like(omega)
    np.divide(1-low, omega, out=f, where=np.abs(omega) > 1e-12)
    generator = u @ (f * mat) @ u.T
    filtered = u @ (low * mat) @ u.T
    identity = np.linalg.norm(a + generator @ h - h @ generator - filtered, 2)
    anti = np.linalg.norm(generator + generator.T, 2)
    comm = lambda x, y: x @ y - y @ x
    second = comm(generator, a) + comm(generator, comm(generator, h)) / 2
    trials = []
    for s in (0.02, 0.01, 0.005):
        unitary = expm(s * generator)
        transformed = unitary @ (h+s*a) @ unitary.T
        error = transformed - h - s*filtered
        trials.append(dict(s=s, remainder_over_s2=float(np.linalg.norm(error, 2)/s**2),
                           second_coefficient_error=float(np.linalg.norm(error/s**2-second, 2))))
    assert identity < 1e-12 and anti < 1e-12
    assert trials[-1]['second_coefficient_error'] < trials[0]['second_coefficient_error']/3
    # The all-positive-gap basis vector has zero contact on every gap.
    no_contact = states.index((1, 1, 1, 1))
    assert sum(int(states[no_contact][j] == 0) for j in range(4)) == 0
    return dict(dimension=len(states), first_order_identity_error=float(identity),
                antihermitian_error=float(anti), trials=trials)


def clock_and_schedule():
    b = 15
    v = np.sqrt(2/(b+1))*np.sin(math.pi*np.arange(1, b+1)/(b+1))
    h = -np.diag(np.ones(b-1), 1)-np.diag(np.ones(b-1), -1)
    e0 = -2*math.cos(math.pi/16)
    residual = np.linalg.norm(h @ v-e0*v)
    h0 = np.zeros((b,b)); h0[0,1] = h0[1,0] = -1
    leakage = h0 @ v-v*np.vdot(v, h0 @ v)
    formula = math.sqrt(v[0]**2+v[1]**2-4*v[0]**2*v[1]**2)
    assert abs(np.linalg.norm(leakage)-formula) < 1e-14
    contact = .5-2/math.pi**2
    weights = (1-contact)*np.cumsum(v*v)[:13]
    assert all(weights[j+1] > weights[j] for j in range(12))
    # Evaluate the prescribed recurrence using integers, with logarithms only
    # for reporting scale comparisons that would underflow in floating point.
    d0, d1 = 1, 1004
    for _ in range(2, 16):
        d0, d1 = d1, 1004*d1-d0
    log2_unit = 4*math.log2(v[0])-2*math.log2(d1)-544
    log2_band_over_J = math.log2(contact)-1216+log2_unit
    # h_color / (J contact lambda^2 E_unit), at K=2^2560.
    offsets = []
    for exponent in (100, 250, 1000):
        logratio = math.log2(weights[-1])+2*exponent-2560-log2_band_over_J
        offsets.append(dict(J_over_t_power=-exponent, log2_color_to_reference=logratio))
    lower_for_budget = (2560+log2_band_over_J-math.log2(weights[-1])-35)/2
    return dict(clock_residual=float(residual), one_contact_code_leakage=formula,
                min_color_field=float(weights[0]), max_color_field=float(weights[-1]),
                color_field_spread=float(weights[-1]-weights[0]),
                log2_E_unit=log2_unit, log2_reference_scale_over_J=log2_band_over_J,
                K_floor_examples=offsets,
                max_exponent_of_t_over_J_for_color_budget_2neg35=lower_for_budget)


def parity_check():
    names = ('L0','L1','R0','R1','DA','DB','XA','XB','E','MA','MB','V0','V1')
    idx = {name:i for i,name in enumerate(names)}
    moves = []
    def add(left, right):
        moves.append((tuple(idx[x] for x in left), tuple(idx[x] for x in right)))
    for a in (0,1):
        add(('R0','R'+str(a)), ('L'+str(a),'R0'))
        for s in ('A','B'):
            add(('D'+s,'R'+str(a)), ('L'+str(a),'D'+s))
            add(('X'+s,'R'+str(a)), ('L'+str(a),'X'+s))
    add(('R0','E'), ('L0','E'))
    for s in ('A','B'):
        add(('D'+s,'E'), ('L0','M'+s))
        add(('X'+s,'E'), ('L0','M'+s))
        for marker in ('E','MA','MB'):
            add((marker,'D'+s), (marker,'X'+s))
    add(('DA','R1'), ('DB','R0'))
    add(('E','R0'), ('MA','DB'))
    add(('E','R0'), ('MB','DA'))
    for a in names[:11]:
        outgoing = 'V1' if a.startswith(('R','D','X')) else 'V0'
        incoming = ('V0','V1') if a in ('R0','E') else (('V1',) if a=='R1' else ('V0',))
        for color in incoming:
            add((color,a), (a,outgoing))
    solutions = []
    for signs in product((0,1), repeat=13):
        if all((sum(signs[i] for i in inp)-sum(signs[i] for i in out)) % 2 == 0
               for inp,out in moves):
            solutions.append(signs)
    assert solutions
    assert all(s[idx['L0']]==s[idx['L1']]==s[idx['R0']]==s[idx['R1']] for s in solutions)
    assert all(s[idx['DA']]==s[idx['DB']] and s[idx['V0']]==s[idx['V1']] for s in solutions)
    return dict(assignments_checked=2**13, conserving_assignments=len(solutions),
                all_assignments_bit_blind=True, literal_moves=len(moves))


def constant_controls():
    for q in (.001,.1,.25,.5):
        first = q*(1+q)/(1-q)**3
        assert first <= 12*q*(1+1e-12)
        assert first-q <= 22*q*q*(1+1e-12)
    m=2**13
    cr=2*(64*math.e)+9
    # Strong support majorant exponent when delta=s^(1/20).
    rows=[]
    for s in (.01,.001,.0001,2.0**-200):
        delta=s**.05
        radius=math.ceil((2*64*math.e+8)*delta**-2)
        q=20*m*radius*s/delta
        rows.append(dict(s=s, delta=delta, radius=radius, q=q,
                         theorem_domain=q<=.5))
    # The final row satisfies the theorem domain; the moderate rows do not.
    assert rows[-1]['theorem_domain']
    # Exact exponent budgets for the illustrative point in PROOF section 6.
    # s=2^-5000, delta=2^-250, K=2^14000, lambda=2^-608.
    assert 20*m*357 < 2**26
    assert 2**20*m*m*357**4 < 2**80
    assert 2**15*m*357**2 < 2**45
    assert 80-5000+2500 == -2420
    assert 10000-14000+2080 == -1920
    return dict(C_R=cr, log2_C_UV=20+2*math.log2(m)+4*math.log2(cr),
                log2_C_exp=28+2*math.log2(m)+4*math.log2(cr), sample_domains=rows,
                illustrative_point=dict(J_over_t_power=-5000, K_power=14000,
                    lambda_power=-608, cutoff_power=-250,
                    UV_over_J_upper_power=-2419,
                    color_over_reference_upper_power=-1920,
                    scope='Only ultraviolet first-moment and leading color budgets'),
                filter_centered_lower_at_delta_half=math.erf(1)-(.5-2/math.pi**2))


result = dict(status='all finite controls passed; mobile initialization remains open',
              phase_filter=phase_filter_check(), clock=clock_and_schedule(),
              parity=parity_check(), constants=constant_controls())
(ROOT/'checks.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
