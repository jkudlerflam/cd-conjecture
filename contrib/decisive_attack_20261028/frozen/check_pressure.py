"""Conservative rational witnesses for the frozen defect-pressure theorem.

Only writes pressure_witnesses.json in this folder. No large diagonalization.
The elementary interval formulas are reproduced from the existing selector
certificate with attribution below. Standard library only. The theorem is
analytic, not inferred from these checks.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
# Same elementary rational interval formulas as
# contrib/counterexample_campaign/route_a/check_global_sectors.py.
def interval_add(x, y):
    return x[0]+y[0], x[1]+y[1]

def interval_scale(c, x):
    z = c*x[0], c*x[1]
    return min(z), max(z)

def interval_sub(x, y):
    return interval_add(x, interval_scale(-1, y))

def atan_interval(x, terms):
    s = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    nxt = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(s, s+nxt), max(s, s+nxt)

def cos_interval_at(x, terms=35):
    s = sum(((-1)**j*x**(2*j)/F(math.factorial(2*j))
             for j in range(terms)), F(0))
    nxt = (-1)**terms*x**(2*terms)/F(math.factorial(2*terms))
    return min(s, s+nxt), max(s, s+nxt)

def cos_interval(x):
    assert 0 < x[0] <= x[1] < 2
    return cos_interval_at(x[1])[0], cos_interval_at(x[0])[1]

pi = interval_sub(interval_scale(16, atan_interval(F(1, 5), 65)),
                    interval_scale(4, atan_interval(F(1, 239), 22)))
aa = {ell: interval_sub((F(2), F(2)),
                         interval_scale(2, cos_interval(
                             interval_scale(F(1, ell+1), pi))))
      for ell in (23, 24, 25)}
delta_next = interval_sub(interval_scale(25, aa[24]),
                            interval_scale(24, aa[25]))
delta_prev = interval_sub(interval_scale(24, aa[23]),
                            interval_scale(23, aa[24]))
selector = interval_scale(F(1, 2), interval_add(delta_prev, delta_next))
eb = interval_scale(F(1, 24), interval_sub(aa[24], selector))
gamma = interval_scale(F(1, 48), interval_sub(delta_prev, delta_next))

def floor_dyadic(x, bits=128):
    return F((x.numerator * 2**bits) // x.denominator, 2**bits)

a = floor_dyadic(-eb[1])
gamma_lower = floor_dyadic(gamma[0])
# sqrt(3) < 2, with every 0 < |lambda| <= 2^-20.
B = F(2, 2**20)
delta_w = gamma_lower - 3*B
c0 = 2+3*B
kappa = floor_dyadic(a*delta_w/(a+c0+delta_w))
delta_i = F(10)-4-5*B
delta = min(2*kappa, delta_i)
assert min(a, delta_w, kappa, delta_i, delta) > 0
assert gamma_lower > 5*B
assert F(10)>4+7*B
assert kappa <= a
# A rational sufficient threshold: ln(2d) < 4 for d=23 because e^4>46.
# The first eight nonnegative series terms already establish the latter.
e4_lower = sum((F(4)**j/F(math.factorial(j)) for j in range(9)), F(0))
assert e4_lower > 46
beta0 = 4/kappa
# For all beta>=beta0, d exp(-beta*kappa) <= 1/2 and likewise for a.

# Independent rational checks of the two-line convex-envelope inequality.
# The analytic proof covers every ell; these cover the crossover and tails.
ells = sorted(set(range(1, 200)) | {int((c0+delta_w)/a)+j for j in range(-3,4)}
              | {10**j for j in range(1,8)})
for ell in ells:
    if ell>=1:
        assert max(delta_w, a*ell-c0) >= kappa*(ell+1)

result = {
    "status": "exact_rational_witnesses_and_analytic_all_length_tail",
    "b":24, "d":23, "K":10, "lambda_max":"2^-20",
    "free_energy_error_prefactor":2*24*23+4*24,
    "rational_lower_bounds":{key:str(value) for key,value in {
        "a":a, "gamma":gamma_lower, "delta_W":delta_w,
        "kappa":kappa, "delta_I":delta_i, "delta":delta}.items()},
    "decimal_displays":{key:float(value) for key,value in {
        "a":a, "gamma":gamma_lower, "delta_W":delta_w,
        "kappa":kappa, "delta_I":delta_i, "delta":delta,
        "sufficient_beta_min":beta0, "sufficient_T_max":1/beta0}.items()},
    "checked_envelope_lengths":len(ells),
    "note":"Decimal displays are not the proof. All positivity assertions use Fraction arithmetic."
}
(HERE/"pressure_witnesses.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({key:value for key,value in result.items()
                  if key!="rational_lower_bounds"},indent=2))
