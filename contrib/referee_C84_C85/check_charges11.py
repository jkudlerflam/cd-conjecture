"""(7) C84 PROOF Section 8: rank of the additive-charge equations of the literal
nonzero two-site rules (one equation q(in_l)+q(in_r) = q(out_l)+q(out_r) per rule),
from my11.py, and the claimed three-parameter solution (identity, marker number,
the charge (1))."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my11

pm = my11.default_params(4, g=.1, td=.2)
rows = []
for il, ir, ol, orr, amp in my11.move_rules(pm) + my11.port_rules(pm):
    r = np.zeros(my11.d)
    r[my11.IX[il]] += 1; r[my11.IX[ir]] += 1; r[my11.IX[ol]] -= 1; r[my11.IX[orr]] -= 1
    rows.append(r)
A = np.array(rows)
rank = int(np.linalg.matrix_rank(A))
ones = np.ones(my11.d)
mark = np.array([l in my11.MARKERS for l in my11.LETTERS], float)
q = np.array([my11.charge(l) for l in my11.LETTERS], float)
sol = np.array([ones, mark, q]).T
out = dict(equations=len(rows), rank=rank, nullity=my11.d - rank,
           claimed_solutions_satisfy=float(np.max(abs(A @ sol))),
           claimed_solutions_span_rank=int(np.linalg.matrix_rank(sol)))
print(out)
json.dump(out, open(os.path.join(HERE, 'check_charges11.json'), 'w'), indent=1)
