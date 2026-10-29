"""Gap-cutoff systematic on defect costs: b = 6 cells, N = 12, V = 8, J = 0.01, G = 5 against G = 8."""
import json
from dmrg_run import ground
from model import selector
mu = selector(6); out = {}
for G in (5, 8):
    e = {k: ground(l, p, 8, 0.01, mu, G=G, chis=(64,))['E'] for k, (l, p) in
         dict(good=([6, 6], [0, 0]), pair=([5, 7], [0, 0]), pin=([6, 6], [1, 0])).items()}
    out[G] = dict(pair=e['pair'] - e['good'], pin=e['pin'] - e['good'])
    print(G, out[G], flush=True)
out['shift'] = {k: out[8][k] - out[5][k] for k in ('pair', 'pin')}
print(out['shift']); json.dump(out, open('check_gapcut.json', 'w'), indent=1)
