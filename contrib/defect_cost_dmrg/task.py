"""Run one task line: tag N V b J G lengths(comma) pins(comma) chis(comma). Writes out/<tag>.json per chi."""
import sys, json, os
os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
from dmrg_run import ground
from model import selector, free_energy_exact, Model
tag, N, V, b, J, G, ls, ps, chis = sys.argv[1:10]
lengths = [int(x) for x in ls.split(',')]; pins = [int(x) for x in ps.split(',')]
V = int(V); b = int(b); J = float(J); G = int(G); chis = [int(x) for x in chis.split(',')]
mu = selector(b)
os.makedirs('out', exist_ok=True)
meta = dict(tag=tag, N=sum(lengths), V=V, b=b, J=J, G=G, lengths=lengths, pins=pins, mu=mu,
            frozen=Model(lengths, pins, V, J, mu, G).frozen(), E0_exact=free_energy_exact(sum(lengths), V))
def save(hist):
    json.dump(dict(meta, hist=hist), open(f'out/{tag}.json', 'w'))
ground(lengths, pins, V, J, mu, G, chis=chis, maxsweeps=8, etol=1e-11, on_chi=save,
       log=lambda x: print(tag, x, flush=True))
