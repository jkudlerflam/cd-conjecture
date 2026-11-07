"""iTEBD upper bounds on e0 for bonds without an exact formula. Usage: python e0_mps.py KIND ARGS
KIND = random NAME | xyz c n | model NAME"""
import sys, json, warnings, time
warnings.filterwarnings('ignore')
from common import *
from e0 import itebd_upper
import models
kind = sys.argv[1]
if kind == 'random' or kind == 'model':
    name = sys.argv[2]; h = models.get(name)['h']; key = name
elif kind == 'xyz':
    c, n = float(sys.argv[2]), int(sys.argv[3]); gam = c / n; h = xyz(1 + gam, 1 - gam, 0.5); key = f'xyz_c{c}_n{n}'
chi = 64
t = time.time()
e, e12, e21 = itebd_upper(h, chi=chi, taus=(0.1, 0.03, 0.01, 0.003, 0.001), steps=400)
out = HERE / 'results' / 'e0_mps' ; out.mkdir(exist_ok=True)
import e0 as E
write_json(out / f'{key}.json', dict(key=key, e0_up=e, e12=e12, e21=e21, chi=chi, seconds=time.time() - t, diag=getattr(E, 'LAST_DIAG', None)))
print(key, e, time.time() - t)
