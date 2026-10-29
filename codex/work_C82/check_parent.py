"""Compare C82's occupied block to C75 using in-memory source loading.

Only the three named source files are read. No module or cache is written
outside this folder. C75's bond17 denotes its intermediate material
alphabet, so it is temporarily loaded under that name for the comparison.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
import sys
import types
from pathlib import Path
import numpy as np
import bond17 as repair


def main():
    saved = {name: sys.modules.get(name) for name in ('bond18', 'bond17', 'bond15')}
    try:
        for name in ('bond18', 'bond17', 'bond15'):
            source = Path('../work_C75')/(name+'.py')
            module = types.ModuleType(name)
            module.__file__ = str(source)
            sys.modules[name] = module
            exec(compile(source.read_text(), str(source), 'exec'), module.__dict__)
        parent, original = sys.modules['bond18'], sys.modules['bond15']
        def translated(a):
            k, s, bit = a
            if k == 'e':
                return ('M', '0', None)
            if k == 'm':
                return ('M', s, None)
            return ({'d': 'D', 'x': 'X'}.get(k, k), s, bit)
        ix = [repair.IX[translated(a)] for a in original.LET]
        pairs = [a*15+b for a in ix for b in ix]
        out = []
        for n in (2, 16):
            p = parent.params(n=n, D=6., g=.07, td=.003, K=31., musel=1.3, field=.017)
            q = repair.params(n=n, D=p['D'], Delta=p['Delta'], g=p['g'],
                              td=p['td'], K=p['K'], musel=p['musel'], field=p['field'])
            h = repair.material_bond(q, .13)
            error = np.max(abs(h[np.ix_(pairs, pairs)]-original.bond(p, .13)))
            assert error < 1e-14
            out.append(dict(n=n, occupied_block_max_error=float(error)))
        with open('parent_checks.json', 'w') as f:
            json.dump(out, f, indent=2)
        print(json.dumps(out))
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


if __name__ == '__main__':
    main()
