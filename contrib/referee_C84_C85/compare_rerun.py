"""(0) Compare the outputs of a scratch rerun with the committed files.
Usage: python3 compare_rerun.py <committed copy> <rerun copy> <tag>
Writes rerun_<tag>.json: per file, identical bytes, or for JSON the list of
differing leaf keys with the largest numeric difference, or for .npz the
array-by-array comparison."""
import sys, os, json
import numpy as np

A, B, tag = sys.argv[1], sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))


def leaves(x, pre=''):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from leaves(v, f'{pre}/{k}')
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, f'{pre}[{i}]')
    else:
        yield pre, x


res = {}
for f in sorted(os.listdir(A)):
    pa, pb = os.path.join(A, f), os.path.join(B, f)
    if not os.path.isfile(pa) or f.endswith('.py') or f.endswith('.md'):
        continue
    if not os.path.exists(pb):
        res[f] = 'missing in rerun'; continue
    ba, bb = open(pa, 'rb').read(), open(pb, 'rb').read()
    if ba == bb:
        res[f] = 'identical'; continue
    if f.endswith('.json'):
        la, lb = dict(leaves(json.load(open(pa)))), dict(leaves(json.load(open(pb))))
        diff = [k for k in set(la) | set(lb) if la.get(k) != lb.get(k)]
        num = [abs(la[k] - lb[k]) for k in diff if isinstance(la.get(k), (int, float)) and isinstance(lb.get(k), (int, float))
               and not isinstance(la.get(k), bool)]
        nontime = [k for k in diff if 'second' not in k.lower() and 'time' not in k.lower()]
        res[f] = dict(differing_leaves=len(diff), differing_non_timing=sorted(nontime)[:12],
                      max_abs_numeric_diff_non_timing=max([abs(la[k] - lb[k]) for k in nontime
                                                           if isinstance(la.get(k), (int, float)) and isinstance(lb.get(k), (int, float))],
                                                          default=None))
    elif f.endswith('.npz'):
        za, zb = np.load(pa, allow_pickle=True), np.load(pb, allow_pickle=True)
        same = all(k in zb.files and np.array_equal(za[k], zb[k]) for k in za.files) and set(za.files) == set(zb.files)
        res[f] = 'arrays identical (container bytes differ)' if same else 'ARRAYS DIFFER'
    else:
        res[f] = 'bytes differ'
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(HERE, f'rerun_{tag}.json'), 'w'), indent=1)
