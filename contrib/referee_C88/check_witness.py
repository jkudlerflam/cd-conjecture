"""(6) Spot check of C88's single-letter creation witnesses (sector_checks.json) with my13:
each input word is legal with F = 0 and lacks the letter, the output is legal with F = 0,
contains the letter, and the amplitude equals my bond entry (K = 0, lambda = 0.17).
Usage: python3 check_witness.py <path to sector_checks.json>"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my13
NAME = {('L', 0): 'L0', ('L', 1): 'L1', ('R', 0): 'R0', ('R', 1): 'R1', ('d', 'A'): 'DA', ('d', 'B'): 'DB',
        ('x', 'A'): 'XA', ('x', 'B'): 'XB', ('e', 0): 'E', ('m', 'A'): 'MA', ('m', 'B'): 'MB', ('V', 0): 'V0', ('V', 1): 'V1'}
d = json.load(open(sys.argv[1]))['deletion']['witnesses']
pm = my13.params(2, K=0.)
h, ls = my13.bond(pm, .17)
ix = {l: i for i, l in enumerate(ls)}; dd = len(ls)
res = {}
for key, wv in d.items():
    let = NAME[tuple(int(x) if x.strip(" '").isdigit() else x.strip(" '") for x in key.strip('()').split(','))]
    a = tuple(NAME[tuple(x)] for x in wv['input']); b = tuple(NAME[tuple(x)] for x in wv['output'])
    i = wv['bond']; j = (i + 1) % len(a)
    amp = h[ix[b[i]] * dd + ix[b[j]], ix[a[i]] * dd + ix[a[j]]]
    same_elsewhere = all(a[k] == b[k] for k in range(len(a)) if k not in (i, j))
    res[let] = dict(input_legal=my13.nbad(a) == 0, output_legal=my13.nbad(b) == 0, F_in=my13.F(a), F_out=my13.F(b),
                    letter_absent_in_input=let not in a, letter_present_in_output=let in b, local=same_elsewhere,
                    amplitude_match=bool(abs(amp - complex(*wv['amplitude'])) < 1e-12))
ok = all(all(v2 in (True, 0) for v2 in v.values()) for v in res.values())
print(json.dumps(dict(all_ok=ok, letters=len(res)), indent=0))
json.dump(dict(all_ok=ok, witnesses=res), open(os.path.join(HERE, 'check_witness.json'), 'w'), indent=1)
