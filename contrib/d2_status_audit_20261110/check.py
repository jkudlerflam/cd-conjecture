"""Read-only arithmetic audit of C190g's completed ground-branch estimates."""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / 'codex/work_C190/gpu/typhon_out/outG'
KEYS = ['p153c', 'sym40_n316', 'rat_n316', 'sym40_n18', 'p153_n316', 'sym80_n316']


def inspect(path):
    data = json.loads(path.read_text())
    if data.get('status') != 'done':
        return None
    hist = data['hist']
    a, b = hist[-1], hist[-2]
    previous = [h for h in hist if h['chi'] < a['chi'] and h.get('Sbar6') is not None]
    c4 = abs(a['Sbar6'] - previous[-1]['Sbar6'])
    c13 = (a['res_max'] <= 1e-8 and abs(a['E']-b['E']) <= 1e-8
           and max(abs(a[k]-b[k]) for k in ('Sbar1','Sbar6','Sbar12')) <= 1e-3
           and abs(a['n_mid']-b['n_mid']) <= 5e-4 and a['trunc'] <= 1e-7)
    return dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                parity=data['cfg']['P'], hz=data['cfg']['hz'], E=a['E_phys'],
                Sbar6=a['Sbar6'], final_canonical_Sbar6=data['Sbar6'],
                chi=a['chi'], passes_C1_C3=c13, C4_entropy_step=c4,
                passes_C4=c4 <= .01)


def main():
    rows = []
    for key in KEYS:
        selected = {}
        for length in (128, 256):
            sectors = []
            for parity in (0, 1):
                candidates = [inspect(p) for p in sorted(DATA.glob(f'{key}_L{length}_P{parity}_*.json'))]
                candidates = [d for d in candidates if d is not None]
                pool = [d for d in candidates if d['passes_C1_C3'] and d['passes_C4']]
                pool = pool or [d for d in candidates if d['passes_C1_C3']]
                if not pool:
                    raise ValueError(f'No admissible state: {key}, L={length}, parity={parity}')
                sectors.append(min(pool, key=lambda d:d['E']))
            selected[length] = min(sectors, key=lambda d:d['E'])
        a, b = selected[128], selected[256]
        rows.append(dict(key=key, L128=a, L256=b,
            c_W6=6*(b['Sbar6']-a['Sbar6'])/math.log(2),
            c_W6_final_canonical=6*(b['final_canonical_Sbar6']-a['final_canonical_Sbar6'])/math.log(2)))
    output = dict(status='Arithmetic audit only; L256 C4 still fails, no counterexample established.', rows=rows)
    (HERE/'results.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    for row in rows:
        print(row['key'], f"c_W6={row['c_W6']:.6f}",
              'parities', row['L128']['parity'], row['L256']['parity'])


if __name__ == '__main__':
    main()
