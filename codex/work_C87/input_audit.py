"""Read all sibling JSON records, retaining only a compact local inventory."""
import json
from pathlib import Path
from fractions import Fraction

records = []
for folder in ['work_C84', 'work_C86']:
    for path in sorted((Path('..') / folder).iterdir()):
        if path.suffix not in ['.md', '.py', '.json']:
            continue
        content = path.read_text()
        record = dict(file=str(path), lines=len(content.splitlines()))
        if path.suffix == '.json':
            obj = json.loads(content)
            record['keys_or_length'] = list(obj) if isinstance(obj, dict) else len(obj)
            if path.name == 'obstruction_certificate.json':
                a = obj['single_rail_certificates']
                b = obj['two_site_source_certificates']
                decode = lambda v: Fraction(int(v[0]), int(v[1]))
                assert all(decode(x['hermitian_minor_determinant_upper']) < 0 for x in a)
                assert all(decode(x['determinant_upper']) < 0 for x in b)
                assert all(decode(x['positive_and_negative_sign_gap_lower']) > 0 for x in a+b)
                record['negative_determinants_checked'] = len(a)+len(b)
            if path.name == 'finite_budget_exact.json':
                val = {k: Fraction(int(v[0]), int(v[1])) for k,v in obj.items()}
                assert val['conversion']+val['penalty']+val['lambda_local'] < val['budget']/8
                record['budget_checked'] = True
        records.append(record)
Path('input_inventory.json').write_text(json.dumps(records, indent=2)+'\n')
print(json.dumps(records, indent=2))
