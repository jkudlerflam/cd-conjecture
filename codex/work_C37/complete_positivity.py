"""Reuse C33 certificates; fill only its previously uncertified sizes.
The positivity algorithm and proof belong to C33, not C37.
"""
from pathlib import Path
import json,math
from positivity_C33 import certificate
out=[]
for n in range(1,24):
    src=Path('../work_C33')/f'certificate_n{n}.json'
    if src.exists():
        r=json.loads(src.read_text());origin='inherited C33 certificate'
        assert r['D']==1000 and r['status']=='proved_positive'
    else:
        r=certificate(n);origin='C33 certifier evaluated at missing size by C37'
        json.dump(r,open(f'positivity_n{n}.json','w'),indent=2)
    assert r['status']=='proved_positive'
    lower=math.log10(2)+2*r['log10_beta']+r['log10_lambda_min_lower']
    # This displayed logarithm summarizes the inherited exact integer test.
    out.append(dict(n=n,origin=origin,source=str(src) if src.exists() else f'positivity_n{n}.json',
        log10_normalized_hopping_lower=lower,certificate=r))
    json.dump(out,open('positivity_summary.json','w'),indent=2)
    print(n,origin,lower,flush=True)
