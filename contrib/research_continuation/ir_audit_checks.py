"""Finite momentum checks for the graph target plus its length mode."""
import json
from pathlib import Path
import numpy as np


def check(n):
    k=np.pi/4+np.arange(1,n+1)*np.pi/(12*(n+1))
    dist=lambda x:np.abs((x+np.pi)%(2*np.pi)-np.pi)
    mixed=min(float(np.min(dist(2*sign*k+2*np.pi*nu)))
              for sign in [-1,1] for nu in np.linspace(5/12,7/12,101))
    assert mixed>=np.pi/6-1e-13
    differences=dist(2*(k[:,None]-k[None,:]))
    differences[np.diag_indices(n)]=np.inf
    payload=float(np.min(differences))
    expected=np.pi/(6*(n+1))
    assert abs(payload-expected)<1e-13
    return dict(payload_modes=n,payload_margin=payload,
                expected_payload_margin=expected,
                sampled_mixed_margin=mixed,proved_mixed_lower_bound=float(np.pi/6))


if __name__=='__main__':
    result={'scope':'Momentum inequalities only, not microscopic or thermodynamic CFT verification.',
            'cases':[check((p+1)*(p*p+p+1)) for p in [2,3,5,7]]}
    Path(__file__).with_name('ir_audit_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
