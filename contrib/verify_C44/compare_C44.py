"""Compare the scale-separated rebuild with C44's saved bonds at C44's own s
(its s differs from ours by eigensolver roundoff in the tunnelling gap)."""
import json
import numpy as np
from precise import precise_bonds

ref = np.load('../../codex/work_C44/source_bonds.npz')
cases = json.load(open('../../codex/work_C44/source_results.json'))['cases']
out = []
for c in cases:
    B, ex = precise_bonds(c['n'], c['wvac'], s=c['s'])
    names = {k: 'n%d_w%s_t%s_%s' % (c['n'], c['wvac'], c['ratio'], k) for k in B}
    diff = max(float(np.abs(B[k] - ref[nm]).max()) for k, nm in names.items())
    d2 = B['11'].shape[0]
    tl = lambda M: M - np.trace(M) / d2 * np.eye(d2)
    E = B['22'] - B['21'] - B['12'] + B['11']
    R = {k: ref[nm] for k, nm in names.items()}
    Er = R['22'] - R['21'] - R['12'] + R['11']
    rec = dict(n=c['n'], wvac=c['wvac'], ratio=c['ratio'], s_C44=c['s'], s_ours=ex['s'] if False else None,
               max_abs_diff=diff, entry_scale=float(np.abs(B['11']).max()),
               payload_traceless=float(np.linalg.norm(tl(B['11']), 2)),
               ELR_scalar=float(np.trace(E).real / d2), ELR_scalar_C44=float(np.trace(Er).real / d2),
               ELR_traceless=float(np.linalg.norm(tl(E), 2)), ELR_traceless_C44=float(np.linalg.norm(tl(Er), 2)),
               EL_traceless=float(np.linalg.norm(tl(B['21'] - B['11']), 2)),
               EL_traceless_C44=float(np.linalg.norm(tl(R['21'] - R['11']), 2)),
               C44_mixed_density_payload_norm=c['mixed_density_payload_norm'])
    out.append(rec)
    print(json.dumps(rec))
json.dump(out, open('compare_C44.json', 'w'), indent=1)
