"""C190g: smoke test of cupy kernels used by z2dmrg_xp (compile path, complex128 linear algebra)."""
import cupy as cp
import numpy as np
a = cp.random.standard_normal((300, 200)) + 1j * cp.random.standard_normal((300, 200))
b = cp.vstack([a, a]); c = cp.hstack([a, a])
q, r = cp.linalg.qr(b)
h = a.conj().T @ a
e, u = cp.linalg.eigh(h)
U, s, V = cp.linalg.svd(a, full_matrices=False)
an = cp.asnumpy(a)
print('qr', float(cp.linalg.norm(q @ r - b)), 'eigh', float(abs(cp.asnumpy(e) - np.linalg.eigvalsh(an.conj().T @ an)).max()),
      'svd', float(abs(cp.asnumpy(s) - np.linalg.svd(an, compute_uv=False)).max()), 'abs', float(cp.abs(a).sum()))
print(cp.cuda.runtime.getDeviceProperties(0)['name'], cp.cuda.get_nvrtc_version() if hasattr(cp.cuda, 'get_nvrtc_version') else '')
