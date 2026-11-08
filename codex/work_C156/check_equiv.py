import numpy as np, scipy.sparse.linalg as spla
from ed import build
K = 1.0; h = 0.3
hv = -h * np.array([1, 1, 0]) / np.sqrt(2)
ex = np.zeros((3, 3)); ex[0, 0] = K
ey = np.zeros((3, 3)); ey[1, 1] = K
Ju = np.zeros((3, 3)); Ju[0, 1] = K   # K S^x_j S^y_{j+1}
L = 8
Ho = build(lambda j: ex if j % 2 == 0 else ey, lambda j: hv, L)
Hu = build(lambda j: Ju, lambda j: hv, L)
eo = np.linalg.eigvalsh(Ho.toarray()); eu = np.linalg.eigvalsh(Hu.toarray())
print("max spectral difference original vs uniform gauge, L=8 PBC:", np.max(np.abs(eo - eu)))
