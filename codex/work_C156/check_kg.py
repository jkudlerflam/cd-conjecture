"""Check that the Kitaev-Gamma chain (alternating x,y bonds) equals a uniform one-site chain J' = J0 Q."""
import numpy as np
from ed import build
rng = np.random.default_rng(1)
for K, G in [(np.cos(0.85*np.pi), np.sin(0.85*np.pi)), (0.3, -0.7)]:
    Jx = np.array([[K, 0, 0], [0, 0, G], [0, G, 0]])
    Jy = np.array([[0, 0, G], [0, K, 0], [G, 0, 0]])
    J0 = -np.diag([K, G, G])                      # rotated x-bond after U6
    Q = np.zeros((3, 3)); Q[2, 0] = 1; Q[1, 2] = 1; Q[0, 1] = 1   # Q e_x=e_z, Q e_z=e_y, Q e_y=e_x
    Ju = J0 @ Q
    L = 6
    eo = np.linalg.eigvalsh(build(lambda j: Jx if j % 2 == 0 else Jy, lambda j: np.zeros(3), L).toarray())
    eu = np.linalg.eigvalsh(build(lambda j: Ju, lambda j: np.zeros(3), L).toarray())
    print(K, G, "max diff", np.max(np.abs(eo - eu)), "J' =", Ju.tolist())
