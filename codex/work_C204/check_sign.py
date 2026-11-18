import numpy as np
from c204lib import *
rng = np.random.default_rng(1)
for trial in range(3000):
    eta = rng.uniform(1e-4, PI/4) if trial % 2 else 10**rng.uniform(-5, np.log10(PI/4))
    al = rng.uniform(1e-5, PI/2 - 1e-5)
    K = Kern(eta); b = PI/(PI-eta)
    phi = np.linspace(0, PI/2, 20001)[1:-1]
    with np.errstate(all='ignore'):
        wp = np.where(al+phi < PI/2, np.log(np.tan(al+phi)/np.tan(al))/b, np.inf)
        wm = np.where(phi < al, -np.log(np.tan(al-phi)/np.tan(al))/b, np.inf)
    Th = lambda w: np.where(np.isinf(w), K.q/2, K.Theta(np.where(np.isinf(w), 0, w)))
    D = Th(wp) - Th(wm)
    if D.max() > 1e-14*eta and D.min() < -1e-14*eta:
        print('eta=%.3e alpha=%.6f cos2a=%.2e  max=%.2e min=%.2e  minor/major=%.2e' % (eta, al, np.cos(2*al), D.max(), D.min(), min(D.max(), -D.min())/max(D.max(), -D.min())))
