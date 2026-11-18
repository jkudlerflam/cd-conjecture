"""C200: D(phi) = Theta(w_+(phi)) - Theta(w_-(phi)) at half filling; total variation / eta, sup|D|/eta, sign changes of D'."""
import numpy as np
from mpmath import mp
def Theta(w, eta):  # int_0^w T, = (1/pi) atan(tanh w / tan eta); Theta(inf) = (pi/2 - eta)/pi = q/2
    return np.arctan(np.tanh(w)/np.tan(eta))/np.pi
def D_on(alpha, eta, N=400001):
    b = np.pi/(np.pi-eta)
    # geometric + uniform mixed grid in phi on (0, pi/2 - alpha)
    top = np.pi/2-alpha
    ph = np.unique(np.concatenate([np.geomspace(1e-14, top, N//2), np.linspace(0, top, N//2),
                                   alpha-np.geomspace(1e-14, alpha, N//4), top-np.geomspace(1e-14, top, N//4)]))
    ph = ph[(ph >= 0) & (ph <= top)]
    with np.errstate(all='ignore'):
        wp = np.where(ph < top, np.log(np.tan(alpha+ph)/np.tan(alpha))/b, np.inf)
        wm = np.where(ph < alpha, np.log(np.tan(alpha)/np.tan(np.maximum(alpha-ph, 0))) /b, np.inf)
    D = Theta(wp, eta)-Theta(wm, eta)
    return ph, D
worst = (0, None); worstsup = (0, None); maxsc = 0
for eta in [np.pi/4, 0.5, 0.2, 0.05, 1e-2, 1e-3, 1e-5]:
    for alpha in np.concatenate([np.geomspace(1e-9, np.pi/4, 60), [np.pi/4*0.999, 0.3, 0.6, 0.78]]):
        ph, D = D_on(alpha, eta)
        dD = np.diff(D)
        tv = np.abs(dD).sum()/eta; sp = np.abs(D).max()/eta
        sg = np.sign(dD[np.abs(dD) > 1e-15*eta]); sc = int((sg[1:] != sg[:-1]).sum())
        maxsc = max(maxsc, sc)
        if tv > worst[0]: worst = (tv, (eta, alpha, sc))
        if sp > worstsup[0]: worstsup = (sp, (eta, alpha))
    print(f'eta={eta:.1e} worst TV/eta so far {worst}, sup|D|/eta {worstsup}, max sign changes {maxsc}', flush=True)
