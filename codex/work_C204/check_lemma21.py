"""Lemma 2.1 / 2.2 and Theorem H, own code."""
import numpy as np
from c204lib import *
rng = np.random.default_rng(1)
worst = 0; worst_sh = 0; nonunimodal = 0
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
    D = np.concatenate([[0], D, [0]])
    TV = np.abs(np.diff(D)).sum()
    bound = abs(np.cos(2*al))*np.sin(2*eta)/PI
    if bound > 1e-14:
        worst = max(worst, TV/bound)
    worst_sh = max(worst_sh, TV/(2*np.abs(D).max()) if np.abs(D).max() > 1e-15 else 1)
    worst_sh = max(worst_sh, TV/(2*np.abs(D).max()+1e-300)) if np.abs(D).max()>1e-13 else worst_sh
    # one sign
    if D.max() > 1e-14*eta and D.min() < -1e-14*eta: nonunimodal += 1
print('Lemma 2.1: max TV/(|cos2a| sin2eta/pi) =', worst, ' max TV/(2 sup|D|) =', worst_sh, ' sign changes:', nonunimodal)
# sharpness: alpha = eta, small eta
for eta in [1e-2, 1e-3, 1e-4]:
    K = Kern(eta); b = PI/(PI-eta); al = eta
    phi = np.linspace(0, PI/2, 400001)[1:-1]
    with np.errstate(all='ignore'):
        wp = np.where(al+phi < PI/2, np.log(np.tan(al+phi)/np.tan(al))/b, np.inf)
        wm = np.where(phi < al, -np.log(np.tan(al-phi)/np.tan(al))/b, np.inf)
    Th = lambda w: np.where(np.isinf(w), K.q/2, K.Theta(np.where(np.isinf(w), 0, w)))
    D = Th(wp) - Th(wm)
    print('  alpha=eta=%g: sup|D|/eta = %.5f (1/pi = %.5f)' % (eta, np.abs(D).max()/eta, 1/PI))
# Theorem H on true roots
print('Theorem H: L sup|Nx-N| (bound 1), eta L sup|sigma-r0| / ((2/pi) eta cot eta)')
for x in [1e-3, 1e-5, 1e-7, 0.29]:
    eta = float(np.arccos(1-x)); K = Kern(eta); b = PI/(PI-eta)
    for L in [2, 4, 10, 64, 256, 1024, 2048]:
        n = L//2
        # initial guess: N(x_j) = (j-1/2)/L
        j = np.arange(1, n+1)
        x0 = np.log(np.tan(PI*(j-0.5)/L))/b
        xs, res, it = bethe_solve(eta, L, n, x0=x0)
        N = np.arctan(np.exp(b*xs))/PI
        KL = np.max(np.maximum(np.abs(j - L*N), np.abs(j-1-L*N)))
        u = np.unique(np.concatenate([xs, (xs[:-1]+xs[1:])/2, np.linspace(-20, 20, 2001)]))
        sig = K.a1(u) + K.T(u[:, None]-xs[None, :]).sum(1)/L
        r0 = b/(2*PI*np.cosh(b*u))
        DL = eta*L*np.max(np.abs(sig-r0))
        print('  1+D=%g L=%d: L sup|Nx-N| = %.4f, ratio to (D)-bound = %.4f, res %.1e' % (x, L, KL, DL/((2/PI)*eta/np.tan(eta)), res))
