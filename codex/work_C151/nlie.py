"""C151: solve the half-filled DdV nonlinear integral equation on the line Im z = eps (0 < eps < eta/2),

  th(x) = L z0(x+i eps) - (1/2 pi i) int R(x-y) [ell(y) - ln2] dy + (1/2 pi i) int R(x-y+2 i eps) [conj ell(y) - ln2] dy,
  ell(y) = Log(1 + (-1)^n exp(2 pi i th(y))),   z0(x) = (1/pi) atan(tanh(b x/2)),  R^ = T^/(1-T^),

by damped fixed-point iteration with FFT convolutions, and compare with Bethe roots:
th against theta(x+i eps) from the roots, the real-axis counting function and its roots, and the energy
F - L e_bulk = (1/pi) Im int eps_inf'(x+i eps) ell(x) dx.  Floating point only."""
import json, sys, time
import numpy as np
from bethe_ferro import solve
from continuum import Kernel, half_filled_bulk
from ddv import DdV

def Rhat(w, eta):
    g = np.pi-eta
    w = np.asarray(w, float); aw = np.abs(w)
    out = np.empty_like(aw)
    small = aw < 1e-8
    out[small] = (1-2*eta/np.pi)/(2*eta/np.pi)
    a = aw[~small]
    A = (np.pi-2*eta)*a/2
    num = 0.5*(1-np.exp(-2*A))
    den = (1+np.exp(-g*a))*0.5*(1-np.exp(-eta*a))
    out[~small] = num*np.exp(A-g*a/2-eta*a/2)/den
    return out

def solve_nlie(eta, L, eps_frac=0.25, omega=0.5, tol=1e-11, maxit=4000, xr=None, guess=None, check=None):
    n = L//2; sgn = (-1.0)**n
    g = np.pi-eta; b = np.pi/g; ep = eps_frac*eta
    X = (np.log(L)/b if xr is None else xr)+16
    h = min(eta/8, 1/(12*L/(2*g)))
    N = int(2*X/h); N += N % 2
    x = -X+h*np.arange(N)
    M = 2*N
    w = 2*np.pi*np.fft.fftfreq(M, d=h)
    Rh = Rhat(w, eta)
    with np.errstate(over='ignore', invalid='ignore'):
        K0 = Rh                                          # R(x-y)
        K2 = np.where(Rh > 0, Rh*np.exp(np.minimum(2*ep*w, 700)), 0.0)   # R(x-y+2i eps)
        Km = np.where(Rh > 0, Rh*np.exp(np.minimum(-ep*w, 700)), 0.0)    # R(x-y-i eps)
        K2 = np.nan_to_num(K2); Km = np.nan_to_num(Km)
    # convolution (f*g)(x_i) = sum_j f(x_i-x_j) g(x_j) h, done as h*IFFT(K^ * FFT(g)) with the
    # phase of the grid origin: g on x_j = -X + j h, result on x_i: shift by the origin.
    def conv(Kh, gvals):
        G = np.zeros(M, complex); G[:N] = gvals
        Gh = np.fft.fft(G)*h*np.exp(1j*w*(-X))   # int e^{iwy} g(y) dy
        out = np.fft.ifft(Kh*Gh)/h                # (1/2pi) int e^{-iwx} ... dw on grid x = j h (j from 0)
        # ifft returns values at x = j h for j = 0..M-1 relative to origin 0: we need x_i = -X + i h
        phase_fix = None
        return out
    # simpler, exact: build the convolution with the origin handled by evaluating at x_i via phase
    def conv2(Kh, gvals):
        G = np.zeros(M, complex); G[:N] = gvals
        Gh = np.fft.fft(G)                        # sum_j g_j e^{-2 pi i m j/M}
        # want sum_j R(x_i - x_j) g_j h = h * sum_j R((i-j)h) g_j : circular conv with kernel samples
        return Gh
    # kernel samples R((k)h), k = -N..N, from its transform
    def samples(Kh):
        # R(u) = (1/2pi) int e^{-iwu} Kh(w) dw, u = k h  ->  (1/(M h)) sum_m e^{-i w_m k h} Kh(w_m)
        return np.fft.fft(Kh)/(M*h)               # index k (mod M): sum_m Kh_m e^{-2 pi i m k/M}
    S0, S2, Sm = samples(K0), samples(K2), samples(Km)
    def cconv(S, gvals):
        G = np.zeros(M, complex); G[:N] = gvals
        return (np.fft.ifft(np.fft.fft(S)*np.fft.fft(G))*h)[:N]
    z0 = lambda z: np.arctan(np.tanh(b*z/2))/np.pi
    drive = L*z0(x+1j*ep)
    th = drive.copy() if guess is None else drive + guess(x)
    ln2 = np.log(2)
    def F(t):
        ell = np.log(1+sgn*np.exp(2j*np.pi*t))
        return drive - cconv(S0, ell-ln2)/(2j*np.pi) + cconv(S2, np.conj(ell)-ln2)/(2j*np.pi)
    it = 0; err = np.inf
    if check is not None:
        tb = check(x+1j*ep)
        check_res = float(np.max(np.abs(F(tb)-tb)[np.abs(x) < X-12]))
    else:
        check_res = None
    for it in range(60 if guess is None else 0):
        new = F(th); err = np.max(np.abs(new-th))
        if err < tol or not np.isfinite(err):
            break
        th = (1-omega)*th+omega*new
    if not err < tol:
        from scipy.optimize import newton_krylov
        th = drive.copy() if not np.all(np.isfinite(th)) else th
        def res(v):
            t = v[:N]+1j*v[N:]
            r = F(t)-t
            return np.concatenate([r.real, r.imag])
        v = newton_krylov(res, np.concatenate([th.real, th.imag]), f_tol=max(tol,1e-10), maxiter=400, method='lgmres')
        th = v[:N]+1j*v[N:]
        err = np.max(np.abs(F(th)-th)); it = -1
    ell = np.log(1+sgn*np.exp(2j*np.pi*th))
    real = L*z0(x) + 2*np.real(-cconv(Sm, ell-ln2)/(2j*np.pi))
    return dict(x=x, th=th, real=real, ell=ell, it=it, err=err, ep=ep, h=h, drive=drive, check_res=check_res)

if __name__ == '__main__':
    res = []
    for eta in [float(e) for e in sys.argv[1].split(',')]:
        for L in [int(v) for v in sys.argv[2].split(',')]:
            t0 = time.time()
            n = L//2
            G, xb, _ = solve(L, n, eta)
            D = DdV(eta, L, xb)
            # continuation in eta from 0.5 (independent of the Bethe roots); check = residual of Bethe theta
            etas = [eta] if eta >= 0.5 else list(np.geomspace(0.5, eta, 1+int(np.ceil(np.log(0.5/eta)/np.log(float(sys.argv[4]) if len(sys.argv)>4 else 1.25)))))
            guess = None
            for e in etas:
                S = solve_nlie(e, L, guess=guess, check=(D.theta if e == etas[-1] else None))
                gx, gd = S['x'].copy(), (S['th']-S['drive']).copy()
                guess = lambda xx, gx=gx, gd=gd: np.interp(xx, gx, gd.real, left=0, right=0)+1j*np.interp(xx, gx, gd.imag, left=0, right=0)
            x = S['x']; sel = np.abs(x) < xb[-1]+3
            dth = np.max(np.abs(S['th'][sel]-D.theta(x[sel]+1j*S['ep'])))
            I = np.arange(1, n+1)-(n+1)/2
            xr = np.interp(I, S['real'], x)
            droot = np.max(np.abs(xr-xb))
            g = np.pi-eta; b = np.pi/g; a = np.pi*np.sin(eta)
            fp = a*b*np.sinh(b*(x+1j*S['ep']))/(2*g*np.cosh(b*(x+1j*S['ep']))**2)
            E = np.imag(np.sum(fp*S['ell'])*S['h'])/np.pi
            K = Kernel(eta)
            Fdev = K.Delta*L/4+G-L*half_filled_bulk(K)
            o = dict(eta=float(eta), L=L, it=int(S['it']), bethe_residual=S['check_res'], conv=float(S['err']), dtheta_line=float(dth), droot=float(droot),
                     LdevE_nlie=float(L*E), LdevE_bethe=float(L*Fdev), sec=round(time.time()-t0, 1))
            print(o, flush=True); res.append(o); json.dump(res, open(f"nlie_{sys.argv[3]}.json", "w"), indent=1)
    json.dump(res, open(f'nlie_{sys.argv[3]}.json', 'w'), indent=1)
