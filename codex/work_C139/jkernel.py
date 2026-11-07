"""Dual (hole) kernel J=(I+K)^{-1}K for 0<gamma<pi/2: is it pointwise positive?
Jhat(w)=(1/2)[1 - tanh(g w/2) coth((pi-g) w/2)]. Diagnostics only."""
import numpy as np, json
def Jhat(w, g):
    w = np.maximum(w, 1e-300)
    return 0.5*(1-np.tanh(g*w/2)/np.tanh((np.pi-g)*w/2))
res = {}
for delta in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999]:
    g = np.arccos(delta)
    W = 40/g
    w = np.linspace(0, W, 400001); dw = w[1]-w[0]
    jw = Jhat(w, g); jw[0] = 0.5*(1-g/(np.pi-g))
    xs = np.concatenate([np.linspace(0, 4*g, 200), np.linspace(4*g, 30, 600)])
    vals = np.array([(np.sum(jw*np.cos(w*x))-0.5*jw[0])*dw/np.pi for x in xs])
    # positive-combination check: J'(x)<=0 for x>0 ?
    dJ = np.diff(vals)
    res[delta] = dict(minJ=float(vals.min()), J0=float(vals[0]), norm_pred=float((np.pi-2*g)/(2*(np.pi-g))),
                      norm_num=float(np.trapezoid(vals, xs)*2), monotone=bool(np.all(dJ <= 1e-12)))
    print(delta, res[delta], flush=True)
json.dump(res, open('jkernel.json', 'w'), indent=1)
