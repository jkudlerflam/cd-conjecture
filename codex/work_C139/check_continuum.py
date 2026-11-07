"""Check the dual-kernel continuum lemmas for 0<Delta<1. Output check_continuum.json."""
import json, numpy as np
from continuum_rep import Model, Restricted
out = []
for delta in [0.1, 0.5, 0.9, 0.99]:
    M = Model(delta)
    for Q in [0.05, 0.3, 1.0, 2.0, 4.0, 7.0]:
        if M.r0(Q)*2*M.g < 1e-9:
            print('skip underresolved', delta, Q); continue
        R = Restricted(M, Q)
        ratio = R.r/M.r0(R.X)
        v = R.epsp/R.r
        order = np.argsort(R.X)
        Var_v = float(np.sum(np.abs(np.diff(v[order]))))
        rQ, ZQ, eQ = R.edge()
        ebpp = eQ/(2*rQ*ZQ**2)
        # hole side check: s >= r0 outside
        xo = np.linspace(Q, Q+15, 400)
        so = R.s_full(xo)/M.r0(xo)
        row = dict(delta=delta, Q=Q, min_r_over_r0=float(ratio.min()), max_r_over_r0=float(ratio.max()),
                   bound_upper=1+4/np.pi, Zmin=float(R.Z.min()), Zmax=float(R.Z.max()), Zlow=1/(1+M.q),
                   VarZ=float(np.sum(np.abs(np.diff(R.Z[order])))), qJ=M.qJ,
                   max_abs_epsp_over_ab_r=float(np.max(np.abs(R.epsp)/(M.a*M.b*R.r))),
                   min_epsp_pos=float(R.epsp[R.X > 0].min()),
                   Var_v=Var_v, Var_v_bound=M.a*M.b*(2+40/np.pi),
                   ebulk2=float(ebpp), ebulk2_bound=M.a*M.b*(1+M.q)**2/2,
                   min_s_out_over_r0=float(so.min()), m=R.m)
        out.append(row)
        print(json.dumps(row), flush=True)
json.dump(out, open('check_continuum.json', 'w'), indent=1)
