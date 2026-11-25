"""Connected matrix of two quartet shapes and candidate neutral probes.

The monomial matrix is retained, with both orientations: no spatial-
symmetry or charge assumption is imposed. 'Candidate neutral' is a probe
label, not a certified emergent charge. Frame is the real ims frame.
"""
import time
import numpy as np


def operators(centered=False):
    out = {
        'Q25m': [('Sm', k) for k in (0, 2, 5, 7)],
        'Q25p': [('Sp', k) for k in (0, 2, 5, 7)],
        'Q36m': [('Sm', k) for k in (0, 3, 6, 9)],
        'Q36p': [('Sp', k) for k in (0, 3, 6, 9)],
        'Z': [('Sigmaz', 0)],
        'ZZ': [('Sigmaz', 0), ('Sigmaz', 1)],
        'Tm': [('Sp', 0), ('Sm', 1)],
        'Tp': [('Sm', 0), ('Sp', 1)],
    }
    if centered:
        for prefix,offsets in (('Q25',(-3,-1,2,4)),('Q36',(-4,-1,2,5))):
            out[prefix+'m']=[('Sm',k) for k in offsets]
            out[prefix+'p']=[('Sp',k) for k in offsets]
        out['Z1']=[('Sigmaz',1)]
    return out


def hermitian_basis(centered=False):
    """Rows map raw monomials to Hermitian quadratures, density and bonds."""
    names = ['Q25_Re', 'Q25_Im', 'Q36_Re', 'Q36_Im', 'Z', 'ZZ', 'T_Re', 'T_Im']
    T = np.zeros((8, 9 if centered else 8), complex)
    for a in (0, 2, 6):
        T[a, a:a+2] = .5
        T[a+1, a:a+2] = [-.5j, .5j]
    T[4, 4] = T[5, 5] = 1.
    if centered:
        T[4,4]=T[4,8]=.5;names[4]='Z_avg'
    return names, T


def pack(a):
    a = np.asarray(a)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def measure(psi, refs, rmin=12, rmax=96, progress=True,centered=False):
    ops = operators(True) if centered else operators()
    L = psi.L
    spans = {k: max(i for _, i in t) for k, t in ops.items()}
    lower={k:-min(i for _,i in t) for k,t in ops.items()}
    exp={}
    for k,t in ops.items():
        exp[k]=np.zeros(L,complex)
        for j in range(lower[k],L-spans[k]):
            exp[k][j]=psi.expectation_value_term([(o,j+i) for o,i in t])
    names = list(ops)
    out = {'L': L, 'chi': int(max(psi.chi)), 'refs': refs, 'operators': ops,'centered':centered,
           'expectation_valid_anchors':{k:[lower[k],L-spans[k]-1] for k in ops},
           'names': names, 'expectations': {k: pack(v) for k, v in exp.items()}, 'rows': {}}
    ts = time.time()
    for i0 in refs:
        assert i0>=max(lower.values()) and i0<L-max(spans.values())
        r = np.arange(max(rmin,max(spans.values())+max(lower.values())+1),min(rmax,L-max(spans.values())-1-i0)+1)
        if not len(r):
            continue
        js = (i0+r).tolist()
        C = np.zeros((len(r), len(names), len(names)), complex)
        for a, left in enumerate(names):
            for b, right in enumerate(names):
                raw = psi.term_correlation_function_right(ops[left], ops[right], i_L=i0, j_R=js)
                C[:, a, b] = raw - exp[left][i0]*exp[right][js]
            if progress:
                print(f'MATRIX ref={i0} left={left} elapsed={time.time()-ts:.1f}s', flush=True)
        out['rows'][str(i0)] = {'r': r.tolist(), 'C': pack(C)}
    out['seconds'] = time.time()-ts
    return out
