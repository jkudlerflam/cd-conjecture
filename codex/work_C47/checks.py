"""C47. Exact finite algebra and small numerical controls, no infrared inference.
[proved] Fraction assertions certify the finite rational identities described in PROOF.md.
[numerical] Floating-point checks are diagnostics, not proofs of an infinite-chain limit.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np
from scipy.linalg import eigh, expm


def zeros(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]


def solve(a, b):
    n, m = len(a), len(b[0])
    c = [list(a[i]) + list(b[i]) for i in range(n)]
    for j in range(n):
        p = next(i for i in range(j, n) if c[i][j])
        c[j], c[p] = c[p], c[j]
        q = c[j][j]
        c[j] = [x/q for x in c[j]]
        for i in range(n):
            if i != j and c[i][j]:
                q = c[i][j]
                c[i] = [x-q*y for x, y in zip(c[i], c[j])]
    return [r[n:n+m] for r in c]


def exact_bond(d, k, jp, jz):
    # Each excited half-cell has syndrome a,u,d and logical bit 0,1.
    states = [(s, c) for s in ('a', 'u', 'd') for c in (0, 1)]
    h = zeros(6, 6)
    for i in range(6): h[i][i] = d
    for x, y, v in [(('a', 0), ('u', 1), k), (('a', 1), ('d', 0), -k)]:
        i, j = states.index(x), states.index(y)
        h[i][j] = h[j][i] = v
    a = zeros(36, 36)
    for i in range(6):
        for j in range(6):
            for l in range(6):
                a[6*i+j][6*l+j] += h[i][l]
                a[6*i+j][6*i+l] += h[j][l]
    src = zeros(36, 4)
    def row(s, c, t, e): return 6*states.index((s, c))+states.index((t, e))
    for c in (0, 1):
        for e in (0, 1):
            j = 2*c+e
            src[row('a', c, 'a', e)][j] = jz
            src[row('u', c, 'd', e)][j] = jp/2
            src[row('d', c, 'u', e)][j] = jp/2
    x = solve(a, src)
    eff = [[-sum(src[r][i]*x[r][j] for r in range(36)) for j in range(4)] for i in range(4)]
    aa = (2*d*d-k*k)/(4*d*(d*d-k*k))
    bb = 2*d/(4*d*d-k*k)
    epar = -jz*jz*aa-jp*jp*bb/2
    eanti = -jz*jz*aa-jp*jp*(1/(2*d)+aa)/4
    t = jp*jz*k*k/(4*d*(d*d-k*k))
    u = 3*jp*jp*k**4/(32*d*(d*d-k*k)*(4*d*d-k*k))
    scalar = (epar+eanti)/2
    target = [[F(0) for _ in range(4)] for _ in range(4)]
    for j in range(4): target[j][j] = scalar+(u if j in (0,3) else -u)
    target[1][2] = target[2][1] = t
    assert eff == target
    # Exact third-order isolated-bond coefficient, since P V P = 0.
    # sqrt(2) times the excited-syndrome compression of sigma-plus.
    m = zeros(6, 6)
    for c in (0, 1):
        m[states.index(('u', c))][states.index(('a', c))] = -1
        m[states.index(('a', c))][states.index(('d', c))] = 1
    zvalues = [F({'a':0, 'u':1, 'd':-1}[ss]) for ss, _ in states]
    vqq = zeros(36, 36)
    for i in range(6):
        for j in range(6):
            for l in range(6):
                for n in range(6):
                    vqq[6*i+j][6*l+n] = jp*(m[i][l]*m[n][j]+m[l][i]*m[j][n])/2
            vqq[6*i+j][6*i+j] += jz*zvalues[i]*zvalues[j]
    vx = [[sum(vqq[i][j]*x[j][c] for j in range(36)) for c in range(4)] for i in range(36)]
    third = [[sum(x[r][i]*vx[r][j] for r in range(36)) for j in range(4)] for i in range(4)]
    assert third[0][0] == third[3][3] and third[1][1] == third[2][2]
    return dict(scalar=scalar, hopping=t, zz=u, anisotropy=(-2*u/t if t else None),third_scalar=(third[0][0]+third[1][1])/2,third_hopping=third[1][2],third_zz=(third[0][0]-third[1][1])/2), eff


def perm_op(op, perm):
    n = len(perm)
    inv = list(np.argsort(perm))
    return op.reshape([2]*(2*n)).transpose(inv+[n+i for i in inv]).reshape(2**n,2**n)


def numerical_checks():
    eye = np.eye(2); sp = np.array([[0.,1.],[0.,0.]])
    sm = sp.T; z = np.diag([1.,-1.]); x = sp+sm
    basis = np.eye(4)
    g = (basis[0]+basis[3])/math.sqrt(2)
    a = (basis[0]-basis[3])/math.sqrt(2)
    u, v = basis[1], basis[2]
    hp = np.kron(np.eye(4)-np.outer(g,g),eye)
    off = .6*(np.kron(np.outer(a,u),sp)-np.kron(np.outer(a,v),sm))
    hp += off+off.T
    hl = np.kron(hp,np.eye(4))
    hr = perm_op(np.kron(hp,np.eye(4)),[4,3,2,0,1])
    cell = hl+hr
    code = np.column_stack([np.kron(np.kron(g,eye[:,c]),g) for c in (0,1)])
    zz = []
    for i in range(32):
        zs = [1-2*((i>>(4-j))&1) for j in range(5)]
        zz.append((zs[0]-zs[1]+zs[2]-zs[3]+zs[4])/2)
    q = np.diag(zz)
    parity = np.diag([(-1)**int(round(.5-w)) for w in zz])
    cconj = x
    for _ in range(4): cconj = np.kron(cconj,x)
    spectrum = eigh(cell,eigvals_only=True)
    # 64-dimensional invariant two-cell space, outer half-cells in g.
    h00 = np.kron(hp,np.eye(8))+np.kron(np.eye(8),hp)
    spp = np.kron(sp,np.eye(4)); smp=spp.T; zp=np.kron(z,np.eye(4))
    vv = 2*(np.kron(spp,smp)+np.kron(smp,spp))-np.kron(zp,zp)
    gg = np.kron(np.kron(g[:,None],eye),np.kron(g[:,None],eye))
    sources = vv@gg
    ew, ev = eigh(h00)
    mask=ew>1e-10
    eff = -sources.T@(ev[:,mask]@((ev[:,mask].T@sources)/ew[mask,None]))
    exact, rat = exact_bond(F(1),F(3,5),F(2),F(-1))
    rational_error=float(np.max(np.abs(eff-np.array(rat,dtype=float))))
    residuals=[]
    for lam in (.02,.01,.005):
        w,vec=eigh(h00+lam*vv)
        overlap=gg.T@vec[:,:4]
        left,sv,right=np.linalg.svd(overlap)
        rotation=left@right
        heff=rotation@np.diag(w[:4])@rotation.T
        residue=heff-lam*lam*eff
        residuals.append(dict(lam=lam,remainder_over_lam3=float(np.linalg.norm(residue,2)/lam**3),min_overlap=float(sv.min())))
    # Independent finite-order local-unitary check in the active two-cell space.
    pp = gg@gg.T
    qq = np.eye(64)-pp
    resolvent = (ev[:,mask]/ew[mask])@ev[:,mask].T
    def diag_blocks(op): return pp@op@pp+qq@op@qq
    def comm(aa,bb): return aa@bb-bb@aa
    ss1 = resolvent@vv@pp-pp@vv@resolvent
    bd1 = diag_blocks(vv)
    offv = vv-bd1
    tt2 = comm(ss1,bd1)+.5*comm(ss1,offv)
    ss2 = resolvent@tt2@pp-pp@tt2@resolvent
    bd2 = diag_blocks(tt2)
    sw=[]
    for lam in (.0008,.0004,.0002):
        unitary = expm(lam*ss1+lam*lam*ss2)
        dd = h00+lam*bd1+lam*lam*bd2
        rem = unitary@(h00+lam*vv)@unitary.T-dd
        block_error=float(np.linalg.norm(gg.T@dd@gg-lam*lam*eff,2))
        assert block_error<1e-12
        sw.append(dict(lam=lam,code_block_error=block_error,remainder_over_lam3=float(np.linalg.norm(rem,2)/lam**3)))
    checks=dict(cell_gap=float(spectrum[2]),zero_modes=int(sum(abs(spectrum)<1e-10)),code_residual=float(np.linalg.norm(cell@code)),charge_commutator=float(np.linalg.norm(cell@q-q@cell)),conjugation_commutator=float(np.linalg.norm(cell@cconj-cconj@cell)),parity_multiplicities=[int(sum(np.diag(parity)>0)),int(sum(np.diag(parity)<0))],rational_second_order_error=rational_error,two_cell_residuals=residuals,finite_order_unitary=sw)
    assert checks['zero_modes']==2 and abs(checks['cell_gap']-.4)<1e-12
    assert rational_error<1e-12
    assert checks['charge_commutator']<1e-12 and checks['conjugation_commutator']<1e-12
    return checks


def main():
    fixtures=[(F(1),F(3,5),F(2),F(-1)),(F(2),F(1,3),F(3,2),F(-2,3)),(F(1),F(0),F(2),F(-1))]
    rows=[]
    for args in fixtures:
        r,_=exact_bond(*args)
        rows.append(dict(parameters=list(map(str,args)),coefficients={k:str(v) for k,v in r.items()}))
    d=F(27,182)
    assert rows[0]['coefficients']['anisotropy']==str(d)
    assert rows[-1]['coefficients']['hopping']=='0' and rows[-1]['coefficients']['zz']=='0'
    hprime=2+4*math.sqrt(60)
    required_weight=1+61*math.log(8*hprime*hprime)+12
    assert required_weight<768
    out=dict(exact_checks=rows,numerical=numerical_checks(),embedding=dict(auxiliary_modes=60,total_modes=61,gap_floor='1/5',input_field_weight=hprime,sufficient_spatial_weight=768,conversion_weight_with_reserve=required_weight))
    Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
