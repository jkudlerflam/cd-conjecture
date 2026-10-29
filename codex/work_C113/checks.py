"""Small deterministic algebra and quadrature checks for C113.
No many-body diagonalization, cluster work, or fitted asymptotic claims.
"""
import json
import numpy as np


def profile(R, eps=0.4):
    x = np.arange(R)
    w = 1 + eps * np.cos(2 * np.pi * (x + 0.5) / R)
    b = w - np.roll(w, 1)
    a = (w + np.roll(w, 1)) / 2
    G = np.sum(b * b / w)
    return x, w, a, b, G


def main():
    eps, g, lam, J = 0.4, 1.0, 0.7, 1.0
    rows = []
    for R in (32, 128, 512, 2048):
        x, w, a, b, G = profile(R, eps)
        # Rationalized expression avoids cancellation in the exact product-spin energy.
        D = np.sum(lam**2 * b**2 / (np.sqrt(g*g*a*a + lam*lam*b*b) + g*a))
        k = 2 * np.pi / R
        nz = np.sin(k*x)
        exchange = J * np.sum(w) * (1 - np.cos(k)) / 4
        field = lam * np.dot(b, nz) / 2
        witness = -exchange - field
        formula = lam*eps*R*np.sin(np.pi/R)/2 - J*R*(1-np.cos(k))/4
        assert abs(witness-formula) < 1e-11
        rows.append(dict(R=R, G=float(G), gapped_Def_over_G=float(D/G),
                         ferromagnet_Def_lower=float(witness),
                         ferromagnet_lower_over_G=float(witness/G)))
    # CFT circle check for w=1+eps cos(2pi x/L), L=1.
    x = np.arange(4096)/4096
    w = 1+eps*np.cos(2*np.pi*x)
    wp = -2*np.pi*eps*np.sin(2*np.pi*x)
    I = np.mean(wp*wp/w)
    ell_eff = np.mean(1/w)
    # The n=1 global conformal deformation leaves the circle vacuum unchanged.
    circle_D = I/(24*np.pi) + np.pi/6*(1/ell_eff-np.mean(w))
    assert abs(circle_D) < 1e-13
    # On-site gapped bond has optimized M=g by extremal product witnesses.
    X=np.array([[0,1],[1,0]],float)
    Z=np.diag([1.,-1.]); Id=np.eye(2)
    h=-g/2*(np.kron(Z,Id)+np.kron(Id,Z))+lam*(np.kron(X,Id)-np.kron(Id,X))
    up=np.array([1.,0.,0.,0.]); down=np.array([0.,0.,0.,1.])
    assert abs(up@h@up+g)<1e-14 and abs(down@h@down-g)<1e-14
    # Ferromagnetic bond: singlet/triplet spectral lower bound and one-magnon z=2.
    Y=np.array([[0,-1j],[1j,0]])
    hf=-J/4*(np.kron(X,X)+np.kron(Y,Y)+np.kron(Z,Z))
    assert np.max(np.abs(np.linalg.eigvalsh(hf)-np.array([-J/4]*3+[3*J/4])))<1e-14
    # Four-state embedding, ordered as A-up,A-down,B-up,B-down.
    ss = -hf / J
    ha = -ss + (np.kron(Z,Id)-np.kron(Id,Z))/2
    hb = 16*ss + 8*np.eye(4)
    hp = 12*np.eye(16, dtype=complex)
    aa = np.array([0,1,4,5]); bb = np.array([10,11,14,15])
    hp[np.ix_(aa,aa)] = ha
    hp[np.ix_(bb,bb)] = hb
    pin_eigs = np.linalg.eigvalsh(hp)
    assert abs(pin_eigs[0]+4)<1e-13 and abs(pin_eigs[-1]-12)<1e-13
    # Two adjacent B bonds are positive, the analytic domain bound's input.
    triple = np.kron(hb, Id) + np.kron(Id, hb)
    assert np.linalg.eigvalsh(triple)[0] > -1e-13
    out=dict(status='numerical checks of proved formulas', rows=rows,
             gapped_limit=lam*lam/(2*g),
             ferromagnet_witness_limit=lam*eps*np.pi/2,
             circle_mobius_defect=float(circle_D),
             gapped_raw_halfspread=float(np.ptp(np.linalg.eigvalsh(h))/2),
             gapped_optimized_M=g,
             four_state_bond_extrema=[float(pin_eigs[0]),float(pin_eigs[-1])],
             two_adjacent_B_bonds_min=float(np.linalg.eigvalsh(triple)[0]))
    with open('checks.json','w') as f:
        json.dump(out,f,indent=2)
        f.write('\n')
    print(json.dumps(out,indent=2))

if __name__ == '__main__':
    main()
