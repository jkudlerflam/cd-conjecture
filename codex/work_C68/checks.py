"""Small deterministic C68 algebra checks. Run only from this directory."""
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def actions(S, D, B):
    """Distinct exchange eigenvalues make this an exhaustive sign test."""
    eig, basis = np.linalg.eigh(S)
    assert min(np.diff(eig)) > 1e-8
    unitary, antiunitary = [], []
    for signs in itertools.product((-1, 1), repeat=3):
        O = basis @ np.diag(signs) @ basis.T
        det = round(np.linalg.det(O))
        if np.linalg.norm(det * O @ D - D) < 1e-10 and np.linalg.norm(O @ B - B) < 1e-10:
            (unitary if det == 1 else antiunitary).append(O.tolist())
    return dict(unitary=unitary, antiunitary=antiunitary)


def main():
    eps = 0.01
    S = np.diag([1 + eps, 1 - eps, .4])
    examples = {
        'T': (np.array([eps, 2*eps, .5]), np.zeros(3)),
        'K': (np.array([0., 0., .5]), np.array([eps, 2*eps, 0.])),
        'N': (np.array([eps, 2*eps, .5]), np.array([3*eps, 4*eps, 5*eps])),
    }
    symmetries = {name: actions(S, D, B) for name, (D, B) in examples.items()}
    assert [len(symmetries[n]['unitary']) for n in ('T', 'K', 'N')] == [1, 1, 1]
    assert [len(symmetries[n]['antiunitary']) for n in ('T', 'K', 'N')] == [1, 1, 0]
    assert np.allclose(symmetries['T']['antiunitary'][0], -np.eye(3))
    assert np.allclose(symmetries['K']['antiunitary'][0], np.diag([1, 1, -1]))

    r = np.sqrt(1.25)
    delta = .4/r
    gamma = np.arccos(delta)
    stiffness = np.pi/(2*(np.pi-gamma))
    velocity = 2*np.pi*r*np.sin(gamma)/gamma
    alpha = np.arctan(.5)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.diag([1., -1.])
    bond = np.kron(X, X)+np.kron(Y, Y)+.4*np.kron(Z, Z)+.5*(np.kron(X, Y)-np.kron(Y, X))
    U = np.kron(np.eye(2), np.diag(np.exp(.5j*alpha*np.array([1., -1.]))))
    target = r*(np.kron(X, X)+np.kron(Y, Y))+.4*np.kron(Z, Z)
    spiral_error = float(np.max(np.abs(U @ bond @ U.conj().T-target)))
    assert spiral_error < 1e-14
    assert np.isclose(np.pi/(2*(np.pi-np.arccos(0.))), 1.)
    assert np.isclose(2*np.pi*r/np.arccos(0.), 4*r)
    vertices = []
    # x=m^2/(4K)+n^2 K. This range exceeds the derived bounds for x<=2.
    for m in range(-4, 5):
        for n in range(-4, 5):
            if (m, n) == (0, 0) or m < 0 or (m == 0 and n < 0):
                continue
            dim = m*m/(4*stiffness)+n*n*stiffness
            if dim <= 2+1e-12:
                Q = (m*(np.pi+alpha)+n*np.pi+np.pi) % (2*np.pi)-np.pi
                vertices.append(dict(m=m, n=n, dimension=float(dim), momentum=float(Q)))
    assert len(vertices) == 5
    assert min(abs(row['momentum']) for row in vertices) > .4
    assert 4*stiffness > 2
    result = dict(
        scope='Algebra and continuum input checks only. No thermodynamic spectrum was computed.',
        epsilon=eps, symmetries=symmetries,
        base=dict(r=float(r), delta=float(delta), K=float(stiffness), v=float(velocity), alpha=float(alpha)),
        relevant_vertices_up_to_conjugation=vertices,
        neutral_umklapp_dimension=float(4*stiffness),
        spiral_bond_error=spiral_error,
        normalization_checks='At Delta=0 the formulas give K=1 and v=4r.',
    )
    # Single requested parameter record, not a bulk data read.
    source = HERE / '../../contrib/d2_followup_220_1282/candidates_params.json'
    point = json.loads(source.read_text())['1282']
    J = np.array(point['J'])
    B = np.array(point['field'])
    S1282 = (J+J.T)/2
    D1282 = np.array([J[1,2]-J[2,1], J[2,0]-J[0,2], J[0,1]-J[1,0]])/2
    result['point_1282'] = dict(J=J.tolist(), B=B.tolist(), symmetry=actions(S1282,D1282,B))
    assert len(result['point_1282']['symmetry']['unitary']) == 1
    assert len(result['point_1282']['symmetry']['antiunitary']) == 0
    (HERE/'checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
