"""Executable two-site bond assembly; no history controller is assumed implicitly.

CoreData supplies actual finite matrices (or sparse column callables).  The
small fixture has one logical fermion, not the proposed high-channel core.
The formulas define exact real/complex coefficients; evaluation uses float64.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
import math
import numpy as np


def plus(theta):
    return np.array([1, np.exp(1j * theta)]) / math.sqrt(2)


def helix(theta):
    v = np.array([0, 1, -np.exp(-1j * theta), 0]) / math.sqrt(2)
    return np.outer(v, v.conj())


def matrix_column(matrix):
    """Return a sparse column callable with no numerical thresholding."""
    def column(index):
        return {int(i): matrix[i, index] for i in np.flatnonzero(matrix[:, index])}
    return column


@dataclass
class Port:
    left: object
    left_adjoint: object
    right: object
    right_adjoint: object
    alpha_left: float
    alpha_right: float
    fermi_momentum: float


@dataclass
class CoreData:
    dimension: int
    length: int
    bulk: object                 # column index c_left*d+c_right
    pin_left: object
    pin_right: object
    upper_charge: tuple
    species_charges: tuple       # charge vector for every local basis state
    charge_offset_per_cell: int
    ports: tuple
    label: str


def one_fermion_core(k=3):
    """Explicit positive core with kernel span{|vac>, |uniform one particle>}.

    The singlet projector moves particles; the |11> penalty excludes every
    sector with >=2 particles from the zero kernel.  This is a test fixture.
    """
    v = np.array([0, 1, -1, 0]) / math.sqrt(2)
    h = np.outer(v, v)
    h[3, 3] += 1
    raising = np.array([[0, 0], [1, 0]], complex)
    return CoreData(
        dimension=2, length=k, bulk=matrix_column(h),
        pin_left=matrix_column(np.zeros((2, 2))),
        pin_right=matrix_column(np.zeros((2, 2))),
        upper_charge=(0, 1), species_charges=((0,), (1,)),
        charge_offset_per_cell=0,
        ports=(Port(matrix_column(raising), matrix_column(raising.conj().T),
                    matrix_column(raising), matrix_column(raising.conj().T),
                    1 / math.sqrt(k), 1 / math.sqrt(k), 7 * math.pi / 24),),
        label="one-fermion positive-core fixture; not a c-d violation")


def charge_ground(k, holes, U, vacancy_hop):
    n = k + holes
    states = list(combinations(range(n), holes))
    index = {s: i for i, s in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for col, locations in enumerate(states):
        occupied = set(locations)
        for j in range(n - 1):
            if (j in occupied) != (j + 1 in occupied):
                h[col, col] += U
                target = occupied ^ {j, j + 1}
                h[index[tuple(sorted(target))], col] -= vacancy_hop
    energies, vectors = np.linalg.eigh(h)
    chi = vectors[:, 0]
    chi *= 1 if chi.sum() > 0 else -1
    # The unique positive ground vector is exactly reflection even.  Impose
    # this known symmetry to remove roundoff mixing across the tiny doublet.
    reflected = np.array([chi[index[tuple(sorted(n-1-j for j in s))]] for s in states])
    chi = chi + reflected
    chi /= np.linalg.norm(chi)
    assert np.linalg.norm(h @ chi - energies[0]*chi) < 1e-12
    assert np.min(chi) > 0
    return states, chi, energies


class MicroscopicBond:
    """Full repeated bond, including its action in all illegal local sectors.

    Site order is (selector, wrapper, helix_A, helix_B).  Selector indices
    are #, L, X, R. Wrapper indices are H0,H1,H2, then (tag,core,filter).
    A column is evaluated without allocating the full d^2 by d^2 matrix.
    """
    allowed = {(0, 1), (0, 2), (1, 1), (1, 2),
               (2, 0), (2, 3), (3, 0), (3, 3)}
    cuts = {(2, 0), (3, 0)}

    def __init__(self, core, s=1e-7, U=1.0, length_hop=0.5):
        self.core = core
        self.k, self.b = core.length, core.length + 1
        if self.k < 3:
            raise ValueError("The selector parameter formulas require k>=3")
        self.dc, self.U, self.s = core.dimension, U, s
        self.occupied_dimension = 2 * self.dc
        self.wrapper_dimension = 3 + 6 * self.dc
        self.dimension = 16 * self.wrapper_dimension
        self.v = U / (2 * self.b)
        self.length_hop = length_hop
        self.theta = 2 * math.pi / (self.k + 2)
        kappa = 1 / ((self.k + 2) * (2 * self.k + 3))
        mismatch = math.sin(math.pi / (self.k + 2)) ** 2
        epsilon = kappa * mismatch / ((1 + math.sqrt(kappa)) ** 2 + mismatch)
        self.filter_scale = 2 * U / epsilon
        self.hfilter = helix(self.theta)
        self.afilter = np.eye(2) - np.outer(plus(0), plus(0).conj())
        p = plus((self.k - 1) * self.theta)
        self.bfilter = np.eye(2) - np.outer(p, p.conj())
        self.charge_states, self.chis, self.charge_energies, self.r = [], [], [], []
        for holes in (1, 2):
            states, chi, energies = charge_ground(self.k, holes, U, self.v)
            self.charge_states.append(states)
            self.chis.append(chi)
            self.charge_energies.append(energies)
            self.r.append(1 - sum(x*x for loc, x in zip(states, chi) if 0 in loc))
        a = lambda ell: 2 * (1 - math.cos(math.pi / ell))
        b = self.b
        ec = a(b + 1) - a(b)
        muc = 2 - ((b + 1) * a(b) - b * a(b + 1))
        eta = min((a(b - 1) - 2*a(b) + a(b + 1))/2,
                  (a(b) - 2*a(b + 1) + a(b + 2))/3)
        self.selector_hop = 2 * U * (1/eta + b/muc + 1/abs(ec))
        e1, e2 = (e[0] for e in self.charge_energies)
        self.mu0 = self.selector_hop * muc + 4*U - (b + 1)*e1 + b*e2
        self.base_density = self.selector_hop * ec + e2 - e1
        self.Ksel = 4 * self.selector_hop + 9 * U
        self.phi = 1 / b
        self.hhelix, self.pin = helix(self.phi), self.afilter
        u1, u2 = plus(-(b - 1)*self.phi), plus(-b*self.phi)
        self.EL = ((np.eye(2) - np.outer(u2, u2.conj())) / self.r[0]
                   + (np.eye(2) - np.outer(u1, u1.conj())) / self.r[1])
        self.EL /= math.sin(self.phi / 2) ** 2
        self.ER = self.EL.conj()
        self.A = sum(math.sqrt(2/b) * math.sin(math.pi*i/b)
                     * math.sqrt(2/(b+1)) * math.sin(math.pi*i/(b+1))
                     for i in range(1, b))
        lookup = dict(zip(self.charge_states[1], self.chis[1]))
        self.beta_wrap = sum(x * lookup[(loc[0], self.k + 1)]
                             for loc, x in zip(self.charge_states[0], self.chis[0])) ** 2
        self.beta_helix = math.cos(self.phi / 2) ** (2*b)
        self.marker_amplitude = length_hop / (self.A**2 * self.beta_wrap * self.beta_helix)
        self.payload_energy = sum((p.fermi_momentum - math.tan(p.fermi_momentum)) / math.pi
                                  for p in core.ports)
        self.zeta = 2 * length_hop / math.pi - self.payload_energy

    def encode(self, selector, wrapper, helix_a, helix_b):
        return ((selector * self.wrapper_dimension + wrapper) * 2 + helix_a) * 2 + helix_b

    def decode(self, i):
        if not 0 <= i < self.dimension:
            raise IndexError(i)
        i, b = divmod(i, 2)
        i, a = divmod(i, 2)
        selector, wrapper = divmod(i, self.wrapper_dimension)
        return selector, wrapper, a, b

    def occupied(self, wrapper):
        if wrapper < 3:
            return None
        tag, record = divmod(wrapper - 3, self.occupied_dimension)
        core, filt = divmod(record, 2)
        return tag, core, filt

    def wrap_encode(self, tag, core, filt):
        return 3 + tag * self.occupied_dimension + 2 * core + filt

    def arrows(self, w):
        if w < 3:
            return w, w
        return ((0, 1), (1, 1), (1, 2))[self.occupied(w)[0]]

    def site_charge(self, i):
        _, w, _, _ = self.decode(i)
        occ = self.occupied(w)
        return self.core.species_charges[occ[1]] if occ else (0,) * len(self.core.species_charges[0])

    def z(self, selector, wrapper):
        occ = self.occupied(wrapper)
        charge = self.core.upper_charge[occ[1]] if occ else 0
        return charge + (self.zeta - self.core.charge_offset_per_cell) * (selector == 0)

    @lru_cache(maxsize=None)
    def wrapper_column(self, x, y):
        out = {}
        def add(a, b, c):
            if c != 0:
                out[a, b] = out.get((a, b), 0) + c
        ox, oy = self.occupied(x), self.occupied(y)
        diagonal = 6*self.U * (self.arrows(x)[1] != self.arrows(y)[0])
        diagonal += self.U * ((x < 3) != (y < 3))
        add(x, y, diagonal)
        for hole_left, tag, hole_right in ((0, 0, 1), (1, 1, 1), (1, 2, 2)):
            if x == hole_left and oy and oy[0] == tag:
                add(y, hole_right, -self.v)
            if y == hole_right and ox and ox[0] == tag:
                add(hole_left, x, -self.v)
        if ox and oy and (ox[0], oy[0]) in ((0, 1), (0, 2), (1, 1), (1, 2)):
            tx, cx, fx = ox
            ty, cy, fy = oy
            for row, value in self.core.bulk(cx*self.dc + cy).items():
                dx, dy = divmod(row, self.dc)
                add(self.wrap_encode(tx, dx, fx), self.wrap_encode(ty, dy, fy), value)
            for row in np.flatnonzero(self.hfilter[:, 2*fx + fy]):
                gx, gy = divmod(int(row), 2)
                add(self.wrap_encode(tx, cx, gx), self.wrap_encode(ty, cy, gy),
                    self.filter_scale*self.hfilter[row, 2*fx + fy])
        if ox and ox[0] == 0:
            tag, c, f = ox
            add(x, y, -2*self.U)
            for target, value in self.core.pin_left(c).items():
                add(self.wrap_encode(tag, target, f), y, value)
            for target in range(2):
                add(self.wrap_encode(tag, c, target), y, self.filter_scale*self.afilter[target, f])
        if oy and oy[0] == 2:
            tag, c, f = oy
            add(x, y, -2*self.U)
            for target, value in self.core.pin_right(c).items():
                add(x, self.wrap_encode(tag, target, f), value)
            for target in range(2):
                add(x, self.wrap_encode(tag, c, target), self.filter_scale*self.bfilter[target, f])
        return {pair: value for pair, value in out.items() if value != 0}

    def column(self, pair_index, component="total"):
        """Return {row: coefficient}; components: base, hopping, move, field.

        Perturbation components include the supplied s.  Hermitian conjugate
        transitions are explicitly generated, not inferred from a valid sector.
        """
        i, j = divmod(pair_index, self.dimension)
        x, y = self.decode(i), self.decode(j)
        sx, wx, ax, bx = x
        sy, wy, ay, by = y
        selector = sx, sy
        out = {}
        def add(left, right, value):
            if value != 0:
                row = self.encode(*left)*self.dimension + self.encode(*right)
                out[row] = out.get(row, 0) + value
        def wants(name):
            return component in ("total", name)
        if wants("base"):
            add(x, y, self.Ksel * (selector not in self.allowed)
                + self.mu0/2 * ((sx == 0) + (sy == 0)))
            if selector == (1, 2):
                add((2, wx, ax, bx), (3, wy, ay, by), -self.selector_hop)
            elif selector == (2, 3):
                add((1, wx, ax, bx), (2, wy, ay, by), -self.selector_hop)
            if selector in self.allowed - self.cuts:
                for (rx, ry), value in self.wrapper_column(wx, wy).items():
                    add((sx, rx, ax, bx), (sy, ry, ay, by), value)
                for row in range(4):
                    ra, rb = divmod(row, 2)
                    add((sx, wx, ra, bx), (sy, wy, rb, by), self.hhelix[row, 2*ax + ay])
                    add((sx, wx, ax, ra), (sy, wy, ay, rb), self.hhelix[row, 2*bx + by])
            elif selector in self.cuts:
                for r in range(2):
                    add((sx, wx, ax, r), y, self.pin[r, bx])
                    add(x, (sy, wy, r, by), self.pin[r, ay])
        if wants("field"):
            add(x, y, self.s/2*(self.z(sx, wx) + self.z(sy, wy)))
        if wants("move"):
            # The marker moves right, transferring the left physical site.
            if selector == (0, 1) and wx == 0:
                coeff = np.exp(1j*self.b*self.phi*(ax + bx - 1))
                add((3, 2, ax, bx), (0, wy, ay, by), -self.s*self.marker_amplitude*coeff)
            elif selector == (3, 0) and wx == 2:
                coeff = np.exp(-1j*self.b*self.phi*(ax + bx - 1))
                add((0, 0, ax, bx), (1, wy, ay, by), -self.s*self.marker_amplitude*coeff)
        ox, oy = self.occupied(wx), self.occupied(wy)
        if wants("hopping") and selector in self.cuts and ox and oy and (ox[0], oy[0]) == (2, 0):
            tx, cx, fx = ox
            ty, cy, fy = oy
            for port in self.core.ports:
                scale = -self.s / (2*math.cos(port.fermi_momentum)*port.alpha_left*port.alpha_right)
                for first, second in ((port.right, port.left_adjoint), (port.right_adjoint, port.left)):
                    for rx, vx in first(cx).items():
                        for ry, vy in second(cy).items():
                            for a in range(2):
                                for b in range(2):
                                    add((sx, self.wrap_encode(tx, rx, fx), a, bx),
                                        (sy, self.wrap_encode(ty, ry, fy), ay, b),
                                        scale*vx*vy*self.ER[a, ax]*self.EL[b, by])
        return {row: value for row, value in out.items() if value != 0}

    def manifest(self):
        return dict(core=self.core.label, core_dimension=self.dc, core_length=self.k,
                    local_dimension=self.dimension, bond_shape=[self.dimension**2]*2,
                    selected_lengths=[self.b, self.b+1], perturbation_strength=self.s,
                    filter_scale=self.filter_scale, selector_hop=self.selector_hop,
                    selector_penalty=self.Ksel, marker_chemical_potential=self.mu0,
                    base_energy_density=self.base_density, endpoint_occupancies=self.r,
                    selector_overlap=self.A, vacancy_overlap=self.beta_wrap,
                    helix_overlap=self.beta_helix, marker_amplitude=self.marker_amplitude,
                    logical_chemical_potential_per_cell=self.zeta,
                    status="Explicit assembly fixture; no thermodynamic or c-d violation claim")
