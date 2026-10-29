"""Exact finite checks of C79's geometry and load accounting, not an RG simulation."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial
from pathlib import Path
import json


def add(*vectors):
    out = {}
    for v in vectors:
        for k, a in v.items():
            out[k] = out.get(k, Q(0)) + a
    return {k: a for k, a in out.items() if a}


def scale(v, a):
    return {k: a * x for k, x in v.items() if a * x}


def delta(k):
    return {k: Q(1)}


class Node:
    serial = 0

    def __init__(self, children=None, terminal_count=0, star=False):
        self.number = Node.serial
        Node.serial += 1
        self.children = children or []
        self.edge_vectors = {}
        self.links = {}
        if not self.children:
            self.terminals = tuple((self.number, j) for j in range(terminal_count))
            for j in range(1, terminal_count):
                a, b = self.terminals[j - 1], self.terminals[j]
                e = ('terminal', self.number, j)
                self.edge_vectors[e] = add(delta(a), scale(delta(b), -1))
            self.height = 0
        else:
            self.terminals = tuple(x for ch in self.children for x in ch.terminals)
            self.height = min(ch.height for ch in self.children) - 1
            for ch in self.children:
                self.edge_vectors.update(ch.edge_vectors)
            for j in range(1, len(self.children)):
                i = 0 if star else j - 1
                left, right = self.children[i], self.children[j]
                a = self.convex(left, j + 1)
                b = self.convex(right, j + 3)
                e = ('bbf', self.number, j)
                self.links[(i, j)] = (e, a, b)
                self.edge_vectors[e] = add(a, scale(b, -1))
        self.owner = {x: i for i, ch in enumerate(self.children) for x in ch.terminals}

    @staticmethod
    def convex(node, seed):
        weights = [1 + ((j + 1) * seed) % 7 for j in range(len(node.terminals))]
        total = sum(weights)
        return {x: Q(w, total) for x, w in zip(node.terminals, weights)}

    def path(self, start, end):
        todo = [(start, [])]
        seen = set()
        while todo:
            now, edges = todo.pop()
            if now == end:
                return edges
            if now in seen:
                continue
            seen.add(now)
            for (i, j), (e, a, b) in self.links.items():
                if now == i:
                    todo.append((j, edges + [(i, j, e, a, b, 1)]))
                elif now == j:
                    todo.append((i, edges + [(j, i, e, b, a, -1)]))
        raise AssertionError('disconnected child tree')

    @lru_cache(None)
    def terminal_difference(self, a, b):
        if a == b:
            return {}
        if not self.children:
            ia, ib = self.terminals.index(a), self.terminals.index(b)
            sign = 1 if ia < ib else -1
            return {('terminal', self.number, j): Q(sign)
                    for j in range(min(ia, ib) + 1, max(ia, ib) + 1)}
        ia, ib = self.owner[a], self.owner[b]
        if ia == ib:
            return self.children[ia].terminal_difference(a, b)
        ans = {}
        current = delta(a)
        for i, j, edge, source, target, sign in self.path(ia, ib):
            ans = add(ans, self.children[i].convex_difference(current, source),
                      {edge: Q(sign)})
            current = target
        return add(ans, self.children[ib].convex_difference(current, delta(b)))

    def convex_difference(self, left, right):
        assert sum(left.values()) == sum(right.values()) == 1
        assert min(left.values()) >= 0 and min(right.values()) >= 0
        out = {}
        for a, u in left.items():
            for b, v in right.items():
                out = add(out, scale(self.terminal_difference(a, b), u * v))
        return out

    def reconstruct(self, coefficients):
        return add(*(scale(self.edge_vectors[e], a) for e, a in coefficients.items()))


def leaf(p):
    return Node(terminal_count=p)


counts = dict(geometric_terminal_pairs=0, geometric_convex_pairs=0,
              minimum_load_surplus=0, factorial_absorption=0,
              weighted_ancestor_sums=0, routing_factorizations=0,
              homogeneity_identities=0, engineering_counter_cases=0, branching_cut_records=0)

shapes = [Node([leaf(2), leaf(4)]),
          Node([Node([leaf(3), leaf(2)]), Node([leaf(4), leaf(2)])]),
          Node([leaf(2), leaf(3), leaf(4), leaf(2)], star=True),
          Node([Node([leaf(2), leaf(3), leaf(2)], star=True),
                Node([leaf(4), leaf(2)]), leaf(3)])]
for node in shapes:
    for a in node.terminals:
        for b in node.terminals:
            coefficients = node.terminal_difference(a, b)
            assert all(abs(c) <= 1 for c in coefficients.values())
            assert node.reconstruct(coefficients) == add(delta(a), scale(delta(b), -1))
            counts['geometric_terminal_pairs'] += 1
    for i, j in product(range(1, 5), repeat=2):
        left, right = node.convex(node, i), node.convex(node, j + 2)
        coefficients = node.convex_difference(left, right)
        assert all(abs(c) <= 1 for c in coefficients.values())
        assert node.reconstruct(coefficients) == add(left, scale(right, -1))
        counts['geometric_convex_pairs'] += 1

# Exhaust every multiplicity pattern on a bounded interval. Birth multiplicity <=2.
for multiplicities in product(range(3), repeat=8):
    b = sum(multiplicities)
    cost = sum(i * x for i, x in enumerate(multiplicities))
    q = ((b - 1) ** 2) // 4 if b else 0
    assert cost >= q
    counts['minimum_load_surplus'] += 1

# Squared form of (b!)^2 2^(-q/2) <= 2^(9b), all arithmetic integral.
for b in range(301):
    q = ((b - 1) ** 2) // 4 if b else 0
    assert factorial(b) ** 4 <= 2 ** (18 * b + q)
    counts['factorial_absorption'] += 1


def all_nodes(node):
    return [node] + [v for ch in node.children for v in all_nodes(ch)]


def sources(node):
    if not node.children:
        return [(e, 0) for e in node.edge_vectors]
    own = [(e, node.height) for e in node.edge_vectors if e[0] == 'bbf' and e[1] == node.number]
    return own + [e for ch in node.children for e in sources(ch)]


# Rational rho suffices to check the double-count identity and geometric bound.
for node in shapes:
    for rho in (Q(1, 2), Q(2, 3), Q(3, 4)):
        sums = [sum((rho ** (t - v.height) for e, t in sources(v)), Q(0))
                for v in all_nodes(node)]
        assert sum(sums) <= len(sources(node)) / (1 - rho)
        counts['weighted_ancestor_sums'] += 1

# Enumerate routes for two branching vertices and verify the product formula.
node = shapes[0]
vertices = [node] + node.children
for loads in product(range(3), repeat=len(vertices)):
    slots = [sources(v) for v, r in zip(vertices, loads) for _ in range(r)]
    rho = Q(1, 2)
    explicit = Q(0)
    birth_scales = [v.height for v, r in zip(vertices, loads) for _ in range(r)]
    for assignment in product(*slots):
        d = sum(t - u for (e, t), u in zip(assignment, birth_scales))
        explicit += rho ** d
    factored = Q(1)
    for v, r in zip(vertices, loads):
        factored *= sum((rho ** (t - v.height) for e, t in sources(v)), Q(0)) ** r
    assert explicit == factored
    counts['routing_factorizations'] += 1

# Polynomial homogeneity, including every local exposed-degree transition.
for m in range(3):
    for a in range(9):
        for k in range(9):
            lhs = (a + k >= m)
            rhs = (k >= max(0, m - a))
            assert lhs == rhs
            r = max(0, m - a)
            assert a + r >= m
            assert r == 0 or a + r <= 2
            counts['homogeneity_identities'] += 1

# Engineering degree is retained separately from localization degree.
for p in (2, 4, 6):
    m = {2: 2, 4: 1, 6: 0}[p]
    allowed = ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1))
    for degrees in product(allowed, repeat=p):
        e = sum(x for x, y in degrees)
        a = sum(y for x, y in degrees)
        r = max(0, m - e - a)
        assert a + r >= max(0, m - e)
        assert e + a + r >= m
        assert r == 0 or e + a + r <= 2
        assert Q(p, 2) - 2 + e + max(0, m - e) >= Q(p, 2) - 2 + m
        counts['engineering_counter_cases'] += 1

# A genuinely branched abstract GN field-survival record. The test is about
# field labels and cut accounting, not nonvanishing of its fermionic amplitude.
parent = {'a': 'u', 'b': 'u', 'c': 'v', 'd': 'v', 'u': 'w', 'v': 'w', 'w': 'root'}
base = {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'u': -3, 'v': -5, 'w': -9, 'root': -12}
# One carrier per quartic child; one terminates at each binary merge.
paths = [('a', 'root'), ('b', 'u'), ('c', 'w'), ('d', 'v')]
for stretch in range(1, 101):
    h = {v: s * stretch for v, s in base.items()}
    flow = {v: 0 for v in parent}
    for start, end in paths:
        here = start
        while here != end:
            flow[here] += 1
            here = parent[here]
    assert all(f >= 1 for f in flow.values())
    direct = sum(h[end] - h[start] for start, end in paths)
    cuts = -sum(flow[v] * (h[v] - h[parent[v]]) for v in parent)
    assert direct == cuts
    assert direct <= -sum(h[v] - h[parent[v]] for v in parent)
    counts['branching_cut_records'] += 1

out = dict(status='passed', arithmetic='exact rational and integer', **counts,
           maximum_factorial_load=300,
           scope='Finite checks of convex routing, factorial surplus, route sums, and cut algebra. No interacting RG simulation.')
Path(__file__).with_name('checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
