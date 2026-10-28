"""Full-space local permutation gates for a reversible one-head partial TM.

Stationary moves are replaced by a right/left excursion. Collision and
movement use fresh/updated head sectors, giving disjoint transpositions with
identity defaults. No tape value or position is decoded by a gate.
"""
from dataclasses import dataclass


class UnitMoves:
    def __init__(self, source):
        self.source = source

    def local_rule(self, control, symbol):
        kind, q = control
        if kind == "BOUNCE":
            if self.source.incoming_direction(q) == 0:
                return ("BASE", q), symbol, -1
            return None
        assert kind == "BASE"
        out = self.source.local_rule(q, symbol)
        if out is None:
            return None
        r, b, d = out
        return (("BOUNCE", r), b, 1) if d == 0 else (("BASE", r), b, d)

    def incoming_direction(self, control):
        kind, q = control
        d = self.source.incoming_direction(q)
        if kind == "BOUNCE":
            return 1 if d == 0 else None
        assert kind == "BASE"
        return -1 if d == 0 else d

    def inverse_rule(self, target, written):
        kind, q = target
        d = self.source.incoming_direction(q)
        if kind == "BASE" and d == 0:
            return ("BOUNCE", q), written, -1
        if kind == "BOUNCE" and d != 0:
            return None
        old = self.source.inverse_rule(q, written)
        if old is None:
            return None
        p, a, move = old
        if kind == "BOUNCE":
            assert move == 0
            return ("BASE", p), a, 1
        assert kind == "BASE" and move != 0
        return ("BASE", p), a, move


@dataclass(frozen=True)
class Site:
    symbol: object
    control: object = None
    phase: int = 0

    def __post_init__(self):
        assert self.phase in (0, 1)
        assert self.control is not None or self.phase == 0


class Gates:
    def __init__(self, unit_moves):
        self.machine = unit_moves

    def collision(self, site):
        if site.control is None:
            return site
        rule = (self.machine.local_rule if site.phase == 0 else self.machine.inverse_rule)
        out = rule(site.control, site.symbol)
        if out is None:
            return site
        q, a, _ = out
        return Site(a, q, 1 - site.phase)

    def movement(self, left, right):
        """A two-site permutation, including zero/two-head input sectors."""
        if (left.control is None) == (right.control is None):
            return left, right
        on_left = left.control is not None
        head = left if on_left else right
        d = self.machine.incoming_direction(head.control)
        if d not in (-1, 1):
            return left, right
        forward = head.phase == 1 and (on_left == (d == 1))
        backward = head.phase == 0 and (on_left == (d == -1))
        if not (forward or backward):
            return left, right
        if on_left:
            return Site(left.symbol), Site(right.symbol, head.control, 1-head.phase)
        return Site(left.symbol, head.control, 1-head.phase), Site(right.symbol)

    def forward_sweep(self, sites):
        """Fixed gate word C on all sites, then G on all bonds left to right.

        On one fresh head with a defined source transition this performs
        exactly one normalized TM step. No source-state-dependent scheduling.
        Undefined transitions are handled by the gates' permutation defaults;
        they must not be interpreted as an absorbing reversible halt.
        """
        out = [self.collision(s) for s in sites]
        for i in range(len(out)-1):
            out[i], out[i+1] = self.movement(out[i], out[i+1])
        return out

    def backward_sweep(self, sites):
        out = list(sites)
        for i in range(len(out)-2, -1, -1):
            out[i], out[i+1] = self.movement(out[i], out[i+1])
        return [self.collision(s) for s in out]
