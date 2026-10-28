"""Finite cyclic endpoint transport with readout safe within every sweep.

This module assumes two verified interior output-location marks and one
head at OUT_L. It does not prepare their quantum records or initialize the
marks. Only the quantum record is swapped; neutral marks and READY remain
at their physical sites. There are no host-side phase handoffs.
"""
from dataclasses import dataclass
import json
from pathlib import Path


MARKS = ("LEFT", "BULK", "RIGHT", "OUT_L", "OUT_R")
CONTROLS = (
    ("LEFT_OUT", -1, "SWAP"),
    ("FIND_RIGHT", 1, "I"),
    ("RIGHT_OUT", 1, "SWAP"),
    ("READ_BACK", -1, "I"),
    ("LEFT_RETURN", 1, "SWAP"),
    ("FIND_END", 1, "I"),
    ("RIGHT_RETURN", -1, "SWAP"),
    ("FIND_START", -1, "I"),
)


@dataclass(frozen=True)
class Site:
    mark: str
    ready: int = 0
    record: int = 0
    control: str | None = None
    phase: int = 0


class Readout:
    def __init__(self):
        self.symbols = tuple((mark, ready) for mark in MARKS for ready in range(2))
        self.controls = tuple(q for q, _, _ in CONTROLS)
        self.direction = {q: d for q, d, _ in CONTROLS}
        self.oracle = {q: u for q, _, u in CONTROLS}
        self.B, self.K = len(self.symbols), len(self.controls)
        self.partial, self.images = {}, {}
        self._build()
        self.table = list(range(self.B * self.K))
        for a, b in self.partial.items():
            self.table[a] = b
        # Static marks are preserved even on malformed inputs. READY need
        # not be conserved, so complete within each static-mark sector.
        for mark in MARKS:
            sector = {a for a in range(self.B * self.K) if self.decode(a)[1][0] == mark}
            domain, image = sector - self.partial.keys(), sector - self.images.keys()
            common = domain & image
            for a, b in zip(sorted(domain-common), sorted(image-common)):
                self.table[a] = b
        assert sorted(self.table) == list(range(self.B * self.K))
        self.inverse = [0] * len(self.table)
        for a, b in enumerate(self.table):
            self.inverse[b] = a

    def encode(self, q, symbol):
        return self.controls.index(q) * self.B + self.symbols.index(symbol)

    def decode(self, a):
        q, s = divmod(a, self.B)
        return self.controls[q], self.symbols[s]

    def add(self, q, mark, ready, target, target_ready=None):
        if target_ready is None:
            target_ready = ready
        a, b = self.encode(q, (mark, ready)), self.encode(target, (mark, target_ready))
        assert a not in self.partial, (q, mark, ready)
        assert b not in self.images, (target, mark, target_ready)
        self.partial[a], self.images[b] = b, a

    def loop(self, q, marks):
        for mark in marks:
            self.add(q, mark, 0, q)

    def _build(self):
        self.loop("LEFT_OUT", ("BULK",))
        self.add("LEFT_OUT", "LEFT", 0, "FIND_RIGHT", 1)
        self.loop("FIND_RIGHT", ("BULK", "OUT_L"))
        self.add("FIND_RIGHT", "OUT_R", 0, "RIGHT_OUT")
        self.loop("RIGHT_OUT", ("BULK",))
        self.add("RIGHT_OUT", "RIGHT", 0, "READ_BACK", 1)
        self.loop("READ_BACK", ("BULK", "OUT_L", "OUT_R"))
        self.add("READ_BACK", "LEFT", 1, "LEFT_RETURN", 0)
        self.loop("LEFT_RETURN", ("BULK",))
        self.add("LEFT_RETURN", "OUT_L", 0, "FIND_END")
        self.loop("FIND_END", ("BULK", "OUT_R"))
        self.add("FIND_END", "RIGHT", 1, "RIGHT_RETURN", 0)
        self.loop("RIGHT_RETURN", ("BULK",))
        self.add("RIGHT_RETURN", "OUT_R", 0, "FIND_START")
        self.loop("FIND_START", ("BULK",))
        self.add("FIND_START", "OUT_L", 0, "LEFT_OUT")

    def initialize(self, n, left, right):
        assert 0 < left < right < n-1
        marks = ["LEFT"] + ["BULK"]*(n-2) + ["RIGHT"]
        marks[left], marks[right] = "OUT_L", "OUT_R"
        return tuple(Site(mark, record=j, control="LEFT_OUT" if j == left else None)
                     for j, mark in enumerate(marks))

    def g(self, a, b):
        """Whole-space two-site involution, including zero/two-head states."""
        if (a.control is None) == (b.control is None):
            return a, b
        head, on_left = (a, True) if a.control is not None else (b, False)
        d = self.direction[head.control]
        forward = head.phase == 0 and ((d == 1) == on_left)
        reverse = head.phase == 1 and ((d == -1) == on_left)
        if not (forward or reverse):
            return a, b
        ra, rb = a.record, b.record
        if self.oracle[head.control] == "SWAP":
            ra, rb = rb, ra
        phase = 1 - head.phase
        if on_left:
            return Site(a.mark, a.ready, ra), Site(b.mark, b.ready, rb, head.control, phase)
        return Site(a.mark, a.ready, ra, head.control, phase), Site(b.mark, b.ready, rb)

    def f(self, a):
        if a.control is None:
            return a
        table = self.table if a.phase else self.inverse
        q, (mark, ready) = self.decode(table[self.encode(a.control, (a.mark, a.ready))])
        return Site(mark, ready, a.record, q, 1-a.phase)

    def micro_sweep(self, tape):
        """Yield after every G bond and every F site, in actual sweep order."""
        tape = list(tape)
        for j in range(len(tape)-1):
            tape[j], tape[j+1] = self.g(tape[j], tape[j+1])
            yield tuple(tape), ("G", j)
        for j in reversed(range(len(tape))):
            tape[j] = self.f(tape[j])
            yield tuple(tape), ("F", j)

    @staticmethod
    def cycle_steps(n, left, right):
        return 4*(n-1)

    def manifest(self):
        return dict(scope="Endpoint transport/readout primitive; assumes verified output marks and prepared input",
                    neutral_symbols_B=self.B, controls_K=self.K,
                    head_overlay_factor=(1+2*self.K)*self.B,
                    prescribed_entries=len(self.partial),
                    symbols=self.symbols, controls=CONTROLS,
                    collision_permutation=self.table)


if __name__ == "__main__":
    path = Path(__file__).with_name("transition_table.json")
    path.write_text(json.dumps(Readout().manifest(), indent=2)+"\n")
