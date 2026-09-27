"""Offline independent Gauss-sum oracle for Clifford-valued parity factors.

This uses polynomial elimination, not the tested tensor contraction algorithm.
It discovers dependencies dynamically and is only a validation tool, never a
Clifft hot executor. Arbitrary non-Clifford local phases are rejected.
"""

from __future__ import annotations

import itertools

import numpy as np
from study_wan_scaling import bits


class Quadratic:
    def __init__(self, rank: int):
        self.linear = [0] * rank
        self.edges = [0] * rank

    def edge(self, a: int, b: int) -> None:
        if a == b:
            self.linear[a] = (self.linear[a] + 2) % 4
        else:
            self.edges[a] ^= 1 << b
            self.edges[b] ^= 1 << a

    def parity(self, mask: int, coefficient: int) -> None:
        coefficient %= 4
        indices = bits(mask)
        for q in indices:
            self.linear[q] = (self.linear[q] + coefficient) % 4
        if coefficient & 1:
            for a, b in itertools.combinations(indices, 2):
                self.edge(a, b)

    def product(self, left: int, right: int) -> None:
        for a in bits(left):
            for b in bits(right):
                self.edge(a, b)

    def mean(self) -> complex:
        active = (1 << len(self.linear)) - 1
        result = 1 + 0j
        while active:
            i = bits(active)[0]
            a, neighbors = self.linear[i], self.edges[i] & active
            active ^= 1 << i
            if a & 1:
                result *= (1 + 1j**a) / 2
                self.parity(neighbors, -a)
            elif not neighbors:
                if a == 2:
                    return 0j
            else:
                # Summing this even-linear variable imposes one affine parity
                # constraint; substitute its neighbor instead of enumerating it.
                j = bits(neighbors)[0]
                active ^= 1 << j
                others = neighbors & active
                adjacent = self.edges[j] & active
                c, b = a // 2, self.linear[j]
                result *= 1j ** (b * c) / 2
                self.parity(others, (-1) ** c * b)
                self.parity(adjacent, 2 * c)
                self.product(others, adjacent)
        return result


def mean(masks: list[tuple[int, ...]], local: np.ndarray, rank: int) -> complex:
    polynomial = Quadratic(rank)
    scalar = 1 + 0j
    for terms, values in zip(masks, local, strict=True):
        size = 1 << len(terms)
        ratios = values[:size] / values[0]
        exponents = np.rint(np.angle(ratios) / (np.pi / 2)).astype(int) % 4
        if np.max(np.abs(ratios - (1j) ** exponents)) > 1e-10:
            raise ValueError("Local factor is not Clifford valued")
        scalar *= values[0]
        polynomial.parity(terms[0], int(exponents[1]))
        if len(terms) == 2:
            polynomial.parity(terms[1], int(exponents[2]))
            cross = int(exponents[3] - exponents[1] - exponents[2]) % 4
            if cross not in (0, 2):
                raise ValueError("Nonquadratic pair phase")
            if cross:
                polynomial.product(*terms)
    return scalar * polynomial.mean()


def validate() -> dict:
    rng = np.random.default_rng(270926)
    error = 0.0
    for rank in range(1, 10):
        for _ in range(24):
            q = Quadratic(rank)
            q.linear = list(map(int, rng.integers(4, size=rank)))
            for a, b in itertools.combinations(range(rank), 2):
                if rng.integers(2):
                    q.edge(a, b)
            expected = 0j
            for assignment in range(1 << rank):
                occupied = bits(assignment)
                exponent = sum(q.linear[a] for a in occupied)
                exponent += 2 * sum(
                    bool(q.edges[a] >> b & 1) for a, b in itertools.combinations(occupied, 2)
                )
                expected += 1j**exponent / (1 << rank)
            error = max(error, abs(q.mean() - expected))
    if error > 1e-12:
        raise AssertionError(("Quadratic oracle differs from exhaustive sum", error))
    return dict(cases=216, maximum_absolute_error=error, reference="exhaustive binary sums")
