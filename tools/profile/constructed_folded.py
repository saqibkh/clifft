"""Constructed terminal folded checks, not an author-supplied MSC protocol.

The code geometry and ideal folded operator extend to odd distances. The cat
schedule is a simple research control; no fault-distance claim is made for it.
"""

from __future__ import annotations

from typing import Any


def make_model(surface: Any) -> dict:
    n, distance = len(surface.points), surface.distance
    cat = list(range(n, n + distance))
    verifier = n + distance
    events: list[dict] = []
    channels: list[dict] = []
    records = 0

    def noise(support: list[int], axis: str | None = None) -> None:
        location = len(channels)
        channels.append(
            dict(location=location, support=support, kind="flip" if axis else "depolarizing")
        )
        if axis:
            channels[-1]["axis"] = axis
        events.append(dict(name="NOISE", location=location))

    def gate(name: str, *support: int, block: int | None = None) -> None:
        event: dict = dict(name=name, support=list(support))
        if block is not None:
            event["native_block"] = block
        events.append(event)
        noise(list(support))

    def measure(q: int) -> None:
        nonlocal records
        noise([q], "X")
        events.append(dict(name="MEASURE", support=[q], axis="Z", inverted=False, record=records))
        records += 1

    for block in range(2):
        # Keep reset faults after all resets, including at the entry boundary.
        events.extend(dict(name="R", support=[q]) for q in [*cat, verifier])
        for q in [*cat, verifier]:
            noise([q], "X")
        gate("H", cat[0])
        for a, b in zip(cat, cat[1:], strict=False):
            gate("CX", a, b)
        gate("CX", cat[0], verifier)
        gate("CX", cat[-1], verifier)
        measure(verifier)
        for j, q in enumerate(surface.fold):
            gate("CSX" if j % 2 == 0 else "CS_DAG_X", cat[j % distance], q, block=block)
        for layer, pairs in enumerate(surface.pairs):
            for j, (a, b) in enumerate(pairs):
                gate("CCZ", cat[(layer + j) % distance], a, b, block=block)
        for a, b in reversed(list(zip(cat, cat[1:], strict=False))):
            gate("CX", a, b)
        gate("T_DAG", cat[0])
        gate("H", cat[0])
        for q in cat:
            measure(q)
    for axis in "XZ":
        for support in surface.checks(axis):
            events.append(
                dict(name="MEASURE_PRODUCT", support=list(support), axis=axis, record=records)
            )
            records += 1
    end = 2 * distance - 1
    return dict(
        provenance="constructed terminal-gadget control; no injection, growth or escape",
        distance=distance,
        events=events,
        channels=channels,
        logical_x={str(surface.index[0, y]): "X" for y in range(0, end, 2)},
        logical_z={str(surface.index[x, 0]): "Z" for x in range(0, end, 2)},
    )
