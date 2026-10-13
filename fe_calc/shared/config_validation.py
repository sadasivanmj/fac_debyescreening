"""Generic checks for FAC configuration strings (standard library only).

A configuration is written as FAC expects it: shells separated by spaces, each
``<n><l><occupation>``, e.g. ``"1s2 2s2 2p5 3d1"``.
"""

from __future__ import annotations

import re

SHELL = re.compile(r"^([1-9])([spdfgh])(\d+)$")
L_VALUE = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5}


def parse_configuration(config: str) -> dict[str, int]:
    """Return {shell: occupation}; raise ValueError for bad syntax or occupancy."""
    shells: dict[str, int] = {}
    for token in config.split():
        m = SHELL.match(token)
        if not m:
            raise ValueError(f"bad shell token {token!r} in {config!r}")
        n, l, occ = int(m.group(1)), m.group(2), int(m.group(3))
        if L_VALUE[l] >= n:
            raise ValueError(f"shell {n}{l} does not exist ({config!r})")
        capacity = 2 * (2 * L_VALUE[l] + 1)
        if not 1 <= occ <= capacity:
            raise ValueError(f"occupation {occ} of {n}{l} outside 1..{capacity} ({config!r})")
        name = f"{n}{l}"
        if name in shells:
            raise ValueError(f"shell {name} repeated in {config!r}")
        shells[name] = occ
    if not shells:
        raise ValueError("empty configuration")
    return shells


def electron_count(config: str) -> int:
    return sum(parse_configuration(config).values())


def canonical(config: str) -> tuple:
    """Order-independent form, used to find duplicates."""
    return tuple(sorted(parse_configuration(config).items()))


def validate_configuration_set(ion: str, electrons: int, configs: list[str]) -> list[str]:
    """Return a list of problems (empty = valid): syntax, occupancies, electron
    count of every configuration, and duplicates."""
    problems: list[str] = []
    seen: dict[tuple, str] = {}
    for c in configs:
        try:
            n = electron_count(c)
        except ValueError as e:
            problems.append(f"{ion}: {e}")
            continue
        if n != electrons:
            problems.append(f"{ion}: {c!r} has {n} electrons, expected {electrons}")
        key = canonical(c)
        if key in seen:
            problems.append(f"{ion}: duplicate configuration {c!r} (same as {seen[key]!r})")
        seen[key] = c
    if not configs:
        problems.append(f"{ion}: no configurations")
    return problems
