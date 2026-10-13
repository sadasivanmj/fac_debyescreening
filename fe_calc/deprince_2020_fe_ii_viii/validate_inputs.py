#!/usr/bin/env python3
"""Check the Paper A configuration data (syntax, occupancies, electron counts, duplicates).

Usage: python validate_inputs.py            (exit code 1 if anything is wrong)
Pure Python; FAC is not imported.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE))

import configurations as cfg                                   # noqa: E402
from config_validation import electron_count, validate_configuration_set   # noqa: E402


def validate(verbose: bool = False) -> list[str]:
    problems: list[str] = []
    for ion, (electrons, configs, source) in cfg.IONS.items():
        problems += validate_configuration_set(ion, electrons, configs)
        if electrons != cfg.Z_IRON - cfg.charge_state(ion):
            problems.append(f"{ion}: electron count / charge state inconsistent")
        if verbose:
            print(f"{ion}: {electrons} electrons, {len(configs)} configurations ({source})")
            for c in configs:
                print("   ", c, f"[{electron_count(c)} e]")
    # every Z_eff = Z - N + 1 must go up by one from ion to ion
    counts = [cfg.IONS[i][0] for i in cfg.ION_ORDER]
    if counts != list(range(25, 17, -1)):
        problems.append(f"electron counts {counts} are not 25 .. 18")
    if set(cfg.CASES.values()) != {0.0, 0.1, 0.25}:
        problems.append("screening cases must be mu = 0, 0.1, 0.25")
    return problems


def main() -> int:
    problems = validate(verbose=True)
    for p in problems:
        print("PROBLEM:", p)
    print("configuration check:", "FAILED" if problems else "passed")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
