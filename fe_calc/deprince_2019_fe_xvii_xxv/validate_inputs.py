#!/usr/bin/env python3
"""Check the Paper B configuration data (syntax, occupancies, electron counts, duplicates).

Usage: python validate_inputs.py            (exit code 1 if anything is wrong)
Pure Python; FAC is not imported.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE))

import configurations as cfg                                   # noqa: E402
from config_validation import canonical, electron_count, validate_configuration_set   # noqa: E402


def validate(verbose: bool = False) -> list[str]:
    problems: list[str] = []
    for ion, (electrons, refs) in cfg.IONS.items():
        problems += validate_configuration_set(ion, electrons, refs)
        if electrons != cfg.Z_IRON - cfg.charge_state(ion):
            problems.append(f"{ion}: electron count / charge state inconsistent")
        for model in cfg.CONFIGURATION_MODELS:
            used = cfg.configurations(ion, model)
            if model == "ground_only" and len(used) != 1:
                problems.append(f"{ion}: ground_only must contain exactly the ground configuration")
            if not set(map(canonical, used)) <= set(map(canonical, refs)):
                problems.append(f"{ion}: {model} contains a configuration that is not a reference configuration")
        if verbose:
            print(f"{ion}: {electrons} electrons, {len(refs)} reference configuration(s)")
            for c in refs:
                print("   ", c, f"[{electron_count(c)} e]")
    counts = [cfg.IONS[i][0] for i in cfg.ION_ORDER]
    if counts != list(range(10, 0, -1)):
        problems.append(f"electron counts {counts} are not 10 .. 1")
    if len(cfg.IONS["Fe XXIV"][1]) != 8:
        problems.append("Fe XXIV (1s+2s+2p)^3 must have 8 configurations")
    for case, (mu_ne, mu_ee) in cfg.CASES.items():
        if mu_ee not in (0.0, mu_ne):
            problems.append(f"case {case}: FAC supports mu_ee = 0 or mu_ee = mu_ne only")
    return problems


def main() -> int:
    problems = validate(verbose=True)
    for p in problems:
        print("PROBLEM:", p)
    print("configuration check:", "FAILED" if problems else "passed")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
