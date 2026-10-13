from pathlib import Path
import subprocess
import sys

from pfac import fac

# ---------------------------------------------------------
# Free-ion ionization potential: Fe II -> Fe III
# No plasma-screening options are set.
# PrintTable(..., 1) writes energies in eV.
# ---------------------------------------------------------

OUT = Path(__file__).resolve().parent

CONFIGS = {
    "FeII": [
        "1s2 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d7 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d6 4s2",
        "1s2 2s2 2p5 3s2 3p6 3d7 4s1",
        "1s2 2s2 2p5 3s2 3p6 3d6 4s2",
        "1s1 2s2 2p6 3s2 3p6 3d7 4s1",
        "1s1 2s2 2p6 3s2 3p6 3d6 4s2",
    ],
    "FeIII": [
        "1s2 2s2 2p6 3s2 3p6 3d6",
        "1s2 2s2 2p6 3s2 3p5 3d7",
        "1s2 2s2 2p6 3s2 3p5 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d5 4s2",
        "1s2 2s2 2p5 3s2 3p6 3d7",
        "1s2 2s2 2p5 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p5 3s2 3p6 3d5 4s2",
        "1s1 2s2 2p6 3s2 3p6 3d7",
        "1s1 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s1 2s2 2p6 3s2 3p6 3d5 4s2",
    ],
}


def calculate_ion(ion):
    """Run one ion calculation in a fresh Python process."""
    group = f"{ion}_ref"
    lev_bin = OUT / f"{ion}_free.lev.b"
    lev_txt = OUT / f"{ion}_free.lev"

    fac.SetAtom("Fe")

    for cfg in CONFIGS[ion]:
        fac.Config(cfg, group=group)

    fac.ConfigEnergy(0)
    fac.OptimizeRadial([group])
    fac.ConfigEnergy(1)

    fac.Structure(str(lev_bin), [group])
    fac.MemENTable(str(lev_bin))

    # Verbose format: energies in eV.
    fac.PrintTable(str(lev_bin), str(lev_txt), 1)

    return lev_txt


def read_e0_eV(filename):
    """Read FAC's E0 reference energy from a verbose level table."""
    with open(filename) as f:
        for line in f:
            if line.strip().startswith("E0"):
                return float(line.split(",")[1])

    raise RuntimeError(f"E0 not found in {filename}")


def main():
    # Each ion gets its own process to avoid carrying FAC state
    # from the first calculation into the second.
    script = str(Path(__file__).resolve())

    for ion in ("FeII", "FeIII"):
        print(f"\nCalculating free-ion {ion}...")
        result = subprocess.run(
            [sys.executable, script, "--ion", ion],
            capture_output=True,
            text=True,
        )

        if result.stdout:
            print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")

        if result.returncode != 0:
            raise RuntimeError(f"{ion} FAC calculation failed.")

    e_feii = read_e0_eV(OUT / "FeII_free.lev")
    e_feiii = read_e0_eV(OUT / "FeIII_free.lev")

    # Fe II (+1) -> Fe III (+2) + e-
    ip_eV = e_feiii - e_feii

    print("\n========== FREE-ION RESULTS ==========")
    print(f"E0(Fe II)  = {e_feii:.8f} eV")
    print(f"E0(Fe III) = {e_feiii:.8f} eV")
    print("--------------------------------------")
    print(f"IP(Fe II)  = {ip_eV:.6f} eV")
    print("======================================")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--ion":
        ion = sys.argv[2]
        if ion not in CONFIGS:
            raise ValueError(f"Unknown ion: {ion}")

        # This process calculates exactly one ion.
        # No plasma screening options are configured.
        calculate_ion(ion)
    else:
        main()


