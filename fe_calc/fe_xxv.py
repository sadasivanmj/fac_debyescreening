#!/usr/bin/env python3
"""FAC calculation of the Fe XXV -> Fe XXVI ionization potential."""

from pathlib import Path
from datetime import datetime
import csv
import math
import subprocess
import sys

Z = 26.0
MU_VALUES = [0.0, 0.1, 0.25]
NPS = 1.0e-6                 # plasma density, Angstrom^-3
RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018
UPS = 0.0
VXF = 1

# Fe XXV: reference configurations plus SD excitations into n=3.
FE25_CONFIGS = [
    "1s2", "1s1 2s1", "1s1 2p1",
    "1s1 3s1", "1s1 3p1", "1s1 3d1",
    "2s1 3s1", "2s1 3p1", "2s1 3d1",
    "2p1 3s1", "2p1 3p1", "2p1 3d1",
    "3s2", "3s1 3p1", "3s1 3d1",
    "3p2", "3p1 3d1", "3d2",
]

# Fe XXVI is hydrogen-like; only its 1s ground configuration is needed.
FE26_CONFIGS = ["1s1"]
CONFIGS = {"FeXXV": FE25_CONFIGS, "FeXXVI": FE26_CONFIGS}
ION_CHARGE = {"FeXXV": 24.0, "FeXXVI": 25.0}


def tag(mu):
    return "mu0" if mu == 0.0 else "mu" + format(mu, ".6g").replace(".", "p")


def temperature_eV(mu):
    """Temperature that gives Debye length lambda_D = 1/mu (in a0)."""
    lambda_d = 1.0 / mu
    density_au = NPS * RBOHR_ANG**3
    return lambda_d**2 * 4.0 * math.pi * density_au * (1.0 + UPS) * HARTREE_EV


def read_e0_eV(filename):
    with open(filename) as f:
        for line in f:
            if line.strip().startswith("E0"):
                return float(line.split(",")[1])
    raise RuntimeError(f"Could not find E0 in {filename}")


def worker(ion, mu, outdir):
    # Import PFAC only in the worker process: FAC keeps global native state.
    from pfac import fac

    outdir = Path(outdir).resolve()
    group = "g1"
    level_bin = outdir / f"{ion}_{tag(mu)}.lev.b"
    level_txt = outdir / f"{ion}_{tag(mu)}.lev"

    print(f"Worker: {ion}, mu={mu:g}", flush=True)
    print(f"Output: {level_txt}", flush=True)
    fac.SetAtom("Fe")

    if mu > 0.0:
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)

    print("Adding configurations...", flush=True)
    for cfg in CONFIGS[ion]:
        fac.Config(cfg, group=group)

    if mu > 0.0:
        lambda_d = 1.0 / mu
        tps = temperature_eV(mu)
        print(
            f"Screening: lambda_D={lambda_d:g} a0, "
            f"nps={NPS:g} Angstrom^-3, T={tps:.8g} eV",
            flush=True,
        )
        # PFAC signature: PlasmaScreen(zps, nps, tps, m, ups, vxf).
        fac.PlasmaScreen(ION_CHARGE[ion], NPS, tps, 1, UPS, VXF)
    else:
        print("Unscreened calculation; PlasmaScreen is not called.", flush=True)

    print("ConfigEnergy(0)...", flush=True)
    fac.ConfigEnergy(0)
    print("OptimizeRadial...", flush=True)
    fac.OptimizeRadial([group])
    print("ConfigEnergy(1)...", flush=True)
    fac.ConfigEnergy(1)
    print("Structure...", flush=True)
    fac.Structure(str(level_bin), [group])
    fac.MemENTable(str(level_bin))
    fac.PrintTable(str(level_bin), str(level_txt), 1)

    e0 = read_e0_eV(level_txt)
    print(f"RESULT {ion} mu={mu:g} E0={e0:.10f} eV", flush=True)


def tail(path, n=50):
    try:
        lines = path.read_text(errors="replace").splitlines()
        return "\n".join(lines[-n:])
    except OSError:
        return "(worker log could not be read)"


def main(mu_values):
    script = Path(__file__).resolve()
    outdir = script.parent / ("fe25_fac_output_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
    outdir.mkdir()  # Never overwrite an existing output directory.
    print(f"Output directory: {outdir}")
    print("Protocol: FAC ConfigEnergy(0), OptimizeRadial, ConfigEnergy(1), Structure.")
    print("Fe XXV uses 18 explicit reference/SD configurations; Fe XXVI uses 1s.")

    energies = {}
    for mu in mu_values:
        for ion in ("FeXXV", "FeXXVI"):
            worker_log = outdir / f"{ion}_{tag(mu)}.log"
            command = [sys.executable, str(script), "--worker", ion, str(mu), str(outdir)]
            print(f"\nRunning {ion}, mu={mu:g}; log: {worker_log}")
            with worker_log.open("w") as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)

            if result.returncode != 0:
                print(f"\nERROR: {ion}, mu={mu:g} failed (exit {result.returncode}).")
                print(f"Full worker log: {worker_log}")
                print("Last lines from worker log:\n" + tail(worker_log))
                raise SystemExit(result.returncode)

            table = outdir / f"{ion}_{tag(mu)}.lev"
            energies[(ion, mu)] = read_e0_eV(table)

    with (outdir / "energies.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ion", "mu_a0_inverse", "E0_eV"])
        for mu in mu_values:
            for ion in ("FeXXV", "FeXXVI"):
                writer.writerow([ion, mu, f"{energies[(ion, mu)]:.10f}"])

    ip0 = energies[("FeXXVI", mu_values[0])] - energies[("FeXXV", mu_values[0])]
    with (outdir / "ip.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["mu_a0_inverse", "E0_FeXXV_eV", "E0_FeXXVI_eV", "IP_eV", "delta_IP_eV"])
        for mu in mu_values:
            e25 = energies[("FeXXV", mu)]
            e26 = energies[("FeXXVI", mu)]
            ip = e26 - e25
            writer.writerow([mu, f"{e25:.10f}", f"{e26:.10f}", f"{ip:.10f}", f"{ip-ip0:.10f}"])
            print(f"mu={mu:g}: IP={ip:.6f} eV, delta IP={ip-ip0:+.6f} eV")

    print(f"\nResults saved in: {outdir}")


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "--worker":
        worker(sys.argv[2], float(sys.argv[3]), sys.argv[4])
    else:
        requested = [float(x) for x in sys.argv[1:]]
        if any(mu < 0.0 for mu in requested):
            raise SystemExit("mu must be non-negative")
        # Always include the unscreened baseline.
        mu_values = sorted(set([0.0] + requested)) if requested else MU_VALUES
        main(mu_values)

