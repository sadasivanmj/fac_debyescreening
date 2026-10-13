# Test fixtures

`validated_ground_only_energies_cases_A_D.csv`: the 40 ground-configuration total energies (eV)
of Fe XVII - Fe XXVI for the screened cases A - D, from the earlier validated runs
(`research/fe-xvii-xxv-ground-ip/results/run_20261013_patched_rerun`, commit `4cfae69a` on
`research/fe-xxvi-screened-ground`, same FAC sources as this framework including the Fe XXVI
one-electron fix). They were produced on this machine with a `-O0 -g` build WITHOUT the
`-ftrivial-auto-var-init=zero` flag, so the build provenance is not identical; a separate rerun
with the flag reproduced all 40 energies exactly (max |dE| = 0 eV), which is why a tight
tolerance (1e-3 eV) is used.

The unscreened energies and the reference_only model have no earlier validated counterpart.

Independent Fe XXVI check (used in test_smoke_fac.py): the 1s Dirac energy of one electron
in V(r) = -Z exp(-mu r)/r for a point nucleus (Z = 26) shifts by +70.5485 eV (mu = 0.1) and
+175.6218 eV (mu = 0.25) relative to mu = 0, from an independent numerical solver
(`supplement_fe_xxvi_dirac_yukawa.py` on branch `research/fe-xvii-xxv-ground-ip`).
