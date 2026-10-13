# Test fixtures

`linux_validated_fe_ii_iv_energies.csv` is a copy of
`research/debye-screening/results/reference-configurations/reference_energies.csv`
(tracked on `integration/debye-screening-complete`): the 18 Fe II - Fe IV total energies
(eV) from the reference configurations, models `ref` and `gref`, mu = 0, 0.1, 0.25,
calculated on Linux with the same three validated C patches (stable `DensityToSZ`,
`radial:hxs_iter0`, Debye `RadialBound`) and the same calculation sequence.

The build provenance is NOT identical to the macOS runs of this framework (different
platform and compiler flags), so the regression test uses a tolerance (1e-3 eV) and
reports the difference in build provenance instead of requiring bit-identity.
