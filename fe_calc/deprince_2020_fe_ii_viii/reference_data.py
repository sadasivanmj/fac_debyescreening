"""Published comparison values for Paper A (comparison only; never used in a calculation)."""

# Deprince et al., A&A 643, A57 (2020), Table 2: computed IP of the named ion (eV) at mu = 0, 0.1, 0.25.
PAPER_TABLE2 = {
    "Fe II":   {"mu0p00": 18.89,  "mu0p10": 14.20,  "mu0p25": 7.78},
    "Fe III":  {"mu0p00": 32.93,  "mu0p10": 25.60,  "mu0p25": 16.42},
    "Fe IV":   {"mu0p00": 55.37,  "mu0p10": 44.65,  "mu0p25": 29.64},
    "Fe V":    {"mu0p00": 74.92,  "mu0p10": 61.56,  "mu0p25": 42.83},
    "Fe VI":   {"mu0p00": 98.20,  "mu0p10": 82.20,  "mu0p25": 59.79},
    "Fe VII":  {"mu0p00": 123.60, "mu0p10": 104.96, "mu0p25": 78.86},
    "Fe VIII": {"mu0p00": 151.14, "mu0p10": 129.85, "mu0p25": 99.98},
}

# NIST isolated-ion IPs (eV) as listed in the same table (Kramida et al. 2019).
NIST_IP = {"Fe II": 16.19920, "Fe III": 30.651, "Fe IV": 54.91, "Fe V": 75.0,
           "Fe VI": 98.985, "Fe VII": 124.98, "Fe VIII": 151.060}
