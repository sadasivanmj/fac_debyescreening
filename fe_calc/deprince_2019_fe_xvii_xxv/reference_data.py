"""Published comparison values for Paper B (comparison only; never used in a calculation).

Deprince et al., A&A 626, A83 (2019).  IPs in eV of the named ion.
"""

# Table 2: fully screened (mu_ne = mu_ee = mu).  NIST applies to the isolated ion only.
NIST_IP = {"Fe XVII": 1262.7, "Fe XVIII": 1357.8, "Fe XIX": 1460.0, "Fe XX": 1575.6, "Fe XXI": 1687.0,
           "Fe XXII": 1798.4, "Fe XXIII": 1950.4, "Fe XXIV": 2045.759, "Fe XXV": 8828.1875}
TABLE2_MU0 = {"Fe XVII": 1260.58, "Fe XVIII": 1357.09, "Fe XIX": 1459.12, "Fe XX": 1573.48, "Fe XXI": 1689.13,
              "Fe XXII": 1797.82, "Fe XXIII": 1950.49, "Fe XXIV": 2044.34, "Fe XXV": 8836.74}
TABLE2_MU0P1 = {"Fe XVII": 1214.75, "Fe XVIII": 1308.56, "Fe XIX": 1407.92, "Fe XX": 1519.58, "Fe XXI": 1632.53,
                "Fe XXII": 1738.54, "Fe XXIII": 1888.61, "Fe XXIV": 1979.79, "Fe XXV": 8768.90}
TABLE2_MU0P25 = {"Fe XVII": 1147.80, "Fe XVIII": 1237.68, "Fe XIX": 1333.09, "Fe XX": 1440.80, "Fe XXI": 1549.75,
                 "Fe XXII": 1651.80, "Fe XXIII": 1798.43, "Fe XXIV": 1885.77, "Fe XXV": 8667.86}

# Table 4: independent screening, columns (mu_ne, mu_ee) = A (0.1, 0), B (0.1, 0.1), C (0.25, 0), D (0.25, 0.25).
TABLE4 = {
    "Fe XVII":  {"A": 1190.72, "B": 1214.75, "C": 1089.22, "D": 1147.80},
    "Fe XVIII": {"A": 1287.21, "B": 1308.56, "C": 1185.55, "D": 1237.68},
    "Fe XIX":   {"A": 1389.10, "B": 1407.92, "C": 1287.41, "D": 1333.09},
    "Fe XX":    {"A": 1503.54, "B": 1519.58, "C": 1401.59, "D": 1440.80},
    "Fe XXI":   {"A": 1619.15, "B": 1632.53, "C": 1517.02, "D": 1549.75},
    "Fe XXII":  {"A": 1727.82, "B": 1738.54, "C": 1625.57, "D": 1651.80},
    "Fe XXIII": {"A": 1880.58, "B": 1888.61, "C": 1778.76, "D": 1798.43},
    "Fe XXIV":  {"A": 1974.41, "B": 1979.79, "C": 1872.58, "D": 1885.77},
    "Fe XXV":   {"A": 8766.28, "B": 8768.90, "C": 8661.14, "D": 8667.86},
}

# Published value per case: the isolated ion from Table 2 (mu = 0), screened cases from Table 4.
# (Table 4 columns B and D are identical to Table 2 at mu = 0.1 and 0.25.)
def published(ion: str, case: str) -> tuple[float, str]:
    if case == "unscreened":
        return TABLE2_MU0[ion], "Table 2 (mu=0)"
    return TABLE4[ion][case], "Table 4"
