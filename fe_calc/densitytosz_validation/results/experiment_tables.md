## Total energies and A/B differences (eV)

| Ion | mu | E(A) | E(B) | dE = E(B) - E(A) |
|---|---|---|---|---|
| Fe II | 0.0 | -34571.2747 | -34571.2747 | 0.0000 |
| Fe II | 0.1 | — | -33635.8239 | — |
| Fe II | 0.25 | — | -32293.7705 | — |
| Fe III | 0.0 | -34555.7544 | -34555.7544 | 0.0000 |
| Fe III | 0.1 | — | -33625.5981 | — |
| Fe III | 0.25 | — | -32289.0833 | — |
| Fe IV | 0.0 | -34527.6602 | -34527.6602 | 0.0000 |
| Fe IV | 0.1 | — | -33605.4643 | — |
| Fe IV | 0.25 | — | -32279.8135 | — |

(A is missing wherever its SCF did not converge; nothing is substituted.)

## Ionization potentials (eV)

| Transition | mu | IP(A) | IP(B) | paper | err A | err B | |err| change A->B | NIST | shift A | shift B | shift paper |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Fe II -> Fe III | 0.0 | 15.5203 | 15.5203 | 18.89 | -3.3697 | -3.3697 | 0.0000 | 16.1992 | 0.0000 | 0.0000 | 0.0 |
| Fe II -> Fe III | 0.1 | — | 10.2258 | 14.2 | — | -3.9742 | — | — | — | -5.2945 | -4.69 |
| Fe II -> Fe III | 0.25 | — | 4.6872 | 7.78 | — | -3.0928 | — | — | — | -10.8331 | -11.11 |
| Fe III -> Fe IV | 0.0 | 28.0942 | 28.0942 | 32.93 | -4.8358 | -4.8358 | 0.0000 | 30.651 | 0.0000 | 0.0000 | 0.0 |
| Fe III -> Fe IV | 0.1 | — | 20.1338 | 25.6 | — | -5.4662 | — | — | — | -7.9604 | -7.33 |
| Fe III -> Fe IV | 0.25 | — | 9.2698 | 16.42 | — | -7.1502 | — | — | — | -18.8244 | -16.51 |

## Supplementary: small-mu ladder where variant A converges, FAC default radial:csi (single SCF step from the free ion)

| mu | dE Fe II | dE Fe III | dE Fe IV | dI Fe II | dI Fe III | A converged (II/III/IV) | shift A Fe II | shift B Fe II |
|---|---|---|---|---|---|---|---|---|
| 0.001 | -0.085 | -0.110 | -0.093 | -0.025 | 0.017 | y/y/y | -0.055 | -0.080 |
| 0.003 | -0.208 | -0.297 | -0.254 | -0.089 | 0.043 | y/y/y | -0.163 | -0.252 |
| 0.005 | -0.214 | -0.441 | -0.389 | -0.227 | 0.052 | y/y/y | -0.271 | -0.498 |
| 0.01 | -0.173 | -0.578 | -0.558 | -0.404 | 0.020 | y/y/y | -0.538 | -0.943 |
| 0.02 | — | -0.587 | -0.657 | — | -0.070 | n/y/y | — | — |
| 0.03 | — | — | -0.625 | — | — | n/n/y | — | — |
| 0.05 | — | — | -0.725 | — | — | n/n/y | — | — |

## Supplementary: the same ladder with radial:csi = 0 in both builds (removes FAC's Coulomb-designed asymptotic-charge adjustment)

| mu | dE Fe II | dE Fe III | dE Fe IV | dI Fe II | dI Fe III | A converged (II/III/IV) |
|---|---|---|---|---|---|---|
| 0.001 | -0.000 | 0.000 | -0.000 | 0.000 | -0.000 | y/y/y |
| 0.003 | -0.001 | -0.000 | -0.000 | 0.001 | 0.000 | y/y/y |
| 0.005 | 0.060 | 0.002 | -0.001 | -0.058 | -0.003 | y/y/y |
| 0.01 | 0.089 | 0.077 | 0.051 | -0.011 | -0.026 | y/y/y |
| 0.02 | — | 0.076 | 0.071 | — | -0.005 | n/y/y |
| 0.03 | — | — | 0.065 | — | — | n/n/y |
| 0.05 | — | — | — | — | — | n/n/n |

## Supplementary: Fe XVII - Fe XXVI, electron-electron screening effect split into Slater and direct-potential parts (eV)

| Ion | mu | dE(e-e screening), A build (Slater only) | dE(e-e screening), B build (Slater + direct) | direct-potential part |
|---|---|---|---|---|
| Fe XX | 0.1 | -56.3436 | -56.4258 | -0.0822 |
| Fe XX | 0.25 | -137.9379 | -137.9506 | -0.0127 |
| Fe XXV | 0.1 | -2.7096 | -2.7095 | 0.0001 |
| Fe XXV | 0.25 | -6.7315 | -6.7328 | -0.0013 |
| Fe XXI | 0.1 | -40.2683 | -40.2941 | -0.0258 |
| Fe XXI | 0.25 | -98.6661 | -98.5876 | 0.0785 |
| Fe XIX | 0.1 | -75.0917 | -75.1039 | -0.0122 |
| Fe XIX | 0.25 | -183.7160 | -183.6383 | 0.0777 |
| Fe XXII | 0.1 | -26.8634 | -26.8909 | -0.0275 |
| Fe XXII | 0.25 | -65.8863 | -65.8637 | 0.0226 |
| Fe XXIV | 0.1 | -8.0851 | -8.0662 | 0.0189 |
| Fe XXIV | 0.25 | -19.9253 | -19.8810 | 0.0443 |
| Fe XVII | 0.1 | -120.5816 | -120.5824 | -0.0008 |
| Fe XVII | 0.25 | -294.6421 | -294.6403 | 0.0018 |
| Fe XXVI | 0.1 | 0.0000 | 0.0000 | 0.0000 |
| Fe XXVI | 0.25 | 0.0000 | 0.0000 | 0.0000 |
| Fe XXIII | 0.1 | -16.1338 | -16.1321 | 0.0017 |
| Fe XXIII | 0.25 | -39.6274 | -39.6025 | 0.0249 |
| Fe XVIII | 0.1 | -96.5054 | -96.4950 | 0.0104 |
| Fe XVIII | 0.25 | -235.9571 | -235.9139 | 0.0432 |
