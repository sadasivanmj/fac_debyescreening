# Radial-solver diagnostics (native SFAC)

| case | binary | rc | status | E0 (eV) | r·V(r>50) max | ips | 3d3/2 (Ha) | last log line | note |
|---|---|---:|---|---:|---:|---:|---:|---|---|
| A1_orig_refs_mu0.1 | installed | 1 | FAILED |  | 14.515 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | original driver input (nps=0.01 A^-3, cold start) |
| A2_ground_mu0.1 | installed | 1 | FAILED |  | 14.514 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | ground configuration only |
| A3_refs1-2_mu0.1 | installed | 1 | FAILED |  | 14.514 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | first 2 Table 1 references |
| A3_refs1-3_mu0.1 | installed | 1 | FAILED |  | 14.514 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | first 3 Table 1 references |
| A3_refs1-4_mu0.1 | installed | 1 | FAILED |  | 14.515 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | first 4 Table 1 references |
| A3_refs1-5_mu0.1 | installed | 1 | FAILED |  | 14.515 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | first 5 Table 1 references |
| A3_refs1-6_mu0.1 | installed | 1 | FAILED |  | 14.515 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | first 6 Table 1 references |
| A3_refs1-7_mu0.1 | installed | 1 | FAILED |  | 14.515 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | first 7 Table 1 references |
| B1_en_only_mu0.1 | installed | 1 | FAILED |  | 24.000 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 -1   0.0000000000E+00 | e-n screening only (radial:ee_screen 0) |
| B2_en_ee_mu0.1 | installed | 1 | FAILED |  | 14.514 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | e-n and e-e screening |
| C1_nps1e-2 | installed | 1 | FAILED |  | 14.514 | 278 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | rps=5.44 a0: screened u(r) frozen beyond rps |
| C2_nps1e-6 | installed | 1 | FAILED |  | 0.000 | 581 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | rps=117 a0: potential decays correctly, still fails |
| D1_fix_cold_mu0.1 | scf_fix | 1 | FAILED |  | 0.000 | 581 | -1e+30 | Error occured in RadialSolver: -6 3 2   0.0000000000E+00 | patched solver, cold start: genuine unbound transient at SCF iter 2 |
| D2_fix_ladder_mu0.1 | scf_fix | 0 | ok | -33635.8239 | 0.000 | 0 | -0.505708232 |  | patched solver, free start + warm ladder |
| D3_fix_ladder_mu0.25 | scf_fix | 0 | ok | -32293.7705 | 0.000 | 0 | -0.330006095 |  | patched solver, ladder to 0.25 |
| D4_fix_ladder_refs_mu0.25 | scf_fix | 0 | ok | -32282.8180 | 0.000 | 0 | -0.466353059 |  | patched solver, all Table 1 refs, ladder to 0.25 |
| D5_unpatched_ladder_mu0.1 | scf | 1 | FAILED |  | 0.000 | 0 | -1e+30 | Error occured in RadialSolver: -2 3 2   0.0000000000E+00 | unpatched solver, ladder: hxs_iter0 unknown option, iteration-0 exchange-free potential |
| D6_fix_ladder_hxs0_mu0.1 | scf_fix | 1 | FAILED |  | 0.000 | 0 | -1e+30 | Error occured in RadialSolver: -6 3 2   0.0000000000E+00 | patched, ladder without ramp offset |
| E1_free_unpatched | scf | 0 | ok | -34561.2842 | -2.000 | 0 | -0.81323386 |  | free ion, unpatched |
| E2_free_patched | scf_fix | 0 | ok | -34561.2842 | -2.000 | 0 | -0.81323386 |  | free ion, patched |
| E3_free_installed_CE | installed | 0 | ok | -34561.2841 | -2.000 | 0 | -0.81323386 |  | original baseline with ConfigEnergy (expect -34561.2841 eV) |
| E4_zps3_mu0.1 | scf_fix | 0 | ok | -33635.8239 | 0.000 | 0 | -0.505708232 |  | zps=3 instead of 1: must not change Debye-mode results |
| E5_nps1e-4_mu0.1 | scf_fix | 0 | ok | -33635.8239 | 0.000 | 0 | -0.505708232 |  | nps=1e-4 (T scaled): lambda_D unchanged |
| E6_mu1e-4_ladder | scf_fix | 0 | ok | -34570.3335 | -1.690 | 0 | -0.668607529 |  | large but finite lambda_D = 1e4 a0 |
| F1_ee_branch_ladder_mu0.1 | ee_fix | 1 | FAILED |  | 24.000 | 0 | -1e+30 | Error occured in RadialSolver: -6 3 2   0.0000000000E+00 | fix-debye-ee-screening + patch: DensityToSZ unscreened |
| F2_ee_branch_free | ee_fix | 0 | ok | -34561.2842 | -2.000 | 0 | -0.81323386 |  | fix-debye-ee-screening free |
