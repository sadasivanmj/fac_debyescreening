# Validation report: Deprince 2019, Fe XVII - Fe XXVI (reduced ground-only / reference-only models)

Run on commit `38a0f45eae` (FAC source `603c0e5358`), branch `experiment/densitytosz-validation`.
Build: cflags=[-O2 -fPIC]; cc=Apple clang version 21.0.0 (clang-2100.3.34.2); ext_sha256=20cf06f7e2a368d8; BUILD_INFO-matches

| Test | Result | Seconds | Note |
|---|---|---|---|
| `AnalysisOnSavedEnergies.test_complete_dataset_independent_recomputation_and_no_fac_import` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_ee_effect_and_shifts_and_published_columns` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_failed_and_duplicated_records_reported` | PASS | 0.1 |  |
| `AnalysisOnSavedEnergies.test_fe_xxv_ip_uses_matching_fe_xxvi_case` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_incompatible_models_or_builds_block_the_ip` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_missing_energy_fails_loudly_and_models_not_substituted` | PASS | 0.1 |  |
| `AnalysisOnSavedEnergies.test_models_are_reported_separately` | PASS | 0.1 |  |
| `Configurations.test_electron_counts_and_charge_states` | PASS | 0.0 |  |
| `Configurations.test_fe_xxiv_is_all_three_electron_occupations_of_1s_2s_2p` | PASS | 0.0 |  |
| `Configurations.test_input_validation_passes` | PASS | 0.0 |  |
| `Configurations.test_ion_name_normalisation` | PASS | 0.0 |  |
| `Configurations.test_models_are_distinct_and_reduced` | PASS | 0.0 |  |
| `Configurations.test_published_reference_lists` | PASS | 0.0 |  |
| `Configurations.test_screening_cases_keep_nuclear_and_ee_separate` | PASS | 0.0 |  |
| `RegressionAndSanity.test_every_run_ok_converged_and_single_build` | PASS | 0.0 |  |
| `RegressionAndSanity.test_ground_only_matches_earlier_validated_energies` | PASS | 0.0 |  |
| `RegressionAndSanity.test_ionization_potentials_positive_and_decrease_with_screening` | PASS | 0.0 |  |
| `RegressionAndSanity.test_reference_only_lowers_energy_and_is_unchanged_only_where_symmetry_forbids_mixing` | PASS | 0.0 |  |
| `RegressionAndSanity.test_screening_ordering_for_every_ion` | PASS | 0.0 |  |
| `PlanningAndPreflight.test_check_only_passes_without_fac` | PASS | 0.0 |  |
| `PlanningAndPreflight.test_preflight_reports_provenance_and_workaround` | PASS | 0.1 |  |
| `PlanningAndPreflight.test_selection_of_ions_cases_models` | PASS | 0.0 |  |
| `ResumeAndFailure.test_resume_overwrite_timeout_and_retry` | PASS | 1.3 |  |
| `SmokeTests.test_fe_xvii_runs_ok_with_raw_outputs` | PASS | 0.0 |  |
| `SmokeTests.test_fe_xxv_ip_from_matching_fe_xxvi_energy` | PASS | 0.0 |  |
| `SmokeTests.test_fe_xxvi_responds_to_nuclear_screening_only` | PASS | 0.0 |  |
| `SmokeTests.test_ground_only_and_reference_only_are_distinct_models` | PASS | 0.0 |  |
| `SmokeTests.test_nuclear_only_versus_fully_screened_physics` | PASS | 0.0 |  |
| `SmokeTests.test_unscreened_path_is_plain_fac` | PASS | 0.0 |  |

**29 passed, 0 failed, 0 skipped.**
