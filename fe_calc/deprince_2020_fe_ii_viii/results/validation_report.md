# Validation report: Deprince 2020, Fe II - Fe IX (reduced reference-configuration model)

Run on commit `38a0f45eae` (FAC source `603c0e5358`), branch `experiment/densitytosz-validation`.
Build: cflags=[-O2 -fPIC]; cc=Apple clang version 21.0.0 (clang-2100.3.34.2); ext_sha256=20cf06f7e2a368d8; BUILD_INFO-matches

| Test | Result | Seconds | Note |
|---|---|---|---|
| `AnalysisOnSavedEnergies.test_allow_incomplete_gives_exit_zero_but_still_reports` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_complete_dataset_independent_recomputation_and_no_fac_import` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_failed_and_duplicated_records_reported` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_incompatible_source_commit_blocks_the_ip` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_ip_pairs_adjacent_ions_only` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_missing_energy_fails_loudly_and_is_not_replaced` | PASS | 0.0 |  |
| `AnalysisOnSavedEnergies.test_models_never_mixed` | PASS | 0.0 |  |
| `Configurations.test_electron_counts_and_charge_states` | PASS | 0.0 |  |
| `Configurations.test_fe_ii_to_iv_equal_validated_lists` | PASS | 0.0 |  |
| `Configurations.test_fe_ix_is_paper_iii_reference_set` | PASS | 0.0 |  |
| `Configurations.test_ground_configuration_is_first_and_has_no_hole` | PASS | 0.0 |  |
| `Configurations.test_hole_notation_translated` | PASS | 0.0 |  |
| `Configurations.test_input_validation_passes` | PASS | 0.0 |  |
| `Configurations.test_ion_name_normalisation` | PASS | 0.0 |  |
| `Configurations.test_no_generated_excitations` | PASS | 0.0 |  |
| `RegressionAgainstValidatedResults.test_all_48_energies_ok_and_finite` | PASS | 0.0 |  |
| `RegressionAgainstValidatedResults.test_fe_ii_to_iv_match_linux_validated_energies` | PASS | 0.0 |  |
| `RegressionAgainstValidatedResults.test_ref_and_gref_energies_are_distinct_with_same_ci_space` | PASS | 0.0 |  |
| `RegressionAgainstValidatedResults.test_single_build_provenance_across_all_runs` | PASS | 0.0 |  |
| `PreflightAndCheckOnly.test_check_only_passes_without_fac` | PASS | 0.0 |  |
| `PreflightAndCheckOnly.test_preflight_reports_provenance_and_workaround` | PASS | 0.1 |  |
| `ResumeAndFailure.test_resume_overwrite_timeout_and_retry` | PASS | 2.6 |  |
| `SmokeTests.test_fe_ii_ip_by_independent_subtraction` | PASS | 0.0 |  |
| `SmokeTests.test_mu_zero_path_is_truly_unscreened` | PASS | 0.0 |  |
| `SmokeTests.test_nonzero_screening_applied_to_both_interactions` | PASS | 0.0 |  |
| `SmokeTests.test_raw_log_and_outputs_exist_and_are_inspected` | PASS | 0.0 |  |
| `SmokeTests.test_ref_and_gref_distinct_with_same_configuration_space` | PASS | 0.0 |  |

**27 passed, 0 failed, 0 skipped.**
