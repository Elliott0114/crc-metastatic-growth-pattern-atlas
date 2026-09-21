# Execution Order

Run the complete route with `reproduce.py full`; the runner supplies the isolated working directory.

| Order | Module | Script |
| --- | --- | --- |
| 1 | iss_ecology | [analyze_hgp_interface_ecology.py](../analysis/analyze_hgp_interface_ecology.py) |
| 2 | iss_composition | [analyze_escrivaconde2026_slpi_state_composition.py](../analysis/analyze_escrivaconde2026_slpi_state_composition.py) |
| 3 | iss_display | [prepare_slpi_iss_publication_maps.py](../analysis/prepare_slpi_iss_publication_maps.py) |
| 4 | iss_adjacency | [analyze_iss_state_specific_liver_adjacency.py](../analysis/analyze_iss_state_specific_liver_adjacency.py) |
| 5 | iss_topology | [analyze_iss_nc3_cohesive_dual_interface.py](../analysis/analyze_iss_nc3_cohesive_dual_interface.py) |
| 6 | mapping | [analyze_cross_platform_state_anchor.py](../analysis/analyze_cross_platform_state_anchor.py) |
| 7 | rec50 | [analyze_ogden_anchor_program.R](../analysis/analyze_ogden_anchor_program.R) |
| 8 | state_pseudobulk | [prepare_ogden_direct_state_pseudobulks.py](../analysis/prepare_ogden_direct_state_pseudobulks.py) |
| 9 | state_contrasts | [analyze_ogden_direct_state_contrasts.R](../analysis/analyze_ogden_direct_state_contrasts.R) |
| 10 | state_enrichment | [run_ogden_direct_state_gsea.py](../analysis/run_ogden_direct_state_gsea.py) |
| 11 | state_source_tables | [summarize_paper_state_evidence.py](../analysis/summarize_paper_state_evidence.py) |
| 12 | c8_enrichment | [analyze_rec_c8_tight_junction_bridge.py](../analysis/analyze_rec_c8_tight_junction_bridge.py) |
| 13 | regulatory_overlap | [analyze_ogden_regulatory_junction_bridge.R](../analysis/analyze_ogden_regulatory_junction_bridge.R) |
| 14 | regulatory_coexistence | [analyze_ogden_regulatory_junction_coexistence.py](../analysis/analyze_ogden_regulatory_junction_coexistence.py) |
| 15 | discovery_coexistence | [analyze_ogden_rec_within_cell_coexistence.py](../analysis/analyze_ogden_rec_within_cell_coexistence.py) |
| 16 | clinical_input | [prepare_rec_clinical_matrix.py](../analysis/prepare_rec_clinical_matrix.py) |
| 17 | clinical | [analyze_rec_public_clinical.R](../analysis/analyze_rec_public_clinical.R) |
| 18 | protein | [analyze_rec_public_proteomics.R](../analysis/analyze_rec_public_proteomics.R) |
| 19 | specificity | [analyze_rec_public_specificity.R](../analysis/analyze_rec_public_specificity.R) |
| 20 | within_cell | [analyze_rec_public_within_cell.py](../analysis/analyze_rec_public_within_cell.py) |
| 21 | paired_states | [analyze_rec_paired_state_correlations.py](../analysis/analyze_rec_paired_state_correlations.py) |
| 22 | hgp_cells | [analyze_rec_program_hgp_scrna_projection.py](../analysis/analyze_rec_program_hgp_scrna_projection.py) |
| 23 | hgp_coexpression | [analyze_rec_epithelial_coexpression.py](../analysis/analyze_rec_epithelial_coexpression.py) |
| 24 | hgp_spatial | [analyze_rec_program_hgp_spatial_projection.py](../analysis/analyze_rec_program_hgp_spatial_projection.py) |
| 25 | hgp_display | [prepare_rec_program_spatial_display.py](../analysis/prepare_rec_program_spatial_display.py) |
| 26 | cross_modal | [analyze_rec_program_cross_modal_concordance.py](../analysis/analyze_rec_program_cross_modal_concordance.py) |
| 27 | regional_projection | [analyze_gse294385_rec_program_extension.py](../analysis/analyze_gse294385_rec_program_extension.py) |
| 28 | cross_context | [analyze_rec_program_two_axis_interpretation.py](../analysis/analyze_rec_program_two_axis_interpretation.py) |
| 29 | bulk_rec50 | [analyze_rec_program_gse151165_bulk_projection.R](../analysis/analyze_rec_program_gse151165_bulk_projection.R) |
| 30 | latacz | [analyze_rec_program_latacz_interface_summary.R](../analysis/analyze_rec_program_latacz_interface_summary.R) |
| 31 | rec44 | [analyze_rec_program_multisource_core.R](../analysis/analyze_rec_program_multisource_core.R) |
| 32 | retained_rec44_spatial | [reproduce_retained_rec44_spatial.py](../analysis/reproduce_retained_rec44_spatial.py) |
| 33 | regional_counts | [prepare_paper_regional_counts.R](../analysis/prepare_paper_regional_counts.R) |
| 34 | historical_rivals | [analyze_rec_external_rival_programs.R](../analysis/analyze_rec_external_rival_programs.R) |
| 35 | rival_clean | [summarize_paper_rivals.R](../analysis/summarize_paper_rivals.R) |
| 36 | historical_bulk | [analyze_paper_hgp_reference.R](../analysis/analyze_paper_hgp_reference.R) |
| 37 | bulk_context | [analyze_paper_bulk_context.R](../analysis/analyze_paper_bulk_context.R) |
| 38 | junction_cells | [analyze_rec_junction_cells.py](../analysis/analyze_rec_junction_cells.py) |
| 39 | junction_origins | [summarize_paper_compartments.py](../analysis/summarize_paper_compartments.py) |
| 40 | junction_context | [analyze_rec_junction_contexts.py](../analysis/analyze_rec_junction_contexts.py) |
| 41 | junction_neighbourhood | [analyze_rec_junction_neighbourhoods.py](../analysis/analyze_rec_junction_neighbourhoods.py) |
| 42 | junction_sensitivity | [analyze_rec_junction_sensitivities.py](../analysis/analyze_rec_junction_sensitivities.py) |
| 43 | junction_measurement | [audit_rec_junction_measurement.py](../analysis/audit_rec_junction_measurement.py) |
| 44 | hgp_spatial_bridge | [analyze_rec_evidence_bridge_spatial.py](../analysis/analyze_rec_evidence_bridge_spatial.py) |
| 45 | hgp_counts | [analyze_rec_evidence_bridge_counts.R](../analysis/analyze_rec_evidence_bridge_counts.R) |
| 46 | hgp_score_bridge | [summarize_paper_hgp_scores.py](../analysis/summarize_paper_hgp_scores.py) |
| 47 | matching | [analyze_rec_reinforcement_cells.py](../analysis/analyze_rec_reinforcement_cells.py) |
| 48 | same_section_input | [prepare_rec_reinforcement_sections.py](../analysis/prepare_rec_reinforcement_sections.py) |
| 49 | same_section | [analyze_rec_reinforcement_sections.R](../analysis/analyze_rec_reinforcement_sections.R) |
| 50 | canonical_bulk | [analyze_rec_reinforcement_bulk.R](../analysis/analyze_rec_reinforcement_bulk.R) |
| 51 | bulk_rec50_omissions | [summarize_paper_rec50_omissions.py](../analysis/summarize_paper_rec50_omissions.py) |
| 52 | member_summary | [summarize_paper_members.py](../analysis/summarize_paper_members.py) |
| 53 | normalization | [diagnose_rec_reinforcement_normalization.R](../analysis/diagnose_rec_reinforcement_normalization.R) |
| 54 | external | [analyze_rec_q2_external.py](../analysis/analyze_rec_q2_external.py) |
| 55 | meta | [meta_rec_q2_external.R](../analysis/meta_rec_q2_external.R) |
| 56 | regional_adjustment | [analyze_rec_q2_regional.R](../analysis/analyze_rec_q2_regional.R) |
| 57 | mrtx_pseudobulk | [prepare_rec_mapk_pseudobulk.py](../analysis/prepare_rec_mapk_pseudobulk.py) |
| 58 | mrtx | [analyze_rec_public_mapk.R](../analysis/analyze_rec_public_mapk.R) |
| 59 | plexin | [analyze_plexinb2_functional_projection.py](../analysis/analyze_plexinb2_functional_projection.py) |
| 60 | gata6 | [analyze_gata6_functional_competitor_projection.R](../analysis/analyze_gata6_functional_competitor_projection.R) |
| 61 | serial_selection | [analyze_gata6_serial_metastatic_selection_projection.R](../analysis/analyze_gata6_serial_metastatic_selection_projection.R) |

Final assembly: `materialize_sources.py`, `plotting/redraw_all.py`, `workbook_sources.py`, then `verify`.

Each module's observed files are in `metadata/module_file_dependencies.json`; the artifact-level map is `metadata/reproduction_map.tsv`.
