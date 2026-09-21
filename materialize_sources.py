"""Rebuild retained display tables from fresh analysis outputs, without old workbooks."""
from pathlib import Path
import json
import shutil

import numpy as np
import pandas as pd

LEGACY = Path("manuscript/rec_discovery_reinforced_2026-09-17/source_data/preserved_figure_inputs")
PHASE = Path("analysis_results/deep_biology_upgrade_2026-08-31")
P2 = PHASE / "phase2_mechanistic_specificity"
REG = P2 / "ogden_regulatory_junction_bridge"
PUBLIC = Path("analysis_results/rec_public_upgrade_2026-09-07")
PAPER = Path("manuscript/rec_text_polish_2026-09-21")


def materialize(root, workspace):
    used = []
    current_sources = []
    producers = {}

    def read(relative):
        current_sources[:] = [str(relative)]
        used.append(str(relative))
        return pd.read_csv(workspace / relative, sep="\t")

    def write(relative, frame):
        path = workspace / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False, encoding="utf-8", na_rep="NA")
        producers[str(relative)] = list(current_sources)

    def save(relative, frame):
        write(LEGACY / relative, frame)

    for name in ("figure1c", "figure1d"):
        source = "patient_composition.tsv" if name == "figure1c" else "patient_adjacency.tsv"
        data = read(Path("analysis_results/hgp_interface_ecology") / source)
        data = data[~data.canonical_patient.isin(["P03", "P05"])].copy()
        if name == "figure1d":
            data = data[data.definition == "liver_epithelial_k10"].copy()
        data["hgp"] = data.hgp.replace({"EHGP": "dHGP", "RHGP": "rHGP"})
        save(f"Figure_1/{name}.tsv", data)
    data = read(Path("analysis_results/slpi_bjc2026_state_composition/patient_state_endpoints.tsv"))
    data = data[~data.canonical_patient.isin(["P03", "P05"])]
    rows = [dict(canonical_patient=r.canonical_patient, hgp={"EHGP": "dHGP", "RHGP": "rHGP"}[r.hgp], state=state, fraction=getattr(r, "prop_" + state)) for r in data.itertuples() for state in ("NC1", "NC2", "NC3", "NC\u03b3")]
    save("Figure_1/figure1b.tsv", pd.DataFrame(rows))
    data = read(Path("analysis_results/slpi_publication_figure_augmentation/representative_spatial_cells.tsv.gz"))
    data["is_tumour"] = data.Tumor.astype(bool)
    data["display_class_new"] = np.select([data.is_tumour & data.Clusters.eq("NC3"), data.is_tumour & data.Clusters.isin(["NC1", "NC2", "NC\u03b3"]), data.Clusters.isin(["HC1", "DHC1", "DHC2", "DHC3", "CC1"])], ["NC3", "Other tumour cells", "Liver epithelium"], default="Other tissue")
    data["display_class_new"] = pd.Categorical(data.display_class_new, categories=["Other tissue", "Liver epithelium", "Other tumour cells", "NC3"], ordered=True)
    save("Figure_1/figure1a.tsv", data.sort_values("display_class_new", kind="stable"))
    data = read(P2 / "iss_cohesive_dual_interface/iss_cohesive_dual_interface_group_descriptive.tsv")
    data = data[(data.k == 10) & (data.state == "NC3") & data.metric.isin(["homotypic_any_neighbour", "dhc_dual_interface", "hc1_dual_interface"]) & (data.value_type == "spatial_null_residual")]
    save("Figure_1/figure_phase2_topology_residual_summary_source_data.tsv", data)
    data = read(Path("analysis_results/hgp_interface_ecology/cluster_effects.tsv"))
    data["p_label"] = data.BH_FDR_across_four_clusters.map(lambda p: f"FDR {p:.3f}")
    save("supplementary_figure_1/figureS1a.tsv", data)
    data = read(Path("analysis_results/cross_platform_state_anchor/state_similarity.tsv"))
    data["label"] = data.pearson_correlation.map(lambda x: f"{x:.2f}")
    save("Figure_2/figure2a.tsv", data)
    aliases = {
        "Figure_2/leave_one_patient_out_mapping.tsv": Path("analysis_results/cross_platform_state_anchor/leave_one_patient_out_mapping.tsv"),
        "Figure_3/ogden_program_paired_contrasts.tsv": PUBLIC / "specificity/ogden_program_paired_contrasts.tsv",
        "Figure_3/rec_c8_tight_junction_bridge_evidence.tsv": P2 / "rec_c8_tight_junction_bridge_evidence.tsv",
        "Supplementary_Figure_7/RNA18_REC50_concordance.tsv": PUBLIC / "specificity/RNA18_REC50_concordance.tsv",
        "supplementary_figure_3/figure4a.tsv": Path("analysis_results/rec_program_hgp_scrna_projection/patient_scores.tsv"),
    }
    for target, source in aliases.items():
        save(target, read(source))
    write(Path("manuscript/rec_reviewer_revision_2026-09-07/source_data/Figure_3/table_phase2_direct_program_gsea.tsv"), read(P2 / "table_phase2_direct_program_gsea.tsv"))
    data = read(Path("analysis_results/rec_program_hgp_spatial_projection/patient_region_scores.tsv"))
    data = data[(data.qc_min_detected_genes == 200) & (data.restriction == "all_tumour_side")].copy()
    near = data[data.region == "near_interface_0_500um"]
    rows = [dict(patient=r.patient, hgp=r.hgp, endpoint=endpoint, raw_value=getattr(r, endpoint)) for r in near.itertuples() for endpoint in ["rec_program_score", "pooled_program_log2_cpm"]]
    abundance = pd.DataFrame(rows)
    abundance["standardized_abundance"] = abundance.groupby("endpoint").raw_value.transform(lambda v: (v-v.mean())/v.std(ddof=1))
    abundance["endpoint_label"] = abundance.endpoint.map({"rec_program_score": "Equal-gene score", "pooled_program_log2_cpm": "Pooled expression"})
    save("supplementary_figure_3/figure4b.tsv", abundance)
    data["region"] = data.region.replace({"near_interface_0_500um": "Near interface", "deep_tumour_gt500um": "Deep tumour"})
    save("supplementary_figure_3/figure4d.tsv", data)
    for name in ["representative_spatial_spots.tsv.gz", "representative_spatial_summary.tsv"]:
        data = read(Path("analysis_results/rec_program_spatial_display") / name)
        data["image_path"] = data["sample"] + "_tissue_hires_image.png"
        save("supplementary_figure_3/" + name, data)
    regulators = ["JUNB", "AP1_COMBINED", "FOSL2", "NFKB_COMBINED", "RELB", "TEAD1", "HNF4A", "CDX2"]
    data = read(REG / "ogden_regulon_tight_junction_ora.tsv").set_index("regulon_id").reindex(regulators).reset_index()
    save("supplementary_figure_6/figure_phase2_regulatory_target_overlap_source_data.tsv", data)
    for source, name, threshold in [("ogden_regulon_paired_state_summary.tsv", "state_activity", None), ("ogden_regulatory_junction_patient_state_correlations.tsv", "coexistence_patient", "n_cells"), ("ogden_regulatory_junction_patient_equal_summaries.tsv", "coexistence_summary", "minimum_cells")]:
        data = read(REG / source)
        data = data[data.regulon_id.isin(regulators)]
        if threshold:
            data = data[(data.method == "matched_residual") & (data.state == "REC") & ((data[threshold] >= 20) if threshold == "n_cells" else (data[threshold] == 20))]
        save(f"supplementary_figure_6/figure_phase2_regulatory_{name}_source_data.tsv", data)
    data = read(Path("analysis_results/rec_program_multisource_core/gene_recurrence_map.tsv"))
    historical = read(Path("analysis_results/retained_rec44_historical/spatial_gene_effects.tsv")).set_index("gene")
    data["spatial_rHGP_vs_dHGP"] = data.gene.map(historical.spatial_rhgp_minus_dhgp)
    hgp_columns = ["bulk_rHGP_vs_dHGP", "interface_rHGP_vs_dHGP", "epithelial_rHGP_vs_dHGP", "spatial_rHGP_vs_dHGP"]
    data["hgp_positive_sources"] = (data[hgp_columns] > 0).sum(axis=1)
    data["all_positive_sources"] = data.hgp_positive_sources + (data.outgrowth_macro_vs_micro > 0)
    data["persistent_five_source_core"] = data.all_positive_sources == 5
    data["hgp_recurrence_class"] = data.hgp_positive_sources.map({4: "positive in all four HGP datasets", 3: "positive in three HGP datasets", 2: "split HGP direction"}).fillna("positive in zero or one HGP dataset")
    data["shared_direction"] = np.where((data.bulk_rHGP_vs_dHGP > 0) & (data.interface_rHGP_vs_dHGP > 0), "Higher in both", "Other direction")
    data["label_gene"] = data.persistent_five_source_core | data.gene.isin(["DUOX2", "TSPAN1", "MUC13"])
    current_sources[:] = ["analysis_results/rec_program_multisource_core/gene_recurrence_map.tsv", "analysis_results/retained_rec44_historical/spatial_gene_effects.tsv"]
    write(Path("manuscript/rec_molecular_oncology_final_2026-09-09/source_data/supplementary_figure_2/figure3e.tsv"), data)
    functional = functional_context(read)
    current_sources[:] = functional.attrs['source_tables']
    save("Supplementary_Figure_9/figure_phase7_functional_context_source_data.tsv", functional)
    save("Supplementary_Figure_9/displayed_Plexin_effects.tsv", functional[functional.context == "Acute Plexin B2"])
    data = read(PUBLIC / "mapk/program_intervention_effects.tsv")
    save("Supplementary_Figure_9/displayed_MRTX_effects.tsv", data[data.selection == "epithelial"])
    notes = pd.read_csv(root / "metadata/source_selection_notes.tsv", sep="\t")
    eligibility = read(Path("analysis_results/rec_q2_integration_2026-09-19/external/patient_eligibility.tsv"))
    for i, row in notes.iterrows():
        selected = eligibility[eligibility.study == row.study]
        if len(selected) != row.local_candidate_patients or selected.n_cells.sum() != row.local_candidate_cells:
            raise ValueError(f"Source selection inventory changed: {row.study}")
    current_sources.append('metadata/source_selection_notes.tsv')
    write(PAPER / "source_data/Source_selection_notes.tsv", notes)
    # These are image/geometry inputs, not fitted model results.
    for row in json.loads((root / "metadata/static_display_inputs.json").read_text(encoding="utf-8")):
        target = workspace / row["workspace_path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / row["path"], target)
    (workspace / "display_source_reads.json").write_text(json.dumps(sorted(set(used)), indent=2), encoding="utf-8")
    (workspace / "display_source_producers.json").write_text(json.dumps(producers, indent=2), encoding="utf-8")


def functional_context(read):
    names = {"FROZEN_REC_TOP50": "REC programme", "REC_PROGRAM_TOP50": "REC programme", "CORE_HRC": "Core HRC", "PARTIAL_EMT": "Partial EMT", "REACTOME_TIGHT_JUNCTION_INTERACTIONS": "Broad tight junction", "JUNB_TIGHT_JUNCTION_TARGETS_5": "Focused junction five", "JUNB_TIGHT_JUNCTION_TARGETS_5_POSTHOC": "Focused junction five", "HALLMARK_E2F_TARGETS": "E2F targets", "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint"}
    specs = [
        ("phase3_external_regulatory_projection/gse151165_regulatory_module_contrasts.tsv", "Patient rHGP\nvs dHGP", "rhgp_minus_dhgp", "15 patient tumours; HGP-labelled", "patient-level"),
        ("phase6_gata6_serial_metastatic_selection_projection/gse290752_programme_contrasts.tsv", "Serial liver-metastatic\nselection", "adjusted_metastatic_minus_primary", "12 libraries; generation-adjusted", "library-level"),
        ("phase5_gata6_functional_competitor_projection/gse290753_programme_contrasts.tsv", "GATA6 knockout", "gata6_ko_minus_control", "3 knockout + 3 control libraries", "library-level"),
        ("phase4_plexinb2_functional_projection/gse267981_endpoint_effects.tsv", "Acute Plexin B2", "mean_difference", "1 library per condition; descriptive", "cell-summary only"),
    ]
    rows = []
    for source, context, effect, evidence, inference in specs:
        frame = read(PHASE / source)
        acute = context == "Acute Plexin B2"
        if acute:
            frame = frame[(frame.qc_definition == "primary_qc") & (frame.endpoint_type == "programme")].rename(columns={"endpoint_id": "set_id"})
        for _, row in frame[frame.set_id.isin(names)].iterrows():
            rows.append(dict(context=context, programme=names[row.set_id], source_set_id=row.set_id, effect=row[effect], lower=np.nan if acute else row.bootstrap_ci_lower, upper=np.nan if acute else row.bootstrap_ci_upper, p_value=np.nan if acute else row.exact_p_two_sided, evidence=evidence, inference=inference, post_inspection=acute and names[row.set_id] == "Focused junction five"))
    frame = pd.DataFrame(rows)
    frame.attrs['source_tables'] = [str(PHASE / source) for source, *_ in specs]
    return frame
