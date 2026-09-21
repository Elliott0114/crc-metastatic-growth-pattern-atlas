"""Fixed junction/HRC bridge in the six HGP-labelled FFPE spatial patients."""
import itertools
import json

import numpy as np
import pandas as pd
from scipy import sparse

from rec_junction_context_common import (
    ROOT, bh, canonical, correlations, definitions, digest, mean_summary,
    rank_design, rank_residual, rng_for,
)
from analyze_rec_program_hgp_spatial_projection import (
    INPUT, MARKERS, build_masks, distance_from_origins, load_metadata,
    marker_score, read_10x_h5, read_positions,
)

OUT = ROOT / "analysis_results/rec_evidence_bridge_2026-09-16"
COMPONENTS = ["Claudin", "Polarity", "Junction", "HRC"]


def write(frame, name):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, encoding="utf-8", na_rep="NA")


def group_summary(values, labels, name):
    values, labels = np.asarray(values, float), np.asarray(labels) == "rHGP"
    assert len(values) == 6 and labels.sum() == 3
    effect = values[labels].mean() - values[~labels].mean()
    null = []
    for selected in itertools.combinations(range(6), 3):
        mask = np.zeros(6, bool)
        mask[list(selected)] = True
        null.append(values[mask].mean() - values[~mask].mean())
    rng = rng_for(name)
    draws = (rng.choice(values[labels], (10000, 3)).mean(axis=1)
             - rng.choice(values[~labels], (10000, 3)).mean(axis=1))
    leave = []
    for i in range(6):
        keep = np.arange(6) != i
        leave.append(values[keep & labels].mean() - values[keep & ~labels].mean())
    low, high = np.quantile(draws, [.025, .975])
    return dict(n=6, effect=effect, low=low, high=high,
                p=float(np.mean(np.abs(null) >= abs(effect) - 1e-12)),
                loo_min=min(leave), loo_max=max(leave))


def main():
    meta = load_metadata()
    sets = {k: v for k, v in definitions().items() if k in COMPONENTS}
    mask_genes = set(MARKERS["epithelial"] + MARKERS["hepatocyte"])
    excluded = set(sets["HRC"]) & mask_genes
    assert excluded == {"CEACAM5"}
    sets["HRC"] = sorted(set(sets["HRC"]) - mask_genes)
    assert all(not (set(v) & mask_genes) for v in sets.values())
    write(pd.DataFrame([dict(component=k, gene=g) for k, vv in sets.items() for g in vv]),
          "spatial/definitions.tsv")

    sources, common, duplicate_rows = {}, None, []
    for sample in meta["sample"]:
        path = INPUT / sample / "outs/filtered_feature_bc_matrix.h5"
        matrix, names, barcodes = read_10x_h5(path)
        names = np.array([canonical(g) for g in names])
        assert np.all(matrix.data >= 0) and np.all(matrix.data == np.floor(matrix.data))
        detected = np.asarray((matrix > 0).sum(axis=1)).ravel()
        unique, inverse, multiplicity = np.unique(names, return_inverse=True, return_counts=True)
        for g, n in zip(unique[multiplicity > 1], multiplicity[multiplicity > 1]):
            duplicate_rows.append(dict(sample=sample, gene=g, n_features=n,
                                       target=g in set().union(*map(set, sets.values()))))
        if len(unique) != len(names):
            projection = sparse.csr_matrix((np.ones(len(names)), (np.arange(len(names)), inverse)),
                                           shape=(len(names), len(unique)))
            matrix = (matrix @ projection).tocsr()
            names = unique
        sources[sample] = (matrix, names, barcodes, detected)
        common = set(names) if common is None else common & set(names)
    common = sorted(common)
    write(pd.DataFrame(duplicate_rows), "spatial/duplicate_symbol_aggregation.tsv")
    coverage = pd.DataFrame([dict(component=k, gene=g, measured=g in common)
                            for k, vv in sets.items() for g in vv])
    write(coverage, "spatial/coverage.tsv")
    programs = {k: sorted(set(v) & set(common)) for k, v in sets.items()}
    region_meta, count_vectors, audits, all_spots, patient_cor = [], [], [], [], []
    provenance = []
    for row in meta.itertuples(index=False):
        sample, patient, hgp = row.sample, row.patient, row.hgp
        matrix, names, barcodes, detected = sources[sample]
        outs = INPUT / sample / "outs"
        positions = read_positions(outs).reindex(barcodes)
        assert positions.index.is_unique and not positions.isna().any().any()
        keep = (positions.in_tissue.to_numpy() == 1) & (detected >= 200)
        matrix, positions, barcodes = matrix[keep].tocsr(), positions.loc[barcodes[keep]], barcodes[keep]
        library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
        assert (library > 0).all()
        normalized = matrix.multiply((1e4 / library)[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        gene_indices = {g: i for i, g in enumerate(names)}
        markers = {k: marker_score(normalized, names, vv)[0]
                   for k, vv in MARKERS.items() if k in ["epithelial", "hepatocyte"]}
        rr, cc = positions.array_row.to_numpy(int), positions.array_col.to_numpy(int)
        masks = build_masks(markers, rr, cc)
        distance = distance_from_origins(rr, cc, masks["liver_side"])
        connected = masks["tumour_side"] & np.isfinite(distance)
        epihigh = masks["tumour_side"] & (markers["epithelial"] >= np.median(markers["epithelial"][masks["tumour_side"]]))
        frame = pd.DataFrame(dict(sample=sample, patient=patient, hgp=hgp, barcode=barcodes,
            array_row=rr, array_col=cc, pixel_row=positions.pxl_row_in_fullres.to_numpy(),
            pixel_col=positions.pxl_col_in_fullres.to_numpy(), tumour_side=masks["tumour_side"],
            liver_side=masks["liver_side"], epithelial_high=epihigh, distance=distance,
            Epithelial=markers["epithelial"], Liver=markers["hepatocyte"],
            log_library=np.log1p(library), log_genes=np.log1p(detected[keep])))
        for k, genes in programs.items():
            frame[k] = np.asarray(normalized[:, [gene_indices[g] for g in genes]].mean(axis=1)).ravel()
        for gene in programs["Junction"]:
            frame["gene_" + gene] = normalized[:, gene_indices[gene]].toarray().ravel()
        all_spots.append(frame)
        for restriction, selected in [("all_tumour", masks["tumour_side"]), ("epithelial_high", epihigh)]:
            d = frame.loc[selected]
            assert len(d) >= 100
            for model, covars in [("unadjusted", []), ("adjusted", ["Epithelial", "Liver", "log_library", "log_genes"])]:
                design = rank_design(d, covars)
                hrc = rank_residual(d.HRC.to_numpy(), design)
                targets = ["Claudin", "Polarity"] + ["gene_" + g for g in programs["Junction"]]
                residuals = rank_residual(d[targets].to_numpy(), design)
                for target, rho in zip(targets, correlations(residuals, hrc)):
                    patient_cor.append(dict(patient=patient, hgp=hgp, restriction=restriction,
                        model=model, component=target, n_spots=len(d), rho=rho,
                        design_rank=np.linalg.matrix_rank(design)))
        for hops, restriction in [(5, "all_tumour"), (5, "epithelial_high"), (3, "all_tumour"), (7, "all_tumour")]:
            selected = np.ones(len(frame), bool) if restriction == "all_tumour" else epihigh
            regions = {"near": connected & (distance >= 1) & (distance <= hops),
                       "deep": connected & (distance > hops), "tumour": masks["tumour_side"]}
            for region, mask in regions.items():
                mask = mask & selected
                identifier = f"{sample}|{hops}|{restriction}|{region}"
                full = np.asarray(matrix[mask].sum(axis=0)).ravel()
                count_vectors.append(pd.Series(full[[gene_indices[g] for g in common]], index=common, name=identifier))
                region_meta.append(dict(id=identifier, sample=sample, patient=patient, hgp=hgp,
                    hops=hops, restriction=restriction, region=region, n_spots=int(mask.sum()),
                    library_size=library[mask].sum(), evaluable=int(mask.sum()) >= 20))
        audits.append(dict(sample=sample, patient=patient, hgp=hgp, retained_spots=len(frame),
            tumour_side=int(masks["tumour_side"].sum()), connected_tumour=int(connected.sum()),
            epithelial_high=int(epihigh.sum()), mask_genes_target_overlap=0))
        provenance.append(dict(sample=sample, counts_sha256=digest(outs / "filtered_feature_bc_matrix.h5"),
                               positions_sha256=digest(outs / "spatial/tissue_positions.csv")))
        print(f"{sample} {patient} {hgp}: {len(frame)} spots; {connected.sum()} connected tumour-side", flush=True)
    counts = pd.concat(count_vectors, axis=1)
    metadata = pd.DataFrame(region_meta)
    assert metadata.id.is_unique and counts.columns.tolist() == metadata.id.tolist()
    write(counts.rename_axis("gene").reset_index(), "spatial/region_counts.tsv.gz")
    write(metadata, "spatial/region_metadata.tsv")
    write(pd.DataFrame(audits), "spatial/sample_qc.tsv")
    write(pd.concat(all_spots, ignore_index=True), "spatial/spot_scores.tsv.gz")
    cor = pd.DataFrame(patient_cor)
    write(cor, "spatial/patient_correlations.tsv")
    cor_sum = []
    for (restriction, model, component), d in cor.groupby(["restriction", "model", "component"]):
        cor_sum.append(dict(restriction=restriction, model=model, component=component,
                           **mean_summary(d.rho, restriction + model + component, fisher=True)))
    cor_sum = pd.DataFrame(cor_sum)
    cor_sum["family"] = np.where(cor_sum.component.str.startswith("gene_"), "genes", "components")
    cor_sum["q"] = cor_sum.groupby(["restriction", "model", "family"])["p"].transform(bh)
    write(cor_sum, "spatial/correlation_summary.tsv")

    scores, effects, genes_out = [], [], []
    for (hops, restriction), mm in metadata.groupby(["hops", "restriction"]):
        assert mm.evaluable.all(), mm.loc[~mm.evaluable].to_dict("records")
        for basis, regions in [("paired_bands", ["near", "deep"]), ("tumour_total", ["tumour"])]:
            m = mm[mm.region.isin(regions)].set_index("id")
            c = counts[m.index]
            lib = m.library_size
            # edgeR-compatible proportional prior; zero counts have a common CPM floor.
            prior = .5 * lib / lib.mean()
            expr = np.log2(c.add(prior, axis=1).div(lib + 2 * prior, axis=1) * 1e6)
            sd = expr.std(axis=1, ddof=0)
            z = expr.sub(expr.mean(axis=1), axis=0).div(sd.replace(0, np.nan), axis=0)
            for construction, matrix in [("mean_logCPM", expr), ("mean_gene_z", z)]:
                for component, gg in programs.items():
                    value = matrix.loc[gg].mean(axis=0)
                    d = m.assign(value=value).reset_index()
                    scores.extend(d.assign(component=component, score=construction, basis=basis).to_dict("records"))
                    if basis == "paired_bands":
                        w = d.pivot(index=["patient", "hgp"], columns="region", values="value").reset_index()
                        delta = w.near - w.deep
                        results = [dict(endpoint="near_minus_deep", **mean_summary(delta, str(hops) + restriction + construction + component)),
                                   dict(endpoint="HGP_near", **group_summary(w.near, w.hgp, "near" + str(hops) + restriction + construction + component)),
                                   dict(endpoint="HGP_localization", **group_summary(delta, w.hgp, "interaction" + str(hops) + restriction + construction + component))]
                    else:
                        results = [dict(endpoint="HGP_tumour", **group_summary(d.value, d.hgp, "tumour" + str(hops) + restriction + construction + component))]
                    for result in results:
                        effects.append(dict(hops=hops, restriction=restriction, score=construction,
                                            component=component, n_genes=len(gg), **result))
            if hops == 5 and restriction == "all_tumour":
                for gene in programs["Junction"]:
                    d = m.assign(value=expr.loc[gene]).reset_index()
                    if basis == "paired_bands":
                        w = d.pivot(index=["patient", "hgp"], columns="region", values="value").reset_index()
                        for endpoint, values in [("HGP_near", w.near), ("HGP_localization", w.near - w.deep)]:
                            genes_out.append(dict(gene=gene, endpoint=endpoint, **group_summary(values, w.hgp, gene + endpoint)))
                    else:
                        genes_out.append(dict(gene=gene, endpoint="HGP_tumour", **group_summary(d.value, d.hgp, gene + "tumour")))
    effects = pd.DataFrame(effects)
    effects["q"] = effects.groupby(["hops", "restriction", "score", "endpoint"])["p"].transform(bh)
    write(pd.DataFrame(scores), "spatial/region_scores.tsv")
    write(effects, "spatial/region_effects.tsv")
    genes_out = pd.DataFrame(genes_out)
    genes_out["q"] = genes_out.groupby("endpoint")["p"].transform(bh)
    write(genes_out, "spatial/gene_effects.tsv")
    (OUT / "spatial/provenance.json").write_text(json.dumps(dict(
        spec_sha256=digest(ROOT / "metadata/rec_evidence_bridge_2026-09-16.initial.md"),
        common_features=len(common), HRC_mask_exclusion=sorted(excluded),
        sources=provenance), ensure_ascii=False, indent=2), encoding="utf-8")
    print(effects[(effects.hops == 5) & (effects.restriction == "all_tumour")
                  & (effects.score == "mean_logCPM")].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
