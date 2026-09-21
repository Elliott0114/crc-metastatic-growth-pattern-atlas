"""Measured Visium programme scores and histological neighbourhoods, patient inference."""
import gzip

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from rec_junction_context_common import (
    ROOT, bh, canonical, correlations, definitions, mean_summary,
    rank_design, rank_residual, write,
)


def main():
    sets = definitions()
    base = ROOT / "data_sources/Liu_2026_GSE294385"
    anno = pd.read_csv(base / "visium_liver_meta_after_qc.tsv.gz", sep="\t")
    manifest = pd.read_csv(ROOT / "metadata/gse294385_liver_sample_manifest.tsv", sep="\t")
    manifest = manifest[manifest.selected_for_paired_extension.astype(str).str.lower() == "true"]
    features = {}
    for sample in manifest["sample"]:
        f = pd.read_csv(base / "extracted" / sample / "filtered_feature_bc_matrix/features.tsv.gz", sep="\t", header=None)
        features[sample] = np.array(list(map(canonical, f[1])))
    common = set.intersection(*[set(v) for v in features.values()])
    programs = {k: sorted(set(v) & common) for k, v in sets.items()}
    write(pd.DataFrame([dict(component=k, gene=g, measured_in_every_section=g in common)
                       for k, vv in sets.items() for g in vv]), "spatial/gene_coverage.tsv")
    all_scores, audits = [], []
    for sample in manifest["sample"]:
        folder = base / "extracted" / sample
        barcodes = pd.read_csv(folder / "filtered_feature_bc_matrix/barcodes.tsv.gz", sep="\t", header=None)[0]
        ann = anno[anno["sample"] == sample].set_index("spot_barcode")
        assert ann.index.is_unique and ann.index.isin(barcodes).all()
        path = folder / "spatial/tissue_positions.csv"
        if path.exists():
            pos = pd.read_csv(path).set_index("barcode")
        else:
            pos = pd.read_csv(folder / "spatial/tissue_positions_list.csv", header=None,
                names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]).set_index("barcode")
        assert ann.index.isin(pos.index).all() and barcodes.isin(pos.index).all()
        pos = pos.loc[barcodes]
        coords = np.column_stack([pos.array_col / 2, pos.array_row * np.sqrt(3) / 2])
        tree = cKDTree(coords)
        neighbourhoods = tree.query_ball_point(coords, 1.01)
        labels = ann.Layer3.reindex(barcodes).fillna("unannotated").to_numpy()
        tumour = np.isin(labels, ["Liver micrometastasis tumor", "Liver macrometastasis tumor"])
        rows = np.flatnonzero(tumour)
        with gzip.open(folder / "filtered_feature_bc_matrix/matrix.mtx.gz", "rb") as f, threadpool_limits(limits=2):
            counts = mmread(f).tocsr().astype(float)
        assert counts.shape == (len(features[sample]), len(barcodes))
        library = np.asarray(counts.sum(axis=0)).ravel()
        detected = np.asarray((counts > 0).sum(axis=0)).ravel()
        assert np.all(library > 0)
        normalized = counts.multiply(1e4 / library).tocsr()
        normalized.data = np.log1p(normalized.data)
        frame = pd.DataFrame({"sample": sample, "patient": ann.patient.iloc[0],
            "barcode": barcodes.iloc[rows].to_numpy(), "region": labels[rows],
            "log_library": np.log1p(library[rows]), "log_genes": np.log1p(detected[rows]),
            "platform_IX": float(sample.startswith("IX_"))})
        for k in ["Claudin", "Polarity", "HRC", "Epithelial", "Liver"]:
            gene_values = []
            for gene in programs[k]:
                indices = np.flatnonzero(features[sample] == gene)
                if len(indices) == 1:
                    gene_values.append(normalized[indices[0], rows].toarray().ravel())
                else:
                    raw = np.asarray(counts[indices][:, rows].sum(axis=0)).ravel()
                    gene_values.append(np.log1p(raw / library[rows] * 1e4))
            frame[k] = np.mean(gene_values, axis=0)
        context = []
        for i in rows:
            neighbours = [j for j in neighbourhoods[i] if j != i and labels[j] != "unannotated"]
            lab = labels[neighbours]
            context.append(dict(n_neighbours=len(neighbours),
                neighbour_liver=np.mean(lab == "Hepatic lobule") if len(lab) else np.nan,
                neighbour_stroma=np.mean(np.char.find(lab.astype(str), "stroma") >= 0) if len(lab) else np.nan))
        frame = pd.concat([frame, pd.DataFrame(context)], axis=1)
        frame["eligible"] = frame.n_neighbours >= 3
        all_scores.append(frame)
        audits.append(dict(sample=sample, patient=ann.patient.iloc[0], n_matrix_spots=len(barcodes),
                           n_annotated=len(ann), n_tumour=len(rows), n_eligible=int(frame.eligible.sum()),
                           max_grid_neighbours=max(len(v)-1 for v in neighbourhoods)))
        print(f"{sample}: {len(rows)} tumour spots, {frame.eligible.sum()} with >=3 annotated neighbours", flush=True)
    frame = pd.concat(all_scores, ignore_index=True)
    write(frame, "spatial/tumour_spot_scores.tsv.gz")
    write(pd.DataFrame(audits), "spatial/coordinate_annotation_audit.tsv")
    frame = frame[frame.eligible].copy()
    frame["macro"] = (frame.region == "Liver macrometastasis tumor").astype(float)
    frame["state"] = frame["sample"] + ":" + frame.region
    covariates = ["HRC", "Epithelial", "Liver", "log_library", "log_genes"]
    patient_rows = []
    for patient, sub in frame.groupby("patient"):
        for neighbour in ["neighbour_liver", "neighbour_stroma"]:
            other = "neighbour_stroma" if neighbour == "neighbour_liver" else "neighbour_liver"
            design = rank_design(sub, covariates + [other], state=True)
            residual = rank_residual(sub[[neighbour, "Claudin", "Polarity"]].to_numpy(float), design)
            rho = correlations(residual[:, 1:], residual[:, 0])
            for i, component in enumerate(["Claudin", "Polarity"]):
                patient_rows.append(dict(patient=patient, neighbour=neighbour, component=component,
                    n_spots=len(sub), n_nonzero_context=int((sub[neighbour] > 0).sum()), rho=rho[i]))
    patient = pd.DataFrame(patient_rows)
    write(patient, "spatial/patient_partial_correlations.tsv")
    summary = pd.DataFrame([dict(neighbour=nb, component=cp, **mean_summary(g.rho, nb+cp, fisher=True))
                          for (nb, cp), g in patient.groupby(["neighbour", "component"])])
    summary["q_within_four_tests"] = bh(summary.p)
    write(summary, "spatial/context_summary.tsv")
    predictions = []
    baseline = covariates + ["macro", "platform_IX"]
    for patient in sorted(frame.patient.unique()):
        train, test = frame[frame.patient != patient], frame[frame.patient == patient]
        weights = 1 / train.patient.map(train.patient.value_counts()).to_numpy()
        for component in ["Claudin", "Polarity"]:
            row = dict(patient=patient, component=component, n_spots=len(test))
            for model, columns in [("baseline", baseline), ("context", baseline + ["neighbour_liver", "neighbour_stroma"])]:
                scaler = StandardScaler().fit(train[columns], sample_weight=weights)
                fit = LinearRegression().fit(scaler.transform(train[columns]), train[component], sample_weight=weights)
                pred = fit.predict(scaler.transform(test[columns]))
                row["rmse_"+model] = np.sqrt(np.mean((pred-test[component].to_numpy())**2))
            row["relative_rmse_improvement"] = 1-row["rmse_context"]/row["rmse_baseline"]
            predictions.append(row)
    pred = pd.DataFrame(predictions)
    write(pred, "spatial/held_out_patient_prediction.tsv")
    summary_pred = pd.DataFrame([dict(component=c, **mean_summary(g.relative_rmse_improvement, "predict"+c))
                                for c,g in pred.groupby("component")])
    summary_pred["q_within_two_components"] = bh(summary_pred.p)
    write(summary_pred, "spatial/prediction_summary.tsv")
    print(summary.to_string(index=False)); print(summary_pred.to_string(index=False))


if __name__ == "__main__":
    main()
