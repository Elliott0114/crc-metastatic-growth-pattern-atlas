"""Parse deposited GEO array data and clinical records without redefining endpoints."""
import csv
import gzip
import json
import numpy as np
import pandas as pd
from rec_public_upgrade_common import DATA, OUT, META, soft_header

dest=OUT / "clinical"; dest.mkdir(parents=True, exist_ok=True)
with gzip.open(DATA / "GSE159216_series_matrix.txt.gz", "rt", encoding="utf-8") as f:
    header=soft_header(f)
    expression=pd.read_csv(f, sep="\t", index_col=0, comment="!")
samples=header["!Sample_geo_accession"][0]
assert samples == expression.columns.tolist()
metadata=pd.DataFrame({"sample_id":samples,"title":header["!Sample_title"][0]})
for row in header["!Sample_characteristics_ch1"]:
    parts=[x.split(": ",1) for x in row]
    assert len({x[0] for x in parts})==1
    metadata[parts[0][0]]=[x[1] if len(x)>1 else "" for x in parts]
metadata.to_csv(dest / "geo_sample_metadata.tsv", sep="\t", index=False)

maps=[]
with (DATA / "GPL17586.soft.txt").open(encoding="utf-8") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):break
    for row in csv.DictReader(f, delimiter="\t"):
        if row["ID"].startswith("!platform_table_end"):break
        if row["ID"] not in expression.index:continue
        symbols=set()
        for a in row.get("gene_assignment", "").split(" /// "):
            fields=a.split(" // ")
            if len(fields)>1 and fields[1].strip() not in ("---","NULL",""):
                symbols.add(fields[1].strip().upper())
        status="unique_gene" if len(symbols)==1 else ("ambiguous_gene" if symbols else "no_symbol")
        maps.append(dict(probe_id=row["ID"],gene=next(iter(symbols)) if len(symbols)==1 else "",
                         status=status,all_symbols=";".join(sorted(symbols))))
mapping=pd.DataFrame(maps)
missing=set(expression.index)-set(mapping.probe_id)
mapping=pd.concat([mapping,pd.DataFrame([dict(probe_id=x,gene="",status="no_platform_annotation",all_symbols="") for x in sorted(missing)])],ignore_index=True)
mapping.to_csv(dest / "probe_mapping.tsv",sep="\t",index=False)
unique=mapping.query("status == 'unique_gene'").set_index("probe_id")
x=expression.loc[unique.index].copy();x.index=unique.gene
x=x.groupby(level=0,sort=True).median()
assert x.index.is_unique and np.isfinite(x.to_numpy()).all()
x.index.name="gene";x.to_csv(dest / "gene_expression_log2.tsv.gz",sep="\t",compression="gzip")
coverage={"n_array_features":len(expression),"n_unique_mapped_probes":len(unique),"n_genes":len(x),"probe_status":mapping.status.value_counts().to_dict(),"processing":"Deposited TAC/RMA log2 values; no second log transform, no count model; ambiguous assignments removed; unique probes median collapsed. This differs from the source paper's custom Brainarray preprocessing."}
(dest / "matrix_qc.json").write_text(json.dumps(coverage,indent=2)+"\n",encoding="utf-8")
print(json.dumps(coverage),flush=True)
