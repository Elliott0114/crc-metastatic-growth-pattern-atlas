"""Retained REC50 patient omissions from the original figure calculation."""
from pathlib import Path
import pandas as pd
N = Path(__file__).resolve().parents[1] / "analysis_results/rec_manuscript_reinforcement_2026-09-17"
def read(base, name):
    return pd.read_csv(base/name, sep="\t")

def main():
    scores=read(N,"bulk/patient_scores.tsv")
    x=read(N,"bulk/logcpm.tsv.gz").set_index("gene");coverage=read(N,"bulk/gene_coverage.tsv")
    # REC50 omission scores are recomputed with the same fixed logCPM matrix.
    rec=coverage[(coverage.component=='REC50')&coverage.measured_variable].gene.tolist()
    rr=scores.set_index('sample_id').r
    recrows=[]
    for omitted in scores.sample_id:
        z=x.drop(columns=omitted);z=z.sub(z.mean(axis=1),axis=0).div(z.std(axis=1),axis=0)
        v=z.loc[rec].mean();recrows.append(v[rr.loc[v.index]==1].mean()-v[rr.loc[v.index]==0].mean())
    pd.DataFrame({"omitted":scores.sample_id,"REC50_restandardized_effect":recrows}).to_csv(N/"bulk/REC50_patient_omissions.tsv",sep="\t",index=False,encoding="utf-8")

if __name__ == "__main__":
    main()
