"""Member-level evidence tables and figures for the focused evidence bridge."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from rec_junction_context_common import ROOT, bh, definitions, digest, mean_summary
from analyze_rec_evidence_bridge_spatial import OUT, group_summary, write

PREVIOUS = ROOT / "analysis_results/rec_junction_context_2026-09-16"
COLORS = {"rHGP": "#B35806", "dHGP": "#2166AC", "Claudin": "#0072B2",
          "Polarity": "#D55E00", "HRC": "#009E73", "Junction": "#7755A5"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    (OUT / "figures").mkdir(exist_ok=True)
    fig.savefig(OUT / f"figures/{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"figures/{name}.png", dpi=320, bbox_inches="tight")
    plt.close(fig)



B=OUT
C=PREVIOUS
N=ROOT / "analysis_results/rec_manuscript_reinforcement_2026-09-17"
def read(base, name):
    return pd.read_csv(base/name, sep="\t")

def main():
    sets = definitions()
    genes = sorted(sets["Claudin"], key=lambda g: int(g[4:])) + sets["Polarity"] + ["F11R"]
    all_gene = pd.read_csv(OUT / "counts/junction_gene_results.tsv", sep="\t")
    nb = all_gene[all_gene.normalization == "TMM"]
    members = pd.DataFrame(index=pd.Index(genes, name="gene"))
    members["component"] = ["Claudin" if g in sets["Claudin"] else "Polarity" if g in sets["Polarity"] else "F11R" for g in genes]
    ogden = pd.read_csv(PREVIOUS / "cells/gene_summary.tsv", sep="\t")
    ogden = ogden[ogden.scope == "All_epithelial"].set_index("gene")
    spatial = pd.read_csv(OUT / "spatial/correlation_summary.tsv", sep="\t")
    sc = spatial[(spatial.restriction == "all_tumour") & (spatial.model == "adjusted")
                 & spatial.component.str.startswith("gene_")].copy()
    sc["gene"] = sc.component.str.removeprefix("gene_")
    sc = sc.set_index("gene")
    members["Ogden_HRC_rho"] = ogden.effect
    members["Spatial_HRC_rho"] = sc.effect
    for context, d in nb.groupby("dataset"):
        for metric in ["logFC", "PValue", "q_junction"]:
            members[context + "_" + metric] = d.set_index("gene")[metric]
    origins = pd.read_csv(PREVIOUS / "cells/compartment_detection_summary.tsv", sep="\t")
    for compartment in ["Epithelial", "Hepatocyte", "Endothelial", "Stromal"]:
        d = origins[(origins.compartment == compartment) & (origins.n_patients >= 3)].set_index("gene")
        members[compartment + "_detection"] = d.mean_detection
    write(members.reset_index(), "integration/all_junction_members.tsv")

    members=read(B,'integration/all_junction_members.tsv').set_index('gene')
    section=read(N,'sections/gene_results.tsv.gz');section=section[section.normalization=='TMM'].set_index('gene')
    members['Same_section_macro_minus_micro_logFC']=section.logFC
    genes=members.index.tolist()
    eligible=read(B,'counts/gene_eligibility.tsv.gz')
    lock=read(C,'cells/gene_control_lock.tsv');measured=set(lock[(lock.role=='target')&lock.source_gene.notna()].gene)
    spcov=read(B,'spatial/coverage.tsv');spmeasured=set(spcov[spcov.measured].gene)
    se=read(N,'sections/gene_eligibility.tsv')
    covcols=['Epithelial_detection','Hepatocyte_detection','Endothelial_detection','Stromal_detection']
    rhocols=['Ogden_HRC_rho','Spatial_HRC_rho']
    contexts=['Bulk_rHGP_minus_dHGP','Spatial_rHGP_minus_dHGP_tumour','Spatial_rHGP_minus_dHGP_near','Liu_macro_minus_micro','Same_section_macro_minus_micro']
    fccols=[k+'_logFC' for k in contexts]
    status=[]
    for gene in genes:
        for col in rhocols+fccols:
            if pd.notna(members.loc[gene,col]):label='estimated'
            elif col in rhocols:
                assay=measured if col==rhocols[0] else spmeasured
                label='insufficient_detection_or_variation' if gene in assay else 'not_measured'
            else:
                ctx=col.removesuffix('_logFC');known=se if ctx.startswith('Same_section') else eligible[eligible.dataset==ctx]
                label='abundance_ineligible' if gene in set(known.gene) else 'not_measured'
            status.append({'gene':gene,'column':col,'status':label,'value':members.loc[gene,col]})
    st=pd.DataFrame(status)
    dest=N/'integration';dest.mkdir(exist_ok=True)
    members.reset_index().to_csv(dest/'all_junction_members.tsv',sep='\t',index=False)
    st.to_csv(dest/'member_measurement_status.tsv',sep='\t',index=False)

if __name__ == "__main__":
    main()
