import pandas as pd
from rec_junction_context_common import OUT, write
def read(path):
    return pd.read_csv(OUT/path,sep="\t")

def main():
    origin=read("cells/patient_compartment_genes.tsv")
    origin=origin[origin.n_cells>=20]
    origin_summary=origin.groupby(["gene","compartment"]).agg(mean_detection=("detection_fraction","mean"),n_patients=("patient","nunique")).reset_index()
    write(origin_summary,"cells/compartment_detection_summary.tsv")

if __name__ == "__main__":
    main()
