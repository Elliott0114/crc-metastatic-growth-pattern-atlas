"""Reproduce the pre-KRT20-correction spatial column retained in Table 4.

This is a historical compatibility branch, not the corrected primary analysis.
The original correction is documented in the 24 August spatial specification.
"""
import json
import pandas as pd

import analyze_rec_program_hgp_spatial_projection as spatial
from analyze_rec_program_cross_modal_concordance import patient_gene_effects


def main():
    spatial.MARKERS['epithelial'].append('KRT20')
    program = pd.read_csv(spatial.ROOT / 'analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv', sep='\t')
    metadata = spatial.load_metadata()
    rows = []
    for sample in metadata['sample']:
        _, _, values, _ = spatial.load_sample(sample, 200, program.gene.tolist())
        rows.extend(values)
    frame = pd.DataFrame(rows).merge(metadata, on='sample', validate='many_to_one')
    frame = frame[(frame.restriction == 'all_tumour_side') & (frame.region == 'near_interface_0_500um')]
    out = spatial.ROOT / 'analysis_results/retained_rec44_historical'
    out.mkdir(parents=True, exist_ok=True)
    patient_gene_effects(frame, 'spatial').to_csv(out / 'spatial_gene_effects.tsv', sep='\t', index=False, encoding='utf-8')
    (out / 'definition.json').write_text(json.dumps(dict(
        role='Historical spatial column in Supplementary Table 4 / REC44_cross_context only',
        epithelial_markers=spatial.MARKERS['epithelial'],
        qc_min_detected_genes=200, boundary_hops=5,
        caveat='KRT20 overlaps REC50. The corrected primary branch excludes it.',
    ), indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
