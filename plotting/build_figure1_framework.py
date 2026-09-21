"""Three resource-to-analysis tracks above the unchanged Figure 1 estimates."""
from pathlib import Path
import hashlib
import json

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd

import figure_common as fc
import build_main_figures as original
import build_review as review

OUT = Path(__file__).resolve().parents[1]
ART = OUT / 'assets/medical_icons.png'
INK = '#20354B'
for module in [fc, original, review]:
    module.OUT = OUT
fc.MANIFEST = OUT / 'panel_manifest.json'
fc.RECORDS = json.loads(fc.MANIFEST.read_text(encoding='utf-8'))

TITLE = 'Study framework, spatial architecture and NC3 reference mapping'
CAPTION = (
    'a, Arrows link resources to analyses; the dashed Ogden link supplies the state reference. AI-generated schematic icons: OpenAI imagegen. '
    'Patient counts refer to analysis-specific subsets. Ogden’s 14 denotes component discovery. '
    'Atlas is an acquisition container; E-MTAB-12022/12043 share four patients; '
    'ten same-section patients are among eleven with paired regions. NC3, source population; REC, reference state; '
    'HRC, primary regenerative readout; REC50, auxiliary projection. '
    'b, P02 (dHGP) and P18 (rHGP) ISS fields share a coordinate scale; the overview locates the P18 crop. '
    'Analyses use complete fields. '
    'c,d, NC3 fraction and the fraction with ≥1 liver epithelial neighbour (ten-neighbour graph), both among neoplastic cells. '
    'Points: 7 dHGP, 8 rHGP patients; bars: means, with equal within-patient ROI weights. '
    'e, Observed-minus-conditional-null adjacency and patient-bootstrap 95% intervals. Homotypic adjacency requires '
    'another NC3 neighbour; dual adjacency requires both NC3 and damaged-hepatocyte (DHC) neighbours around an NC3 cell. '
    'The null preserves local abundance, anatomical strata and coordinates. pp, percentage points. '
    'f, ISS × Ogden centroid Pearson correlations across 94 genes. '
    'NC3–REC (r = 0.387) is outlined. Colono., Colonocyte; Interm., Intermediate; Stem-N., Stem NOTUM. '
    'Robustness: Supplementary Fig. 1 and Table 2.'
)
review.LEGENDS = {'Figure_1': (TITLE, CAPTION)}


def icon_bounds(image):
    """Identify six view windows; the generated artwork pixels are not edited."""
    h, w = image.shape[:2]
    records = []
    for row in range(2):
        for col in range(3):
            x0, x1 = round(col*w/3), round((col+1)*w/3)
            y0, y1 = round(row*h/2), round((row+1)*h/2)
            visible = np.min(image[y0:y1, x0:x1, :3], axis=2) < .94
            if image.shape[2] == 4:
                visible &= image[y0:y1, x0:x1, 3] > .05
            yy, xx = np.where(visible)
            assert len(xx), (row, col)
            records.append([max(x0, x0+int(xx.min())-8), max(y0, y0+int(yy.min())-8),
                            min(x1, x0+int(xx.max())+9), min(y1, y0+int(yy.max())+9)])
    return records


def framework(fig):
    records = []
    height = fig._height_mm
    picture = plt.imread(ART)
    windows = icon_bounds(picture)

    def label(x, y, value, bold=False, size=8, color=INK, **kwargs):
        obj = fc.text(fig, x, y, value, fontsize=size, color=color,
                      fontweight='bold' if bold else 'normal', **kwargs)
        records.append({'text': value, 'font_pt': size, 'x_mm': x, 'top_mm': y})
        return obj

    def box(x, y, w, h, fill, edge='none'):
        fig.patches.append(FancyBboxPatch((x/170, 1-(y+h)/height), w/170, h/height,
            boxstyle='round,pad=0,rounding_size=0.004', transform=fig.transFigure,
            facecolor=fill, edgecolor=edge, linewidth=.5, zorder=-2))

    def path(points, dashed=False, arrow=False):
        xy = [(x/170, 1-y/height) for x, y in points]
        p = MplPath(xy, [MplPath.MOVETO] + [MplPath.LINETO]*(len(xy)-1))
        if arrow:
            obj = FancyArrowPatch(path=p, arrowstyle='-|>', mutation_scale=7,
                linewidth=.5 if dashed else .7, color='#889AA8' if dashed else '#61788A',
                linestyle=(0, (2.4, 2)) if dashed else '-',
                transform=fig.transFigure, zorder=1)
        else:
            obj = PathPatch(p, fill=False, linewidth=.6, color='#D1DCE2',
                            transform=fig.transFigure, zorder=-1)
        fig.add_artist(obj)

    def icon(index, x, y, w, h):
        a = fc.ax(fig, x, y, w, h)
        a.imshow(picture, interpolation='none')
        x0, y0, x1, y1 = windows[index]
        a.set(xlim=(x0, x1), ylim=(y1, y0))
        a.set_aspect('equal'); a.set_axis_off()
        a.set_gid('schematic_icon_' + str(index))

    label(5, 1, 'A', True, 11)
    label(11, 1.2, 'Multi-cohort study design', True, 9.5)
    label(5, 7, 'Data resources', True, color='#536A7B')
    label(106, 7, 'Patient-level analyses', True, color='#536A7B')

    # Leave dataset labels light; only the three resource headings carry emphasis.
    label(7, 11.8, 'Spatial pathology context', True, 8.5)
    label(7, 26, 'Epithelial transcription', True, 8.5)
    label(7, 40.8, 'Tissue context', True, 8.5)
    path([(5, 25), (101, 25)])
    path([(5, 39.8), (101, 39.8)])
    secondary = '#526577'
    icon(0, 7.5, 16, 17, 8.2)
    label(28, 16.2, 'ISS', size=8.5)
    label(28, 20.2, '17 patients', color=secondary)
    # Identical three-level baselines make discovery and validation roles explicit.
    role_color = '#3E766D'
    for index, x_icon, x_label, role, title, count in [
        (1, 7, 22, 'Discovery', 'Ogden scRNA-seq', '14 patients'),
        (2, 54.8, 68.2, 'Independent validation', 'External scRNA-seq', '5 studies · 28 patients')]:
        icon(index, x_icon, 31.2, 12, 7.5)
        label(x_label, 29.5, role, True, 7.2, color=role_color)
        label(x_label, 32.8, title)
        label(x_label, 36.4, count, color=secondary)
    # Three spacious source tiles, with artwork above rather than squeezed beside text.
    for index, center, title, count in [
        (3, 21, 'Bulk RNA-seq', '15 patients'),
        (4, 53, 'HGP spatial RNA', '6 patients'),
        (5, 85, 'Micro/macro', '11 patients ·\npaired regions')]:
        icon(index, center-7, 44.2, 14, 5.5)
        label(center, 50, title, ha='center')
        label(center, 53.5, count, ha='center', color=secondary, linespacing=1.05)

    for y, h, fill in [(11, 13.8, '#F3F6F9'), (26, 14.5, '#EDF6F3'), (42, 15.5, '#F6F5F8')]:
        box(106, y, 59, h, fill)
    for y, title, line1, line2 in [
        (13.1, 'NC3 spatial context', 'NC3–REC mapping', 'Conditional adjacency'),
        (28, 'Components & members', 'Within-patient comparison', 'Meta-analysis · 30 junction genes'),
        (44.3, 'Tissue contrasts', 'HGP · Micro/macro regions', 'Source-proxy adjustment')]:
        label(109, y, title, True, 8.5)
        label(109, y+4.1, line1)
        label(109, y+7.9, line2)
    for y in [15, 33.4, 49.8]:
        path([(101.2, y), (105.7, y)], arrow=True)
    # Ogden has a separate reference role; this route avoids both resource labels.
    path([(48, 34.1), (54.1, 34.1), (54.1, 20.7), (103.4, 20.7), (103.4, 18.8), (105.7, 18.8)],
         dashed=True, arrow=True)
    label(76, 17.2, 'State reference', size=6.8, ha='center', color=secondary)
    box(5, 60, 160, 3.5, '#F8FAFB')
    label(8, 60.55, 'Supporting analyses: protein · clinical · perturbation', size=6.6, color=secondary)
    label(94, 60.5, 'Primary readouts: HRC · Claudin · Polarity', size=6.8, color=INK)
    path([(88, 60.55), (88, 62.8)])
    return records, windows


def assemble(fig, name):
    # Frozen implementation constructs the empirical artists and source tables.
    # Only positions, display labels and the schematic are changed here.
    for patch in fig.patches:
        patch.set_visible(False)
    for obj in list(fig.texts):
        if (1-obj.get_position()[1])*175 < 37:
            obj.remove()
    empirical = [a for a in fig.axes if a.get_position(original=True).height*175 < 174]
    def snapshot():
        def digest(values):
            array = np.ma.asarray(values, dtype=float)
            return {'shape': list(array.shape), 'sha256': hashlib.sha256(
                array.filled(np.nan).tobytes()+np.ma.getmaskarray(array).tobytes()).hexdigest()}
        return [{'xlim': list(a.get_xlim()), 'ylim': list(a.get_ylim()),
                 'lines': [(digest(line.get_xdata()), digest(line.get_ydata())) for line in a.lines],
                 'collections': [{'offsets': digest(c.get_offsets()),
                                  'values': None if c.get_array() is None else digest(c.get_array())}
                                 for c in a.collections],
                 'rectangles': [(list(p.get_xy()), p.get_width(), p.get_height())
                                for p in a.patches if isinstance(p, fc.Rectangle)]}
                for a in empirical]
    empirical_before = snapshot()
    for axis in list(fig.axes):
        position = axis.get_position(original=True)
        x, top = position.x0*170, (1-position.y1)*175
        if position.height*175 > 174:
            axis.remove(); continue
        if round(top) == 53:
            xx, ww = (5, 65) if x < 50 else (85, 54)
            yy, hh = 75, 26
        elif round(top) == 54:
            xx, yy, ww, hh = 139, 75, 12, 26
        elif round(top) == 101:
            xx, yy, ww, hh = x, 114, position.width*170, 19.5
        elif round(top) == 143:
            xx, yy, ww, hh = x, 147, position.width*170, 17
        elif round(top) == 138:
            xx, yy, ww, hh = x, 143, position.width*170, 1.5
        else:
            raise AssertionError((x, top))
        axis.set_position([xx/170, 1-(yy+hh)/175, ww/170, hh/175])
        if axis.get_ylabel() == 'Liver neighbour (%)':
            axis.set_ylabel('Liver neighbour\n(%)')
    moves = {41:66, 93:108, 133:140, 138:143}
    for obj in fig.texts:
        value = obj.get_text(); x, y = obj.get_position(); old = (1-y)*175
        if value == 'Original ISS coordinate fields':
            obj.set_text('Original ISS fields')
        if value == 'Pearson r':
            obj.set_position((120/170, 1-143.8/175)); obj.set_va('center'); continue
        if value in ['P02 · dHGP', 'P18 · rHGP detail', 'Overview']:
            xx = {'P02 · dHGP':37.5, 'P18 · rHGP detail':112, 'Overview':145}[value]
            obj.set_position((xx/170, 1-71/175)); obj.set_ha('center'); continue
        if value in ['NC3–NC3', 'NC3 + DHC']:
            yy = 151.25 if value == 'NC3–NC3' else 159.75
            obj.set_text('NC3\nhomotypic' if value == 'NC3–NC3' else 'NC3 + DHC\n(dual)')
            obj.set_position((2/170, 1-yy/175)); obj.set_ha('left'); continue
        for before, after in moves.items():
            if abs(old-before) < .01:
                obj.set_position((x, 1-after/175)); break
            if abs(old-before-.2) < .01:
                obj.set_position((x, 1-(after+.2)/175)); break
    fig.legends[0].set_bbox_to_anchor((.03, 1-102/175), transform=fig.transFigure)
    records, windows = framework(fig)
    assert empirical_before == snapshot(), 'Empirical artists changed during layout.'
    (OUT/'review/empirical_artist_checks.json').write_text(json.dumps(
        {'status': 'passed', 'axes': empirical_before, 'unchanged': 'coordinates, statistical artists, matrix, crop and limits'},
        indent=2)+'\n', encoding='utf-8')
    (OUT/'review/Figure_1A_typography.json').write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8')
    (OUT/'assets/icon_view_windows.json').write_text(json.dumps(windows, indent=2)+'\n', encoding='utf-8')
    fc.save(fig, name)


def main():
    fig = fc.canvas(64)
    framework(fig)
    (OUT/'study_design').mkdir(exist_ok=True)
    for extension in ['pdf', 'svg', 'png']:
        fig.savefig(OUT/'study_design'/f'Figure_1A_study_design.{extension}', dpi=600)
    plt.close(fig)
    original.save = assemble
    original.figure1()
    dimensions_path = OUT/'review/figure_dimensions.json'
    dimensions = json.loads(dimensions_path.read_text(encoding='utf-8'))
    current_dimensions = next(r for r in dimensions if r['figure'] == 'Figure_1')
    current_dimensions['figure_body_export'] = 'vector text, statistical lines, heatmap and schematic connectors; biomedical icons and dense ISS points are raster layers'
    current_dimensions['standalone_panel_A_height_mm'] = 64
    dimensions_path.write_text(json.dumps(dimensions, indent=2)+'\n', encoding='utf-8')
    hashes = json.loads((OUT/'baseline_panel_hashes.json').read_text(encoding='utf-8'))
    mode = json.loads((OUT/'validation_mode.json').read_text(encoding='utf-8'))['mode']
    for name, digest in hashes.items():
        if mode == 'quick':
            assert fc.sha(OUT/'source_data/panels'/name) == digest, name
        else:
            pd.testing.assert_frame_equal(
                pd.read_csv(OUT/'reference_panels'/name, sep='\t'),
                pd.read_csv(OUT/'source_data/panels'/name, sep='\t'),
                check_dtype=False, check_exact=False, atol=1e-8, rtol=1e-6)
    (OUT/'Figure_1_legend.md').write_text(TITLE+'\n\n'+CAPTION+'\n', encoding='utf-8')
    links = [
        ('ISS', 'NC3 spatial context', 'solid', 'NC3 abundance, spatial adjacency and source state for mapping'),
        ('Ogden epithelial states', 'NC3–REC mapping', 'dashed', 'State reference; 14 refers to component discovery, not the reference coverage'),
        ('Ogden component discovery; five external studies', 'Components & members', 'solid', 'Within-patient contrasts, independent meta-analysis and all 30 members'),
        ('Bulk HGP; spatial HGP; paired micro/macro regions', 'Tissue contrasts', 'solid', 'HGP and regional expression; source-proxy adjustment'),
    ]
    pd.DataFrame(links, columns=['input', 'analysis', 'line_style', 'meaning']).to_csv(
        OUT/'Study_design_connections.tsv', sep='\t', index=False, encoding='utf-8')
    manifest = json.loads(fc.MANIFEST.read_text(encoding='utf-8'))
    for row in manifest:
        if row['figure_id'] != 'Figure_1':
            continue
        row['panel_data_sha256'] = fc.sha(OUT/row['panel_source_data'])
        row['output_sha256'] = fc.sha(OUT/row['output_path'])
        row['legend'] = CAPTION
        if row['panel_id'] == 'A':
            row.update(illustration_artifact='assets/medical_icons.png', illustration_sha256=fc.sha(ART),
                       connections='Study_design_connections.tsv', source_counts='Study_design_sources.tsv',
                       notes='Three resource groups; separate Ogden reference input. Biomedical icon views are raster; all labels, grouping and arrows are vector.')
    fc.MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    pd.DataFrame(manifest).to_csv(OUT/'panel_manifest.tsv', sep='\t', index=False, encoding='utf-8')
    print('Figure 1 framework complete; six panel tables verified.', flush=True)


if __name__ == '__main__':
    main()
