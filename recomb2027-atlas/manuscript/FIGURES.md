# Figure sources

The figures use white backgrounds, open axes, restrained colors and explicit
symbols. They are drawn with Matplotlib using a custom style; they were not
generated in GraphPad Prism. The retained PDFs preserve vector lines and text; the scripts also produce SVG and PNG exports.

From this directory, regenerate the figures with Python, NumPy and Matplotlib:

```sh
python scripts/geometry_pipeline.py
python scripts/journal_figures.py
```

The scripts read the included `figure_data/` files. They do not rerun the
experiments. `scripts/prepare_figure_data.py` documents how these inputs were
extracted from the restored scientific archive and public extension records.
`figure_data/provenance.json` records the source hashes and scientific revision.

| Figure | Input | Interpretation |
| --- | --- | --- |
| Main 1: boundary, envelopes and count states | `scripts/geometry_pipeline.py` | Analytic schematic, not a tissue observation. |
| Main 2: LDHB geometry and intervals | `ldhb_geometry.json`, `ldhb_mask.npz`, `ldhb_intervals.csv` | Exact step locations; all eight archived budget rows agree. |
| Main 3: donor replication | `donor_regions_plot.csv`, `donor_budget_plot.csv`, `donors.csv` | Regional points and pooled donor fractions. Eligibility is fixed over the complete grid. |
| Main 4: coverage and reference-count reach | `width_summary.csv`, `union_reach.csv` | Fixed reference denominators; the displayed union comparison has no area cap. |
| Supplement: width and budget | `width_budget_plot.csv` | Each atlas has its own positive-source or editable pair population. |
| Supplement: count-chain runtime | `chain_scaling_plot.csv` | Three processes per constructed size; open points and median. |
| Supplement: Proseg context | `proseg_context_plot.csv` | Paired input windows for the same core molecules, one no-factorization run per setting. |

To edit a graph in Prism, import the relevant CSV into an XY or column table,
using the unit-labelled columns. Use the provided percentages directly for donor
budget curves: their denominator is the full-grid eligible population, not the
number eligible at an individual budget. Pooled donor fractions divide summed
counts; they are not means of the three regional percentages.

The original data and exact fractions remain in the accompanying JSON/CSV files.
Lines between sampled budgets and between paired runs are visual guides.
Only the LDHB interval curves use exact step thresholds. No error bars or
significance tests have been added to unreplicated measurements.

The native-mask and donor-location maps retain the original source coordinates.
The full analysis and its data rights are described in the paper and supplement.

## Final-size donor map

`python scripts/figure_typography.py` reletters the saved donor map with 10-point
labels at a width of 6.5 inches. It requires PyMuPDF and the Nimbus Sans fonts
from URW Base35. The input `figure_data/biological_regions_original.pdf` is the
original coordinate map; its points, region boxes and axes are preserved. Axis
labels use millimetres instead of micrometres with a 10³ multiplier. The caption
identifies the 200 µm cores and 800 µm context windows. This typography pass uses
no new tissue data. The raw-data plotting script remains available separately.

The final PDF audit measures ordinary figure text after LaTeX scaling, not just
the plot's native font setting. Every ordinary label is at least 10 PDF points;
mathematical subscripts and superscripts use conventional smaller type.

After regenerating plots, run `python scripts/normalize_embedded_fonts.py`
(PyMuPDF and fontTools). This corrects older plotting-backend output that embeds
OpenType/CFF data under a TrueType FontFile2 declaration. It converts those font
programs to TrueType while preserving glyph order, advance widths and PDF text
coordinates. Diagram geometry and numerical data are unchanged. The final PDFs
have embedded fonts and produce no Poppler font-embedding warnings.
