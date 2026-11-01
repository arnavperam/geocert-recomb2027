# GeoCert

Exact expression-count intervals under a shared cell-boundary deformation budget.

Code and data for **GeoCert: Certified Expression Intervals under Segmentation Uncertainty in Spatial Transcriptomics**, by Arnav Peram.

[Manuscript](recomb2027-atlas/manuscript/main.pdf) · [Supplement](recomb2027-atlas/manuscript/supplement.pdf) · [LaTeX and figure sources](recomb2027-atlas/manuscript/) · [Data provenance](recomb2027-atlas/docs/DATA_PROVENANCE.md)

GeoCert computes the lowest and highest attainable cell-by-gene counts within a specified boundary family and returns a boundary witness for each endpoint. The intervals are conditional on that family; they do not estimate anatomical segmentation accuracy.

## Run and verify

The prepared-input solver uses Python's standard library. The scientific workflows test Python 3.12 and 3.13. Keep assertions enabled.

```sh
python scripts/verify_repository.py
python recomb2027-atlas/reassemble.py --extract atlas-replay
cd atlas-replay/recomb2027-atlas
python atlas_cli.py verify-manifest
python release/restore_results.py
python release/check_interval_hashes.py
python atlas_cli.py validate
python atlas_cli.py showcase --output run/showcase.json
```

The example reconstructs the count-17 LDHB witness in breast format-test cell 81. Reassembly checks every archive part and the complete archive. Restoration checks compressed and expanded results without downloading new data. Allow approximately 1.2 GB for expanded base records.

To recompute the primary tables, run from the extracted package:

```sh
python atlas_cli.py reproduce --dataset lung --output-dir run/lung
python atlas_cli.py reproduce --dataset breast --output-dir run/breast
python analysis/reach/reproduce.py --dataset lung --width 1 --output-dir run/lung-width1
```

## Additional experiments

These extensions run from this checkout and use the restored base package where indicated.

| Study | Instructions and records |
|---|---|
| Finite union of complete atlases | [Atlas union](recomb2027-atlas/extensions/atlas_union/README.md) |
| Two breast cancer donors and Proseg context | [Biological replication](recomb2027-atlas/extensions/external_validation/README.md) |
| Contiguous 400 and 600 µm regions | [Scaling study](recomb2027-atlas/extensions/scaling/README.md) |
| Mandatory area and source feasibility | [Composition checks](recomb2027-atlas/extensions/composition_audit/README.md) |

The two [scientific workflows](.github/workflows/) verify archived inputs, algorithm tests, witnesses and recorded experimental outputs. Their logs distinguish fresh replays from stored timing measurements. New data acquisition and new Proseg executions are separate from these checks.

## Paper and reproducibility

See [manuscript build instructions](recomb2027-atlas/manuscript/README.md), [figure regeneration](recomb2027-atlas/manuscript/FIGURES.md), and [archive restoration](recomb2027-atlas/ARCHIVE.md). Optional analysis dependencies are in [requirements-analysis.txt](recomb2027-atlas/requirements-analysis.txt).

The paper, source archive, inputs and results are all distributed here. [PROVENANCE.json](PROVENANCE.json) identifies distribution 1.1.1 and the document and archive checksums. [MANIFEST.sha256](MANIFEST.sha256) covers the checkout.

## License and disclosure

Software is [MIT licensed](LICENSE). Source-derived 10x data retain CC BY 4.0 attribution and the breast format-test restriction. See [NOTICE.md](recomb2027-atlas/NOTICE.md) and [data provenance](recomb2027-atlas/docs/DATA_PROVENANCE.md). Proseg's executable is not redistributed.

Generative AI was used solely for software implementation, testing, and drafting.
