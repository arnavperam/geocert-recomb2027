# Scientific source archive

Distribution 1.1.1 contains the scientific implementation, prepared inputs,
study protocols and numerical results. The archive includes the files needed to
replay the base experiments.

| Archive | Value |
|---|---|
| Filename | `GeoCert_RECOMB_1.1.1_Source.zip` |
| Size | 74,385,746 bytes |
| Parts | 18 |
| SHA-256 | `a14c8986551a4c568a822991981198db88bb583b22b23329a5ed09c844904f03` |

[artifact_parts/manifest.json](artifact_parts/manifest.json) specifies assembly order,
part hashes and total size. [ARCHIVE_SHA256SUMS.txt](ARCHIVE_SHA256SUMS.txt) identifies
the individual files inside the archive. All parts are included here.

## Extract and verify

Run from the repository root using Python 3.12 or 3.13:

```sh
python recomb2027-atlas/reassemble.py --extract atlas-replay
cd atlas-replay/recomb2027-atlas
python atlas_cli.py verify-manifest
python release/restore_results.py
python release/check_interval_hashes.py
python atlas_cli.py validate
python atlas_cli.py showcase --output run/showcase.json
```

Keep assertions enabled. These commands use the standard library and require no
external repository or new data download. Allow about 1.2 GB for restored records.
Extraction verifies part and ZIP integrity and requires an absent or empty output
directory. Result restoration checks compressed and expanded hashes.

The base archive contains the primary crops and width study. The [extensions](extensions/)
contain the atlas union, donor replication, contiguous-region experiments and
composition checks, with their own source records and replay commands.
[Current papers](manuscript/) and their [release checksums](release/PUBLICATION.json)
are maintained separately from the base archive.

The [scientific workflow](../.github/workflows/atlas-final.yml) verifies the public
code against the archive and replays the scientific checks on Python 3.12 and 3.13.
It does not rerun Proseg acquisition or replace recorded historical timings.

## Repackage an edited source tree

Use [package_archive.py](../scripts/package_archive.py) after assigning a new
version in the expanded package's `VERSION` file:

```sh
python scripts/package_archive.py path/to/recomb2027-atlas --parts recomb2027-atlas/artifact_parts
```

The script regenerates the internal file manifest, archive, parts and checksums.
Do this before restoring compressed results or generating bytecode in that tree.
The repository manifest and publication record must then be refreshed for the new
archive identity.
