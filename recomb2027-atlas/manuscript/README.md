# Manuscript sources

From this directory, run each command twice:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
```

The main document has one title page, ten body pages and two reference pages.
The supplement has 28 pages. Figures are included as PDFs; [FIGURES.md](FIGURES.md)
covers regeneration from the retained plotting scripts and numerical inputs.

The tracked PDFs are the verified publication bytes. A rebuild can change PDF
metadata and hashes. [PUBLICATION.json](../release/PUBLICATION.json) records the
published hashes. Both documents link to this repository. The abstract is [ABSTRACT.txt](../ABSTRACT.txt).
