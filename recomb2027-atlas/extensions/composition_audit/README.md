# Complete-atlas composition

A tissue-wide deformation budget includes every patch of the chosen atlas.
A patch that cannot change the queried count can still require positive area.
The composition must retain that cost and a boundary attaining it.

`complete.py` applies the existing min-plus routine to a declared complete
atlas. It requires a complete spectrum for every patch, including zero-score
patches, and returns all selected local witnesses. Missing or unfinished
spectra are errors. An infeasible local spectrum makes the global family empty.
Local spectra must use count changes relative to the same source, physical area
in the same units, and witnesses verified by the local solver. Geometry and
local optimality remain the responsibility of that solver.

The scalar count chain already supports positive mandatory area. In contrast,
the original tissue adapter uses a zero-cost source contract for its closure
fast path. The saved tissue inputs satisfy that contract, including RNA-free
patches skipped by the count loop. This extension leaves the core algorithms
and recorded tissue results unchanged.

## Checks

Run from the repository root:

```sh
python recomb2027-atlas/extensions/composition_audit/test_composition.py
python recomb2027-atlas/extensions/composition_audit/audit_sources.py --data-root /path/to/restored/recomb2027-atlas
```

The source audit uses the base scientific archive restored by `reassemble.py`
and `release/restore_results.py`, plus the external-validation and scaling
archives already present in this repository. It verifies archive hashes and
atlas identities before directly testing the zero-area source boundary using
rational inequalities. It imports no GeoCert optimization or geometry code.

The regression constructs two patches: transferring the target RNA costs 9/16,
while a second patch with no target contribution requires area 1/16. At budget
3/5 the exact interval is [1, 1]; omitting the second patch would incorrectly
give [0, 1]. Analytic triangle areas provide an independent expected answer.
The test also checks empty feasible sets, fixed counts with positive minimum
area, complete witnesses, and missing or unfinished inputs.

Recorded outputs are in `results/`. All 28 distinct saved atlases pass the
source check: 26,349 patches, including 4,793 RNA-free patches. This establishes
that the clarified composition condition changes none of the reported spectra.
It is a check of scalar composition and its data assumptions, not an independent
proof review of the coupled-interface theorem or an anatomical calibration.
