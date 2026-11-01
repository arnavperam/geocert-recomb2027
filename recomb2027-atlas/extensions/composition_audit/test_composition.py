"""Regression cases for mandatory area in patches with zero target score.

Local spectra come from the unchanged critical-height solver. Expected costs
come from triangular areas calculated below, independently of the count DP.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from atlas.chain_spectrum import compile_chain, spectrum
from atlas.model import interval, make_patch, verify_witness
from mixed import MixedLayers, TaggedBox, verify_mixed
from complete import compose_complete


def local(x, score, source=0, pins=None):
    problem = MixedLayers(
        [[(0, 0), (4, 0)]],
        [TaggedBox(x, x, 2, 2, 0, False)],
        [source], [1], [], [], clearance=F(1, 4),
        bands=[(-1, 1)], mode="nominal",
    )
    answer = spectrum(compile_chain(problem, pins), [score])
    reference_count = score[source]
    result = {
        "complete": answer["complete"],
        "minima": {q - reference_count: a for q, a in answer["minima"].items()},
        "witnesses": {q - reference_count: w for q, w in answer["witnesses"].items()},
    }
    for w in result["witnesses"].values():
        verify_mixed(problem, w["owners"], w["boundaries"], w["area"])
    return result


def expect_rejected(action):
    try:
        action()
    except (ValueError, KeyError):
        return
    raise AssertionError("invalid composition input was accepted")


def main(output):
    # Flipping this target needs a downward tent of height 3/4, with area 9/16.
    target = local(F(-1, 2), [1, 0])
    assert target["minima"] == {0: F(0), -1: F(9, 16)}
    # At x=0, either owner needs a tent of height 1/4, hence area 1/16.
    # This patch has no target-gene contribution, but cannot be left at zero.
    mandatory = local(F(0), [0, 0])
    assert mandatory["minima"] == {0: F(1, 16)}
    full = compose_complete(["target", "other"], {"target": target, "other": mandatory}, 1)
    assert {q: r["area"] for q, r in full.items()} == {0: F(5, 8), 1: F(1, 16)}
    for q, result in full.items():
        assert sum(F(w["area"]) for w in result["witnesses"].values()) == result["area"]
        assert 1 + sum(result["changes"].values()) == q
    states = {q: (v["area"], v["changes"]) for q, v in full.items()}
    for budget in [F(0), F(1, 32)]:
        assert interval(states, budget)["status"] == "empty_uncertainty_set"
    for budget in [F(1, 16), F(3, 5)]:
        got = interval(states, budget)
        assert (got["lower"], got["upper"]) == (1, 1)
    got = interval(states, F(5, 8))
    assert (got["lower"], got["upper"]) == (0, 1)

    # Omitting the other patch would falsely admit a removal at budget 3/5.
    dropped = {1 + d: (a, {}) for d, a in target["minima"].items()}
    wrong = interval(dropped, F(3, 5))
    assert (wrong["lower"], wrong["upper"]) == (0, 1)
    # Completeness is explicit, so this erroneous omission fails at the API.
    expect_rejected(lambda: compose_complete(["target", "other"], {"target": target}, 1))
    expect_rejected(lambda: compose_complete(["target", "target"], {"target": target}, 1))
    unfinished = dict(mandatory, complete=False)
    expect_rejected(lambda: compose_complete(["other"], {"other": unfinished}, 1))

    # A positive zero-score cost also matters when the count is fixed everywhere.
    fixed = compose_complete(["other"], {"other": mandatory}, 37)
    assert set(fixed) == {37} and fixed[37]["area"] == F(1, 16)
    # A zero-cost source is exactly the special case used by the tissue pipeline.
    free = local(F(-1, 2), [0, 0])
    simple = compose_complete(["target", "other"], {"target": target, "other": free}, 1)
    assert {q: r["area"] for q, r in simple.items()} == {0: F(9, 16), 1: F(0)}
    # A pinned owner outside the band makes even an unscored patch infeasible.
    empty = local(F(2), [0, 0], source=1, pins={0: 0})
    assert empty["minima"] == {}
    assert compose_complete(["target", "other"], {"target": target, "other": empty}, 1) == {}
    assert compose_complete([], {}, 37)[37]["area"] == 0

    # The existing tissue adapter explicitly rejects the source-infeasible point.
    patch = {
        "points": [{"gene": "OTHER", "source": 0, "x": "0", "t": "2"}],
        "nuclear_constraints": [], "clearance_pixels": "1/4",
        "domain_pixels": [0, 4], "band_halfwidth_pixels": "1",
        "halfwidth_pixels": "2", "slope": "1", "area_scale_um2": "1",
        "cell_labels": [3, 4],
    }
    expect_rejected(lambda: make_patch(patch))
    expect_rejected(lambda: verify_witness(patch, [0], [[(0, 4, 0, 0)]], 0))
    result = {
        "scope": "Complete scalar-atlas composition and source-feasible adapter contract",
        "target_change_cost": "9/16", "unscored_patch_minimum_area": "1/16",
        "complete_minima": {str(q): str(r["area"]) for q, r in sorted(full.items())},
        "budget_3_5_exact_interval": [1, 1], "budget_3_5_if_patch_dropped": [0, 1],
        "tested_budgets": ["0", "1/32", "1/16", "3/5", "5/8"],
        "zero_target_positive_area": True, "infeasible_unscored_patch": True,
        "missing_or_incomplete_patch_rejected": True, "zero_cost_source_case": True,
        "complete_attaining_witnesses": True, "tissue_adapter_rejects_inadmissible_source": True,
        "failures": 0,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O: regression assertions must be enabled.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results/regression.json")
    main(parser.parse_args().output)
