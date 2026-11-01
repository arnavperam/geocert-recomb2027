"""Compose complete local spectra without assuming a zero-cost source.

The caller declares the complete atlas, including patches that cannot change
the queried count. Each local spectrum is indexed by count change and includes
an attaining witness with an ``area`` in the same physical units as its cost.
Geometry and local optimality are certified by the local solver, not here.
"""
from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from atlas.model import compose


def compose_complete(patch_ids, spectra, source_count):
    """Return global minima and complete tuples of local attaining witnesses.

    A zero-score patch may have a positive cost at change zero. An empty local
    feasible set makes the global set empty; missing or unfinished results are
    rejected rather than treated as zero-cost patches.
    """
    patch_ids = tuple(patch_ids)
    if len(set(patch_ids)) != len(patch_ids):
        raise ValueError("atlas patch identifiers must be unique")
    if set(spectra) != set(patch_ids):
        raise ValueError("one local spectrum is required for every atlas patch")
    if type(source_count) is not int or source_count < 0:
        raise ValueError("source count must be a nonnegative integer")

    parts = []
    witnesses = {}
    for pid in patch_ids:
        local = spectra[pid]
        if local.get("complete") is not True:
            raise ValueError("incomplete local spectrum: " + str(pid))
        minima = {}
        witness = {}
        for key, value in local["minima"].items():
            delta = int(key)
            if str(delta) != str(key) or delta in minima:
                raise ValueError("count changes must be distinct integers")
            area = Fraction(value)
            if area < 0:
                raise ValueError("deformation area must be nonnegative")
            boundary = local["witnesses"][key]
            if Fraction(boundary["area"]) != area:
                raise ValueError("local witness area differs from its minimum")
            minima[delta] = area
            witness[delta] = boundary
        parts.append((pid, minima))
        witnesses[pid] = witness

    states = compose(parts, source_count)
    result = {}
    for count, (area, trace) in states.items():
        if set(trace) != set(patch_ids):
            raise AssertionError("incomplete attaining trace")
        result[count] = {
            "area": area,
            "changes": trace,
            "witnesses": {pid: witnesses[pid][delta] for pid, delta in trace.items()},
        }
    return result
