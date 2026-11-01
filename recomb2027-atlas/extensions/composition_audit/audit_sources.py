"""Check zero-area source witnesses in every reported scalar atlas.

This independent check uses rational inequalities, without importing GeoCert's
compiler, solver, closure or witness routines. It checks the recorded chart
constraints; the existing geometry and source-data audits establish that those
charts describe the input tissue. No count optimization is rerun here.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import gzip
import hashlib
import io
import json
from pathlib import Path
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def check_patch(patch):
    """Verify h(t)=0, with all recorded source owners, including RNA-free patches."""
    left, right = map(F, patch["domain_pixels"])
    band, width = F(patch["band_halfwidth_pixels"]), F(patch["halfwidth_pixels"])
    rho = F(patch["clearance_pixels"])
    require(left == 0 and right > left, "invalid source trace domain")
    require(0 <= band < width and F(patch["slope"]) >= 0, "zero trace violates chart bounds")
    require(rho > 0 and F(patch["area_scale_um2"]) > 0, "invalid clearance or area scale")
    # h=0 is continuous, has slope zero and fixes both endpoints. Its integral
    # against the reference h=0 is zero, including after the physical Jacobian.
    for point in patch["points"]:
        x, t, side = F(point["x"]), F(point["t"]), point["source"]
        require(side in (0, 1) and left < t < right, "invalid RNA owner or coordinate")
        require((x + rho <= 0) if side == 0 else (x - rho >= 0),
                "source owner violates RNA clearance")
    for side, value, a, b in patch["nuclear_constraints"]:
        value, a, b = F(value), F(a), F(b)
        require(side in (0, 1) and left <= a <= b <= right, "invalid structural interval")
        require((value <= 0) if side == 0 else (value >= 0),
                "zero trace violates nuclear ownership")


def audit_atlas(raw, name):
    atlas = json.loads(raw)
    body = {k: v for k, v in atlas.items() if k != "atlas_id"}
    identity = digest(json.dumps(body, sort_keys=True, separators=(",", ":")).encode())
    require(identity == atlas["atlas_id"], "atlas content digest mismatch: " + name)
    seen = set()
    for patch in atlas["patches"]:
        require(patch["id"] not in seen, "duplicate patch identifier")
        seen.add(patch["id"])
        try:
            check_patch(patch)
        except ValueError as error:
            raise ValueError(f"{name}, patch {patch['id']}: {error}") from error
    return {
        "input": name, "input_sha256": digest(raw), "atlas_id": identity,
        "roi": atlas["roi"], "patches": len(seen),
        "rna_free_patches": sum(not p["points"] for p in atlas["patches"]),
        "rna_constraints": sum(len(p["points"]) for p in atlas["patches"]),
        "nuclear_constraints": sum(len(p["nuclear_constraints"]) for p in atlas["patches"]),
        "minimum_required_area_um2": "0", "failures": 0,
    }


def base_atlases(root):
    # Explicit study inventory avoids depending on which compressed files a
    # restoration command has already expanded.
    names = [f"data/atlas_roi_{i:02}.json" for i in (1, 2)]
    names += [f"external/data/atlas_breast_{i:02}.json" for i in (1, 2)]
    for tissue, stem in [("lung", "roi"), ("breast", "breast")]:
        for width in ("r0p5", "r1", "r2", "r3"):
            names += [f"reach/{tissue}/{width}/data/atlas_{stem}_{i:02}.json" for i in (1, 2)]
    for name in names:
        path = root / name
        if path.exists():
            raw = path.read_bytes()
        else:
            raw = gzip.decompress(Path(str(path) + ".gz").read_bytes())
        yield raw, name


def atlas_members(raw, hashes, label):
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = sorted(n for n in archive.namelist()
                       if Path(n).parent.name == "data"
                       and Path(n).name.startswith("atlas_") and n.endswith(".json"))
        require(len(names) == 1, "expected one atlas in " + label)
        for name in names:
            content = archive.read(name)
            require(digest(content) == hashes[name], "atlas file hash mismatch")
            yield content, label + "::" + name


def external_atlases(public_root):
    folder = public_root / "extensions/external_validation/artifacts"
    manifest = json.loads((folder / "manifest.json").read_text())
    require(len(manifest["archives"]) == 8, "external atlas inventory changed")
    for record in manifest["archives"]:
        raw = (folder / record["name"]).read_bytes()
        require(len(raw) == record["bytes"] and digest(raw) == record["sha256"],
                "external archive checksum mismatch")
        yield from atlas_members(raw, manifest["files"], record["name"])


def scaling_atlases(public_root):
    folder = public_root / "extensions/scaling"
    manifest = json.loads((folder / "artifact_manifest.json").read_text())["artifacts"]
    require(set(manifest) == {"400", "600"}, "scaling atlas inventory changed")
    for size in ("400", "600"):
        record = manifest[size]
        chunks = []
        for part in record["parts"]:
            raw = (folder / "artifacts" / part["filename"]).read_bytes()
            require(len(raw) == part["bytes"] and digest(raw) == part["sha256"],
                    "scaling archive part checksum mismatch")
            chunks.append(raw)
        raw = b"".join(chunks)
        require(len(raw) == record["bytes"] and digest(raw) == record["sha256"],
                "scaling archive checksum mismatch")
        yield from atlas_members(raw, record["files_sha256"], record["filename"])


def negative_controls():
    patch = {
        "domain_pixels": [0, 4], "halfwidth_pixels": "2", "band_halfwidth_pixels": "1",
        "slope": "1", "clearance_pixels": "1/4", "area_scale_um2": "1",
        "points": [{"x": "-1/2", "t": "2", "source": 0}], "nuclear_constraints": [],
    }
    check_patch(patch)
    bad_rna = deepcopy(patch)
    bad_rna["points"][0]["x"] = "0"
    bad_nucleus = deepcopy(patch)
    bad_nucleus["points"] = []
    bad_nucleus["nuclear_constraints"] = [[0, "1/4", "1", "3"]]
    for bad in (bad_rna, bad_nucleus):
        try:
            check_patch(bad)
        except ValueError:
            continue
        raise ValueError("source-feasibility negative control was accepted")


def main(data_root, public_root, output):
    negative_controls()
    records = []
    for generator in (base_atlases(data_root), external_atlases(public_root), scaling_atlases(public_root)):
        for raw, name in generator:
            records.append(audit_atlas(raw, name))
    unique = {r["atlas_id"]: r for r in records}
    result = {
        "scope": "Zero-area source witnesses in all reported scalar tissue atlases",
        "method": "Direct rational evaluation of h=0 against RNA clearance, nuclear ownership, band, slope and fixed-trace constraints; no optimizer imports",
        "archive_records": len(records), "distinct_atlases": len(unique),
        "totals_over_distinct_atlases": {
            key: sum(r[key] for r in unique.values())
            for key in ("patches", "rna_free_patches", "rna_constraints", "nuclear_constraints")
        },
        "negative_controls_rejected": 2, "failures": 0,
        "conclusion": "Every recorded patch has a feasible zero-area source witness. Omitting count-invariant patches therefore changes none of the reported spectra.",
        "limitation": "This verifies the source-feasibility premise on saved scalar inputs. It neither recalibrates geometry nor validates the general coupled-interface theorem.",
        "atlases": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "atlases"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True, help="Restored base scientific archive")
    parser.add_argument("--public-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results/source_feasibility.json")
    args = parser.parse_args()
    main(args.data_root, args.public_root, args.output)
