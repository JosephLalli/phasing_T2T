#!/usr/bin/env python3
"""Validate public whole-genome notebook inputs and an execution artifact."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path, PurePosixPath

NOTEBOOK = "notebooks/notebooks_whole_genome/make_plots.ipynb"
NOTEBOOK_HASH = "e8dfdaa8fc57779e7f2c3cc9733b4394307232c406754443d3aee5aca55fd770"
VERSIONS = {
    "numpy": "2.3.4",
    "polars": "1.29.0",
    "pandas": "2.2.3",
    "scipy": "1.16.2",
    "matplotlib": "3.9.2",
    "seaborn": "0.13.2",
    "pyarrow": "17.0.0",
}
IDENTITY = {
    "format_version": 1,
    "release": "NG-TR68126R",
    "canonical_notebook": NOTEBOOK,
    "file_count": 188,
    "notebook_sha256": NOTEBOOK_HASH,
}


def fail(message):
    raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save(directory, name, value):
    Path(directory).mkdir(parents=True, exist_ok=True)
    (Path(directory) / name).write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n"
    )


def environment():
    found = {}
    for name, expected in VERSIONS.items():
        found[name] = __import__(name).__version__
        if found[name] != expected:
            fail(f"{name} version {found[name]}, expected {expected}")
    if (
        subprocess.check_output(
            ["fc-match", "-f", "%{family}", "Arial"], text=True
        ).strip()
        != "Arial"
    ):
        fail("Arial font is unavailable")
    return found


def input_path(root, relative):
    if not isinstance(relative, str) or "\\" in relative:
        fail(f"invalid manifest path: {relative!r}")
    rel = PurePosixPath(relative)
    if (
        rel.is_absolute()
        or ".." in rel.parts
        or not (
            relative.startswith("intermediate_data_whole_genome/")
            or relative.startswith("imputation_statistics_whole_genome/")
        )
    ):
        fail(f"invalid manifest path: {relative!r}")
    path = root.joinpath(*rel.parts)
    if (
        any(part.is_symlink() for part in (root, *path.parents))
        or path.is_symlink()
        or not path.is_file()
    ):
        fail(f"missing or symlinked input: {relative}")
    return path


def preflight(repo, inputs):
    root = Path(inputs)
    manifest_file = root / "INPUT_MANIFEST.json"
    if manifest_file.is_symlink() or not manifest_file.is_file():
        fail("missing or symlinked INPUT_MANIFEST.json")
    manifest = json.loads(manifest_file.read_text())
    if {key: manifest.get(key) for key in IDENTITY} != IDENTITY:
        fail("manifest identity differs from canonical contract")
    files = manifest.get("files")
    if not isinstance(files, list) or len(files) != 188:
        fail("manifest file list must contain 188 files")
    total, seen = 0, set()
    categories = [0, 0]
    for entry in files:
        relative = entry.get("path") if isinstance(entry, dict) else None
        if relative in seen:
            fail(f"duplicate manifest path: {relative!r}")
        path = input_path(root, relative)
        seen.add(relative)
        categories[
            0 if relative.startswith("intermediate_data_whole_genome/") else 1
        ] += 1
        if path.stat().st_size != entry.get("bytes") or sha256(path) != entry.get(
            "sha256"
        ):
            fail(f"input checksum mismatch: {relative}")
        total += path.stat().st_size
    if manifest.get("total_bytes") != total or categories != [12, 176]:
        fail("manifest totals or categories differ from contract")
    source = Path(repo) / NOTEBOOK
    if sha256(source) != NOTEBOOK_HASH:
        fail("canonical notebook hash mismatch")
    return {
        "status": "PASS",
        "manifest_identity": IDENTITY,
        "manifest_sha256": sha256(manifest_file),
        "source_sha256": sha256(source),
        "environment_versions": environment(),
    }


def inventory(directory):
    root = Path(directory)
    return [
        {
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(p for p in root.rglob("*") if p.is_file())
    ]


def recipe(args):
    sources = {
        key: os.environ[key]
        for key in (
            "MOUNT_REPO_SOURCE",
            "MOUNT_INPUT_SOURCE",
            "MOUNT_INTERMEDIATE_SOURCE",
            "MOUNT_IMPUTATION_SOURCE",
            "MOUNT_FIGURES_SOURCE",
            "MOUNT_TABLES_SOURCE",
            "MOUNT_SCRATCH_SOURCE",
        )
    }
    save(
        args.scratch,
        "run_recipe.json",
        {
            "canonical_notebook": NOTEBOOK,
            "image_id": args.image_id,
            "cpus": os.environ["RUN_CPUS"],
            "memory": os.environ["RUN_MEMORY"],
            "network": "none",
            "read_only_root": True,
            "mount_sources": sources,
        },
    )


def verify(args):
    import nbformat

    cached = json.loads((Path(args.scratch) / "preflight_receipt.json").read_text())
    source_path = Path(args.repo) / NOTEBOOK
    if (
        cached.get("status") != "PASS"
        or cached.get("manifest_identity") != IDENTITY
        or sha256(source_path) != cached.get("source_sha256")
    ):
        fail("preflight receipt identity mismatch")
    result = {
        **cached,
        "status": "PASS" if args.exit_status == 0 else "FAIL",
        "exit_status": args.exit_status,
        "started_unix": args.started,
        "finished_unix": time.time(),
        "image_id": args.image_id,
    }
    result["duration_seconds"] = result["finished_unix"] - args.started
    try:
        source, executed = (
            nbformat.read(source_path, as_version=4),
            nbformat.read(args.executed, as_version=4),
        )
        nbformat.validate(executed)
        if [(c.cell_type, c.source) for c in source.cells] != [
            (c.cell_type, c.source) for c in executed.cells
        ]:
            fail("executed notebook cell sequence differs")
        code = [
            c for c in executed.cells if c.cell_type == "code" and str(c.source).strip()
        ]
        if (
            len(executed.cells) != 89
            or len(code) != 48
            or any(c.execution_count is None for c in code)
        ):
            fail("not all expected code cells executed")
        if any(o.output_type == "error" for c in code for o in c.outputs):
            fail("executed notebook contains an error output")
        result.update(
            {
                "observed_cells": len(executed.cells),
                "observed_nonempty_code_cells": len(code),
                "executed_notebook_sha256": sha256(args.executed),
            }
        )
    except Exception as exc:
        result.update(status="FAIL", verification_error=str(exc))
    result["environment_versions"] = environment()
    result["output_inventory"] = {
        "figures_whole_genome": inventory(args.figures),
        "tables_whole_genome": inventory(args.tables),
    }
    for path in (
        "/sys/fs/cgroup/memory.peak",
        "/sys/fs/cgroup/memory/memory.max_usage_in_bytes",
    ):
        if os.path.exists(path):
            result["peak_memory_bytes"] = Path(path).read_text().strip()
            break
    save(args.scratch, "run_receipt.json", result)
    return result["status"] == "PASS"


def main():
    parser = argparse.ArgumentParser()
    for name in ("repo", "inputs"):
        parser.add_argument(
            f"--{name}", default=f"/{'work/repo' if name == 'repo' else 'inputs'}"
        )
    parser.add_argument("--scratch", required=True)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--write-run-recipe", action="store_true")
    parser.add_argument("--executed")
    parser.add_argument("--figures")
    parser.add_argument("--tables")
    parser.add_argument("--exit-status", type=int, default=0)
    parser.add_argument("--started", type=float, default=0)
    parser.add_argument("--image-id", default="")
    args = parser.parse_args()
    try:
        if args.preflight:
            save(
                args.scratch,
                "preflight_receipt.json",
                preflight(args.repo, args.inputs),
            )
        if args.write_run_recipe:
            recipe(args)
        return 0 if not args.verify else int(not verify(args))
    except Exception as exc:
        save(
            args.scratch,
            "preflight_receipt.json",
            {"status": "FAIL", "failure": str(exc)},
        )
        save(args.scratch, "run_receipt.json", {"status": "FAIL", "failure": str(exc)})
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
