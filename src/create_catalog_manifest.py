#!/usr/bin/env python3
"""Create a checksummed manifest for one completed regulatory-element catalog."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
from pathlib import Path
import tempfile

from short_read_processing.artifacts import (
    CATALOG_FILE_FIELDS,
    read_catalog_manifest,
    sha256_file,
)
from short_read_processing.cli import cli_main


RELATIVE_FILES = {
    "catalog": "activity/catalog/master_elements_long.tsv.gz",
    "metrics": "activity/catalog/regulatory_element_metrics.json",
    "provenance": "activity/catalog/regulatory_element_provenance.json",
    "resolved_config": "provenance/configs/report.resolved_config.json",
}


def _write_if_changed(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_text(encoding="utf-8") == content:
        return
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{path.name}.",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    root = args.result_root.resolve()
    output = args.output.resolve()
    artifacts = {field: root / RELATIVE_FILES[field] for field in CATALOG_FILE_FIELDS}
    missing = [path for path in artifacts.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Catalog bundle is incomplete; missing: "
            + ", ".join(str(path) for path in missing)
        )
    resolved = json.loads(artifacts["resolved_config"].read_text(encoding="utf-8"))
    metrics = json.loads(artifacts["metrics"].read_text(encoding="utf-8"))
    contexts = resolved.get("activity", {}).get("contexts", [])
    row = {
        "genome": resolved.get("reference", {}).get("name", ""),
        "method": metrics.get("method", ""),
        "contexts": ",".join(str(context) for context in contexts),
        "source_project": resolved.get("project", ""),
        "source_run_id": resolved.get("run_id", ""),
    }
    for field, path in artifacts.items():
        row[field] = Path(os.path.relpath(path, output.parent)).as_posix()
        row[f"{field}_sha256"] = sha256_file(path)
    fieldnames = [
        "genome",
        "method",
        "contexts",
        "source_project",
        "source_run_id",
    ]
    fieldnames.extend(
        item for field in CATALOG_FILE_FIELDS for item in (field, f"{field}_sha256")
    )
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerow(row)
    _write_if_changed(output, buffer.getvalue())
    read_catalog_manifest(output)
    print(f"Wrote regulatory-element catalog manifest to {output}")
    return 0


if __name__ == "__main__":
    cli_main(main)
