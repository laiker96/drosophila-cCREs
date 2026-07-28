#!/usr/bin/env python3
"""Build a portable IGV session for final ATAC-condition and ChIP outputs."""

from __future__ import annotations

import argparse
import os
import xml.etree.ElementTree as ET
from pathlib import Path
import tempfile


def relative_path(path: Path, session_path: Path) -> str:
    return Path(os.path.relpath(path, session_path.parent)).as_posix()


def add_track(
    resources: ET.Element,
    panel: ET.Element,
    *,
    path: Path,
    output: Path,
    name: str,
    color: str,
    signal: bool,
) -> None:
    track_id = relative_path(path, output)
    ET.SubElement(resources, "Resource", path=track_id)
    ET.SubElement(
        panel,
        "Track",
        id=track_id,
        name=name,
        color=color,
        renderer="BAR_CHART" if signal else "BASIC_FEATURE",
        height="52" if signal else "24",
        expand="false",
        visible="true",
        windowFunction="mean" if signal else "count",
    )


def _require(paths: list[Path]) -> None:
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing IGV inputs:\n" + "\n".join(map(str, missing)))


def _write_xml_if_changed(session: ET.Element, output: Path) -> None:
    ET.indent(session, space="  ")
    content = ET.tostring(session, encoding="utf-8", xml_declaration=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.read_bytes() == content:
        return
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{output.name}.",
            dir=output.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, output)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def build_session(
    atac_root: Path,
    output: Path,
    genome: str,
    locus: str,
    chip_root: Path | None = None,
    final_atac_only: bool = False,
    chip_one_per_context: bool = False,
    master_bed: Path | None = None,
) -> tuple[int, int, int]:
    conditions_root = atac_root / "conditions"
    conditions = sorted(path for path in conditions_root.iterdir() if path.is_dir())
    if not conditions:
        raise ValueError(f"No ATAC condition outputs found under {conditions_root}")

    output.parent.mkdir(parents=True, exist_ok=True)
    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="true",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name="ATAC conditions and ChIP")
    track_count = 0
    if master_bed is None:
        candidate = atac_root / "master" / "master_dhs.bed"
        master_bed = candidate if candidate.is_file() else None
    if master_bed is not None:
        _require([master_bed])
        add_track(
            resources,
            panel,
            path=master_bed,
            output=output,
            name="Master DHS registry",
            color="106,27,154",
            signal=False,
        )
        track_count += 1
    for condition_root in conditions:
        condition = condition_root.name
        label = condition.upper()
        peaks = condition_root / "peaks"
        tracks = condition_root / "tracks"
        consensus = peaks / f"{condition}.replicate-supported.bed"
        qpois = tracks / f"{condition}.qpois.bw"
        if qpois.is_file():
            inputs = [
                tracks / f"{condition}.MACS3-pileup.unscaled.bw",
                qpois,
                consensus,
            ]
            if not final_atac_only:
                inputs[2:2] = [
                    peaks / f"{condition}.candidates.narrowPeak",
                    peaks / f"{condition}.qpois-refined.bed",
                ]
            _require(inputs)
            specifications = [
                (inputs[0], f"{label} | MACS3 insertion pileup", "31,120,180", True),
                (inputs[1], f"{label} | qpois signal", "117,112,179", True),
            ]
            if not final_atac_only:
                specifications.extend(
                    [
                        (inputs[2], f"{label} | lenient candidates", "105,105,105", False),
                        (inputs[3], f"{label} | qpois-refined peaks", "230,85,13", False),
                    ]
                )
            specifications.append(
                (consensus, f"{label} | replicate-supported peaks", "0,145,130", False)
            )
        else:
            inputs = [
                tracks / f"{condition}.CPM.bw",
                peaks / f"{condition}.hmmratac.narrowPeak",
                consensus,
            ]
            _require(inputs)
            specifications = (
                (inputs[0], f"{label} | pooled CPM signal", "31,120,180", True),
                (inputs[1], f"{label} | pooled HMMRATAC peaks", "117,112,179", False),
                (inputs[2], f"{label} | replicate-supported peaks", "0,145,130", False),
            )
        for path, name, color, signal in specifications:
            add_track(
                resources,
                panel,
                path=path,
                output=output,
                name=name,
                color=color,
                signal=signal,
            )
            track_count += 1

    chip_samples = 0
    if chip_root is not None:
        bigwigs = sorted((chip_root / "tracks").glob("*.CPM.bw"))
        if chip_one_per_context:
            selected: dict[str, Path] = {}
            for bigwig in bigwigs:
                sample = bigwig.name.removesuffix(".CPM.bw")
                context = sample.split("_", 1)[0]
                selected.setdefault(context, bigwig)
            bigwigs = [selected[context] for context in sorted(selected)]
        for bigwig in bigwigs:
            sample = bigwig.name.removesuffix(".CPM.bw")
            peak_dir = chip_root / "peaks" / sample
            peaks = list(peak_dir.glob(f"{sample}_peaks.*Peak"))
            if len(peaks) != 1:
                raise FileNotFoundError(f"Expected one MACS3 peak file for {sample}")
            for path, name, color, signal in (
                (bigwig, f"{sample} | ChIP CPM signal", "221,126,32", True),
                (peaks[0], f"{sample} | MACS3 peaks", "202,61,52", False),
            ):
                add_track(
                    resources,
                    panel,
                    path=path,
                    output=output,
                    name=name,
                    color=color,
                    signal=signal,
                )
                track_count += 1
            chip_samples += 1

    _write_xml_if_changed(session, output)
    return len(conditions), chip_samples, track_count


def build_catalog_session(
    *,
    context: str,
    genome: str,
    atac_bigwig: Path,
    h3k27ac_bigwig: Path,
    context_dhs_bed: Path,
    master_dhs_bed: Path,
    active_elements_bed: Path,
    output: Path,
    locus: str = "All",
) -> int:
    """Build a portable five-track IGV session for one atlas context."""

    paths = [
        atac_bigwig,
        h3k27ac_bigwig,
        context_dhs_bed,
        master_dhs_bed,
        active_elements_bed,
    ]
    _require(paths)
    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="true",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name=f"{context} regulatory atlas")
    label = context.upper()
    specifications = [
        (
            atac_bigwig,
            f"{label} | mean ATAC Tn5 signal (background-TMM)",
            "31,120,180",
            True,
        ),
        (
            h3k27ac_bigwig,
            f"{label} | mean H3K27ac fragment coverage (background-TMM)",
            "221,126,32",
            True,
        ),
        (
            context_dhs_bed,
            f"{label} | context DHSs (master coordinates)",
            "0,145,130",
            False,
        ),
        (
            active_elements_bed,
            f"{label} | active regulatory elements",
            "202,61,52",
            False,
        ),
        (
            master_dhs_bed,
            "Master DHS registry",
            "106,27,154",
            False,
        ),
    ]
    for path, name, color, signal in specifications:
        add_track(
            resources,
            panel,
            path=path,
            output=output,
            name=name,
            color=color,
            signal=signal,
        )
    _write_xml_if_changed(session, output)
    return len(specifications)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("atac_root", type=Path, help="Run's results/.../atac directory")
    parser.add_argument("--chip-root", type=Path, help="Optional ChIP run result root")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome", default="dm6")
    parser.add_argument("--locus", default="All")
    parser.add_argument(
        "--master-bed",
        type=Path,
        help="Master DHS BED; defaults to ATAC_ROOT/master/master_dhs.bed when present",
    )
    parser.add_argument(
        "--final-atac-only",
        action="store_true",
        help="Include only pooled signal and replicate-supported ATAC peaks",
    )
    parser.add_argument(
        "--chip-one-per-context",
        action="store_true",
        help="Include only the first sorted ChIP replicate for each context",
    )
    args = parser.parse_args()
    condition_n, chip_n, track_n = build_session(
        args.atac_root.resolve(),
        args.output.resolve(),
        args.genome,
        args.locus,
        args.chip_root.resolve() if args.chip_root else None,
        args.final_atac_only,
        args.chip_one_per_context,
        args.master_bed.resolve() if args.master_bed else None,
    )
    print(
        f"Wrote {condition_n} ATAC conditions, {chip_n} ChIP samples, "
        f"and {track_n} tracks to {args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
