"""Build context-specific element-promoter graphs and gene-candidate tables."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from contextlib import contextmanager
import csv
from dataclasses import dataclass
import gzip
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Iterable

from .artifacts import sha256_file
from .contact_metadata import verify_reported_checksum
from .contacts import canonical_chromosome, write_json_atomic


ATTRIBUTE_RE = re.compile(
    r'(?:^|;\s*)(gene_id|gene_name|transcript_id)\s+"([^"]+)"'
)


PROMOTER_FIELDS = (
    "promoter_id",
    "gene_id",
    "gene_name",
    "chrom",
    "tss",
    "strand",
    "start",
    "end",
    "transcript_ids",
    "transcript_n",
    "selection_rule",
)

EDGE_FIELDS = (
    "context",
    "source_node_id",
    "target_node_id",
    "master_dhs_id",
    "promoter_id",
    "gene_id",
    "gene_name",
    "distance_bp",
    "promoter_active",
    "promoter_activity_score",
    "promoter_atac_accessible",
    "promoter_atac_signal",
    "promoter_h3k27ac_posterior",
    "contact_strategy",
    "contact_assay",
    "contact_match",
    "contact_resolution_bp",
    "contact_status",
    "same_contact_bin",
    "observed_balanced_contact",
    "powerlaw_expected_contact",
    "contact_weight",
    "observed_over_expected",
    "candidate_link_score",
)

GENE_FIELDS = (
    "context",
    "master_dhs_id",
    "element_chrom",
    "element_start",
    "element_end",
    "element_summit",
    "regulatory_class",
    "element_h3k27ac_posterior",
    "element_combined_activity",
    "element_blacklist_overlap",
    "gene_id",
    "gene_name",
    "promoter_count",
    "active_promoter_count",
    "best_promoter_id",
    "best_active_promoter_id",
    "minimum_distance_bp",
    "best_promoter_active",
    "best_promoter_atac_signal",
    "best_promoter_h3k27ac_posterior",
    "contact_strategy",
    "contact_assay",
    "contact_match",
    "contact_resolution_bp",
    "best_contact_status",
    "best_observed_balanced_contact",
    "best_powerlaw_expected_contact",
    "best_contact_weight",
    "best_observed_over_expected",
    "candidate_link_score",
    "candidate_gene_rank",
    "active_candidate_gene_rank",
    "contact_gene_rank",
    "nearest_gene_rank",
)

NODE_FIELDS = (
    "context",
    "node_id",
    "node_type",
    "master_dhs_id",
    "promoter_id",
    "gene_id",
    "gene_name",
    "chrom",
    "start",
    "end",
    "anchor",
    "regulatory_class",
    "atac_accessible",
    "atac_signal",
    "h3k27ac_posterior",
    "combined_activity",
    "active",
    "blacklist_overlap",
)


@dataclass(frozen=True)
class Promoter:
    promoter_id: str
    gene_id: str
    gene_name: str
    chrom: str
    tss: int
    strand: str
    start: int
    end: int
    transcript_ids: str


@dataclass(frozen=True)
class Element:
    master_dhs_id: str
    chrom: str
    start: int
    end: int
    summit: int
    regulatory_class: str
    activity_state: str
    atac_signal: float
    h3k27ac_posterior: float | None
    combined_activity: float
    blacklist_overlap: int


@dataclass(frozen=True)
class PromoterActivity:
    dhs_ids: str
    accessible: int
    atac_signal: float
    h3k27ac_posterior: float | None
    combined_activity: float
    activity_score: float
    active: int


@dataclass(frozen=True)
class ContactEvidence:
    status: str
    same_bin: int
    observed: float | None
    expected: float
    weight: float
    enrichment: float | None


def _open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8", newline="")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8", newline="")
    )


def _attributes(value: str) -> dict[str, str]:
    return {key: item for key, item in ATTRIBUTE_RE.findall(value)}


def _chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or int(fields[1]) < 1:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            sizes[canonical_chromosome(fields[0])] = int(fields[1])
    return sizes


def _centered_window(center: int, width: int, chromosome_size: int) -> tuple[int, int]:
    start = max(0, center - width // 2)
    end = min(chromosome_size, center + (width - width // 2))
    return start, end


def _atomic_tsv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with _output_text(temporary, compressed=path.suffix == ".gz") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fields),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _output_text(path: Path, *, compressed: bool):
    """Open a deterministic UTF-8 output stream for one temporary path."""

    if not compressed:
        with path.open("w", encoding="utf-8", newline="") as handle:
            yield handle
        return
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as compressed_handle:
            with io.TextIOWrapper(
                compressed_handle, encoding="utf-8", newline=""
            ) as text_handle:
                yield text_handle


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def build_promoter_table(
    *,
    annotation: Path,
    chrom_sizes: Path,
    canonical_chromosomes: list[str],
    promoter_width: int,
    annotation_checksum: str | None = None,
    output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Create one promoter node per distinct gene/TSS from the current GTF."""

    if promoter_width < 1:
        raise ValueError("Promoter width must be positive")
    if annotation_checksum:
        verify_reported_checksum(annotation, annotation_checksum)
    canonical = {canonical_chromosome(value) for value in canonical_chromosomes}
    sizes = _chrom_sizes(chrom_sizes)
    missing = canonical - set(sizes)
    if missing:
        raise ValueError(f"Chromosome sizes lack canonical sequences: {sorted(missing)}")

    grouped: dict[tuple[str, str, int, str], dict[str, Any]] = {}
    with _open_text(annotation) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            chromosome = canonical_chromosome(fields[0])
            if chromosome not in canonical:
                continue
            attributes = _attributes(fields[8])
            gene_id = attributes.get("gene_id", "")
            transcript_id = attributes.get("transcript_id", "")
            strand = fields[6]
            if not gene_id or not transcript_id or strand not in {"+", "-"}:
                continue
            tss = int(fields[3]) - 1 if strand == "+" else int(fields[4]) - 1
            key = (gene_id, chromosome, tss, strand)
            record = grouped.setdefault(
                key,
                {
                    "gene_name": attributes.get("gene_name", gene_id),
                    "transcripts": set(),
                },
            )
            record["transcripts"].add(transcript_id)
    if not grouped:
        raise ValueError(f"No canonical transcript promoters found in {annotation}")

    chromosome_order = {
        canonical_chromosome(chromosome): index
        for index, chromosome in enumerate(canonical_chromosomes)
    }
    ordered = sorted(
        grouped,
        key=lambda key: (chromosome_order[key[1]], key[2], key[0], key[3]),
    )
    rows = []
    for index, (gene_id, chromosome, tss, strand) in enumerate(ordered, start=1):
        start, end = _centered_window(tss, promoter_width, sizes[chromosome])
        record = grouped[(gene_id, chromosome, tss, strand)]
        transcripts = sorted(record["transcripts"])
        rows.append(
            {
                "promoter_id": f"DM6PROM{index:08d}",
                "gene_id": gene_id,
                "gene_name": record["gene_name"],
                "chrom": chromosome,
                "tss": tss,
                "strand": strand,
                "start": start,
                "end": end,
                "transcript_ids": ";".join(transcripts),
                "transcript_n": len(transcripts),
                "selection_rule": "all_distinct_gene_tss",
            }
        )
    _atomic_tsv(output, PROMOTER_FIELDS, rows)
    metrics = {
        "schema_version": 1,
        "annotation": str(annotation),
        "annotation_sha256": sha256_file(annotation),
        "annotation_reported_checksum": annotation_checksum,
        "chrom_sizes": str(chrom_sizes),
        "chrom_sizes_sha256": sha256_file(chrom_sizes),
        "canonical_chromosomes": canonical_chromosomes,
        "promoter_width_bp": promoter_width,
        "promoter_count": len(rows),
        "gene_count": len({row["gene_id"] for row in rows}),
        "output": str(output),
        "output_sha256": sha256_file(output),
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def read_promoters(path: Path) -> list[Promoter]:
    promoters = []
    seen: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(PROMOTER_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} lacks promoter fields: {sorted(missing)}")
        for row in reader:
            promoter_id = row["promoter_id"]
            if promoter_id in seen:
                raise ValueError(f"Duplicate promoter ID {promoter_id}")
            seen.add(promoter_id)
            promoters.append(
                Promoter(
                    promoter_id=promoter_id,
                    gene_id=row["gene_id"],
                    gene_name=row["gene_name"],
                    chrom=row["chrom"],
                    tss=int(row["tss"]),
                    strand=row["strand"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    transcript_ids=row["transcript_ids"],
                )
            )
    if not promoters:
        raise ValueError(f"No promoters in {path}")
    return promoters


def read_context_elements(path: Path, context: str) -> list[Element]:
    required = {
        "master_dhs_id",
        "chrom",
        "start",
        "end",
        "summit",
        "context",
        "context_membership",
        "regulatory_class",
        "atac_normalized_cpm_per_kb",
        "mixture_high_posterior_probability",
        "activity_state",
        "combined_activity_max_500",
        "blacklist_overlap",
    }
    elements = []
    seen: set[str] = set()
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} lacks catalog fields: {sorted(missing)}")
        for row in reader:
            if row["context"] != context or row["context_membership"] != "1":
                continue
            master_id = row["master_dhs_id"]
            if master_id in seen:
                raise ValueError(f"{path} repeats {master_id}")
            seen.add(master_id)
            posterior_text = row["mixture_high_posterior_probability"]
            posterior = float(posterior_text) if posterior_text else None
            elements.append(
                Element(
                    master_dhs_id=master_id,
                    chrom=row["chrom"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    summit=int(row["summit"]),
                    regulatory_class=row["regulatory_class"],
                    activity_state=row["activity_state"],
                    atac_signal=float(row["atac_normalized_cpm_per_kb"]),
                    h3k27ac_posterior=posterior,
                    combined_activity=float(row["combined_activity_max_500"]),
                    blacklist_overlap=int(row["blacklist_overlap"]),
                )
            )
    if not elements:
        raise ValueError(f"No context-member elements for {context} in {path}")
    return elements


def promoter_activities(
    promoters: list[Promoter],
    elements: list[Element],
    posterior_threshold: float,
) -> dict[str, PromoterActivity]:
    """Summarize context-member DHSs overlapping each fixed promoter window."""

    by_chromosome: dict[str, list[Element]] = defaultdict(list)
    for element in elements:
        by_chromosome[element.chrom].append(element)
    starts: dict[str, list[int]] = {}
    for chromosome, values in by_chromosome.items():
        values.sort(key=lambda value: (value.start, value.end, value.master_dhs_id))
        starts[chromosome] = [value.start for value in values]

    result = {}
    for promoter in promoters:
        values = by_chromosome.get(promoter.chrom, [])
        stop = bisect_left(starts.get(promoter.chrom, []), promoter.end)
        overlaps = [value for value in values[:stop] if value.end > promoter.start]
        posteriors = [
            element.h3k27ac_posterior
            for element in overlaps
            if element.h3k27ac_posterior is not None
        ]
        posterior = max(posteriors) if posteriors else None
        combined = max((element.combined_activity for element in overlaps), default=0.0)
        activity_score = max(
            (
                element.combined_activity * element.h3k27ac_posterior
                for element in overlaps
                if element.h3k27ac_posterior is not None
            ),
            default=0.0,
        )
        result[promoter.promoter_id] = PromoterActivity(
            dhs_ids=";".join(element.master_dhs_id for element in overlaps),
            accessible=int(bool(overlaps)),
            atac_signal=max((element.atac_signal for element in overlaps), default=0.0),
            h3k27ac_posterior=posterior,
            combined_activity=combined,
            activity_score=activity_score,
            active=int(
                bool(overlaps)
                and posterior is not None
                and posterior >= posterior_threshold
            ),
        )
    return result


def powerlaw_expected(distance: int, gamma: float, scale: float, minimum: int) -> float:
    return scale * max(distance, minimum) ** gamma


class ObservedContact:
    """Return balanced contact and explicit distance-model components."""

    def __init__(
        self,
        path: Path,
        *,
        gamma: float,
        scale: float,
        pseudocount_fraction: float,
    ):
        try:
            import cooler  # type: ignore
        except ImportError as error:
            raise RuntimeError("cooler is required for observed contact links") from error
        self.cooler = cooler.Cooler(str(path))
        if self.cooler.binsize is None:
            raise ValueError(f"Variable-bin Cooler is unsupported: {path}")
        self.resolution = int(self.cooler.binsize)
        self.gamma = gamma
        self.scale = scale
        self.pseudocount_fraction = pseudocount_fraction
        self.aliases = {
            canonical_chromosome(chromosome): chromosome
            for chromosome in self.cooler.chromnames
        }
        self._cached_chromosome = ""
        self._cached_matrix = None

    def __call__(
        self, chromosome: str, first_position: int, second_position: int
    ) -> ContactEvidence:
        distance = abs(first_position - second_position)
        expected = powerlaw_expected(
            distance, self.gamma, self.scale, self.resolution
        )
        source_chromosome = self.aliases.get(canonical_chromosome(chromosome))
        same_bin = int(first_position // self.resolution == second_position // self.resolution)
        if source_chromosome is None:
            return ContactEvidence(
                "chromosome_absent_from_matrix",
                same_bin,
                None,
                expected,
                self.pseudocount_fraction * expected,
                None,
            )
        if self._cached_chromosome != source_chromosome:
            self._cached_matrix = self.cooler.matrix(
                balance=True, sparse=True
            ).fetch(source_chromosome).tocsr()
            self._cached_chromosome = source_chromosome
        matrix = self._cached_matrix
        first_bin = first_position // self.resolution
        second_bin = second_position // self.resolution
        if first_bin >= matrix.shape[0] or second_bin >= matrix.shape[1]:
            return ContactEvidence(
                "position_outside_matrix",
                same_bin,
                None,
                expected,
                self.pseudocount_fraction * expected,
                None,
            )
        if same_bin:
            neighbors: list[float] = []
            if first_bin > 0:
                neighbors.append(float(matrix[first_bin - 1, first_bin]))
            if first_bin + 1 < matrix.shape[0]:
                neighbors.append(float(matrix[first_bin, first_bin + 1]))
            finite = [
                value for value in neighbors if math.isfinite(value) and value >= 0
            ]
            if not finite:
                return ContactEvidence(
                    "unbalanced_or_masked_adjacent_bins",
                    same_bin,
                    None,
                    expected,
                    self.pseudocount_fraction * expected,
                    None,
                )
            observed = max(finite)
            status = "adjacent_bin_proxy"
        else:
            observed = float(matrix[first_bin, second_bin])
            if not math.isfinite(observed) or observed < 0:
                return ContactEvidence(
                    "unbalanced_or_masked_bin",
                    same_bin,
                    None,
                    expected,
                    self.pseudocount_fraction * expected,
                    None,
                )
            status = "observed_matrix_pixel"
        return ContactEvidence(
            status,
            same_bin,
            observed,
            expected,
            observed + self.pseudocount_fraction * expected,
            observed / expected if expected > 0 else None,
        )


def _rank(values: list[dict[str, Any]], key) -> dict[str, int]:
    return {
        row["gene_id"]: index
        for index, row in enumerate(sorted(values, key=key), start=1)
    }


def _element_gene_rows(
    element: Element,
    edge_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in edge_rows:
        by_gene[row["gene_id"]].append(row)
    result = []
    for gene_id, rows in by_gene.items():
        best = max(
            rows,
            key=lambda row: (
                float(row["candidate_link_score"]),
                float(row["contact_weight"]),
                -int(row["distance_bp"]),
                row["promoter_id"],
            ),
        )
        active_rows = [row for row in rows if row["promoter_active"] == 1]
        best_active = (
            max(
                active_rows,
                key=lambda row: (
                    float(row["candidate_link_score"]),
                    float(row["contact_weight"]),
                    -int(row["distance_bp"]),
                    row["promoter_id"],
                ),
            )
            if active_rows
            else None
        )
        result.append(
            {
                "context": best["context"],
                "master_dhs_id": element.master_dhs_id,
                "element_chrom": element.chrom,
                "element_start": element.start,
                "element_end": element.end,
                "element_summit": element.summit,
                "regulatory_class": element.regulatory_class,
                "element_h3k27ac_posterior": _format(element.h3k27ac_posterior),
                "element_combined_activity": _format(element.combined_activity),
                "element_blacklist_overlap": element.blacklist_overlap,
                "gene_id": gene_id,
                "gene_name": best["gene_name"],
                "promoter_count": len(rows),
                "active_promoter_count": len(active_rows),
                "best_promoter_id": best["promoter_id"],
                "best_active_promoter_id": (
                    best_active["promoter_id"] if best_active else ""
                ),
                "minimum_distance_bp": min(int(row["distance_bp"]) for row in rows),
                "best_promoter_active": best["promoter_active"],
                "best_promoter_atac_signal": best["promoter_atac_signal"],
                "best_promoter_h3k27ac_posterior": best[
                    "promoter_h3k27ac_posterior"
                ],
                "contact_strategy": best["contact_strategy"],
                "contact_assay": best["contact_assay"],
                "contact_match": best["contact_match"],
                "contact_resolution_bp": best["contact_resolution_bp"],
                "best_contact_status": best["contact_status"],
                "best_observed_balanced_contact": best[
                    "observed_balanced_contact"
                ],
                "best_powerlaw_expected_contact": best[
                    "powerlaw_expected_contact"
                ],
                "best_contact_weight": best["contact_weight"],
                "best_observed_over_expected": best["observed_over_expected"],
                "candidate_link_score": best["candidate_link_score"],
                "_contact_rank_value": max(float(row["contact_weight"]) for row in rows),
                "_active_score": (
                    float(best_active["candidate_link_score"])
                    if best_active
                    else None
                ),
            }
        )
    candidate_ranks = _rank(
        result,
        lambda row: (
            -float(row["candidate_link_score"]),
            int(row["minimum_distance_bp"]),
            row["gene_id"],
        ),
    )
    contact_ranks = _rank(
        result,
        lambda row: (
            -float(row["_contact_rank_value"]),
            int(row["minimum_distance_bp"]),
            row["gene_id"],
        ),
    )
    nearest_ranks = _rank(
        result,
        lambda row: (int(row["minimum_distance_bp"]), row["gene_id"]),
    )
    active = [row for row in result if row["_active_score"] is not None]
    active_ranks = _rank(
        active,
        lambda row: (
            -float(row["_active_score"]),
            int(row["minimum_distance_bp"]),
            row["gene_id"],
        ),
    )
    for row in result:
        gene_id = row["gene_id"]
        row["candidate_gene_rank"] = candidate_ranks[gene_id]
        row["active_candidate_gene_rank"] = active_ranks.get(gene_id, "")
        row["contact_gene_rank"] = contact_ranks[gene_id]
        row["nearest_gene_rank"] = nearest_ranks[gene_id]
        row.pop("_contact_rank_value")
        row.pop("_active_score")
    return sorted(result, key=lambda row: int(row["candidate_gene_rank"]))


def build_context_links(
    *,
    context: str,
    context_elements: Path,
    promoters_path: Path,
    powerlaw_path: Path,
    strategy: str,
    contact_assay: str,
    contact_match: str,
    configured_resolution: int,
    maximum_distance: int,
    pseudocount_fraction: float,
    posterior_threshold: float,
    observed_contact_path: Path | None,
    nodes_output: Path,
    edges_output: Path,
    genes_output: Path,
    metrics_output: Path,
) -> dict[str, Any]:
    """Write one weighted graph layer and its element-to-gene projection."""

    if strategy not in {"observed", "powerlaw"}:
        raise ValueError(f"Unsupported contact strategy: {strategy}")
    if maximum_distance < 1 or configured_resolution < 1:
        raise ValueError("Contact distance and resolution must be positive")
    if not 0 <= pseudocount_fraction <= 1 or not 0 <= posterior_threshold <= 1:
        raise ValueError("Contact pseudocount and posterior threshold must be in [0, 1]")
    promoters = read_promoters(promoters_path)
    elements = read_context_elements(context_elements, context)
    activities = promoter_activities(promoters, elements, posterior_threshold)
    with powerlaw_path.open(encoding="utf-8") as handle:
        powerlaw = json.load(handle)
    if strategy == "observed":
        if observed_contact_path is None:
            raise ValueError(f"{context}: observed strategy requires a contact map")
        fit = powerlaw["contexts"].get(context)
        if fit is None:
            raise ValueError(f"{context}: observed contact decay fit is absent")
        contact_reader = ObservedContact(
            observed_contact_path,
            gamma=float(fit["gamma"]),
            scale=float(fit["scale"]),
            pseudocount_fraction=pseudocount_fraction,
        )
        resolution = contact_reader.resolution
        if resolution != configured_resolution:
            raise ValueError(
                f"{context}: contact resolution {resolution} differs from configured "
                f"{configured_resolution}"
            )
        gamma = float(fit["gamma"])
        scale = float(fit["scale"])
    else:
        if observed_contact_path is not None:
            raise ValueError(f"{context}: powerlaw strategy cannot receive a contact map")
        atlas_fit = powerlaw["atlas_powerlaw"]
        gamma = float(atlas_fit["gamma"])
        scale = float(atlas_fit["scale"])
        resolution = configured_resolution
        contact_reader = None

    promoters_by_chromosome: dict[str, list[Promoter]] = defaultdict(list)
    promoter_tss: dict[str, list[int]] = {}
    for promoter in promoters:
        promoters_by_chromosome[promoter.chrom].append(promoter)
    for chromosome, values in promoters_by_chromosome.items():
        values.sort(key=lambda promoter: (promoter.tss, promoter.promoter_id))
        promoter_tss[chromosome] = [promoter.tss for promoter in values]

    node_rows: list[dict[str, Any]] = []
    for element in elements:
        node_rows.append(
            {
                "context": context,
                "node_id": f"element:{element.master_dhs_id}",
                "node_type": "element",
                "master_dhs_id": element.master_dhs_id,
                "promoter_id": "",
                "gene_id": "",
                "gene_name": "",
                "chrom": element.chrom,
                "start": element.start,
                "end": element.end,
                "anchor": element.summit,
                "regulatory_class": element.regulatory_class,
                "atac_accessible": 1,
                "atac_signal": _format(element.atac_signal),
                "h3k27ac_posterior": _format(element.h3k27ac_posterior),
                "combined_activity": _format(element.combined_activity),
                "active": int(
                    element.h3k27ac_posterior is not None
                    and element.h3k27ac_posterior >= posterior_threshold
                ),
                "blacklist_overlap": element.blacklist_overlap,
            }
        )
    for promoter in promoters:
        activity = activities[promoter.promoter_id]
        node_rows.append(
            {
                "context": context,
                "node_id": f"promoter:{promoter.promoter_id}",
                "node_type": "promoter",
                "master_dhs_id": activity.dhs_ids,
                "promoter_id": promoter.promoter_id,
                "gene_id": promoter.gene_id,
                "gene_name": promoter.gene_name,
                "chrom": promoter.chrom,
                "start": promoter.start,
                "end": promoter.end,
                "anchor": promoter.tss,
                "regulatory_class": "promoter",
                "atac_accessible": activity.accessible,
                "atac_signal": _format(activity.atac_signal),
                "h3k27ac_posterior": _format(activity.h3k27ac_posterior),
                "combined_activity": _format(activity.combined_activity),
                "active": activity.active,
                "blacklist_overlap": "",
            }
        )
    _atomic_tsv(nodes_output, NODE_FIELDS, node_rows)

    edge_count = 0
    gene_count = 0
    observed_edge_count = 0
    active_edge_count = 0

    def edge_and_gene_rows():
        nonlocal edge_count, gene_count, observed_edge_count, active_edge_count
        for element in elements:
            chromosome_promoters = promoters_by_chromosome.get(element.chrom, [])
            tss_values = promoter_tss.get(element.chrom, [])
            left = bisect_left(tss_values, element.summit - maximum_distance)
            right = bisect_right(tss_values, element.summit + maximum_distance)
            element_edges = []
            for promoter in chromosome_promoters[left:right]:
                distance = abs(element.summit - promoter.tss)
                activity = activities[promoter.promoter_id]
                if contact_reader is None:
                    expected = powerlaw_expected(distance, gamma, scale, resolution)
                    evidence = ContactEvidence(
                        "distance_model",
                        int(element.summit // resolution == promoter.tss // resolution),
                        None,
                        expected,
                        expected,
                        None,
                    )
                else:
                    evidence = contact_reader(
                        element.chrom, element.summit, promoter.tss
                    )
                candidate_score = evidence.weight * activity.activity_score
                row = {
                    "context": context,
                    "source_node_id": f"element:{element.master_dhs_id}",
                    "target_node_id": f"promoter:{promoter.promoter_id}",
                    "master_dhs_id": element.master_dhs_id,
                    "promoter_id": promoter.promoter_id,
                    "gene_id": promoter.gene_id,
                    "gene_name": promoter.gene_name,
                    "distance_bp": distance,
                    "promoter_active": activity.active,
                    "promoter_activity_score": _format(activity.activity_score),
                    "promoter_atac_accessible": activity.accessible,
                    "promoter_atac_signal": _format(activity.atac_signal),
                    "promoter_h3k27ac_posterior": _format(
                        activity.h3k27ac_posterior
                    ),
                    "contact_strategy": strategy,
                    "contact_assay": contact_assay,
                    "contact_match": contact_match,
                    "contact_resolution_bp": resolution,
                    "contact_status": evidence.status,
                    "same_contact_bin": evidence.same_bin,
                    "observed_balanced_contact": _format(evidence.observed),
                    "powerlaw_expected_contact": _format(evidence.expected),
                    "contact_weight": _format(evidence.weight),
                    "observed_over_expected": _format(evidence.enrichment),
                    "candidate_link_score": _format(candidate_score),
                }
                element_edges.append(row)
                edge_count += 1
                observed_edge_count += int(evidence.observed is not None)
                active_edge_count += activity.active
            for row in element_edges:
                yield "edge", row
            gene_rows = _element_gene_rows(element, element_edges)
            gene_count += len(gene_rows)
            for row in gene_rows:
                yield "gene", row

    # Stream the expensive pair construction once into two atomic outputs.
    edges_output.parent.mkdir(parents=True, exist_ok=True)
    genes_output.parent.mkdir(parents=True, exist_ok=True)
    edge_descriptor, edge_temporary_name = tempfile.mkstemp(
        prefix=f".{edges_output.name}.", suffix=".tmp", dir=edges_output.parent
    )
    gene_descriptor, gene_temporary_name = tempfile.mkstemp(
        prefix=f".{genes_output.name}.", suffix=".tmp", dir=genes_output.parent
    )
    os.close(edge_descriptor)
    os.close(gene_descriptor)
    edge_temporary = Path(edge_temporary_name)
    gene_temporary = Path(gene_temporary_name)
    try:
        with _output_text(edge_temporary, compressed=True) as edge_handle, _output_text(
            gene_temporary, compressed=True
        ) as gene_handle:
            edge_writer = csv.DictWriter(
                edge_handle,
                fieldnames=EDGE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            gene_writer = csv.DictWriter(
                gene_handle,
                fieldnames=GENE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            edge_writer.writeheader()
            gene_writer.writeheader()
            for row_type, row in edge_and_gene_rows():
                if row_type == "edge":
                    edge_writer.writerow(row)
                else:
                    gene_writer.writerow(row)
        os.replace(edge_temporary, edges_output)
        os.replace(gene_temporary, genes_output)
    finally:
        edge_temporary.unlink(missing_ok=True)
        gene_temporary.unlink(missing_ok=True)

    metrics = {
        "schema_version": 1,
        "context": context,
        "contact_strategy": strategy,
        "contact_assay": contact_assay,
        "contact_match": contact_match,
        "contact_resolution_bp": resolution,
        "maximum_distance_bp": maximum_distance,
        "pseudocount_fraction": pseudocount_fraction,
        "promoter_posterior_threshold": posterior_threshold,
        "promoter_activity_method": "maximum over context-member master DHSs "
        "overlapping the fixed promoter window",
        "link_score": "contact_weight * maximum over overlapping DHSs of "
        "(combined_activity * h3k27ac_posterior)",
        "element_node_count": len(elements),
        "promoter_node_count": len(promoters),
        "active_promoter_count": sum(activity.active for activity in activities.values()),
        "element_promoter_edge_count": edge_count,
        "element_gene_candidate_count": gene_count,
        "observed_edge_count": observed_edge_count,
        "active_promoter_edge_count": active_edge_count,
        "powerlaw_gamma": gamma,
        "powerlaw_scale": scale,
        "inputs": {
            "context_elements": {
                "path": str(context_elements),
                "sha256": sha256_file(context_elements),
            },
            "promoters": {
                "path": str(promoters_path),
                "sha256": sha256_file(promoters_path),
            },
            "powerlaw": {
                "path": str(powerlaw_path),
                "sha256": sha256_file(powerlaw_path),
            },
            "observed_contact": (
                {
                    "path": str(observed_contact_path),
                    "sha256": sha256_file(observed_contact_path),
                }
                if observed_contact_path is not None
                else None
            ),
        },
        "outputs": {
            "nodes": {"path": str(nodes_output), "sha256": sha256_file(nodes_output)},
            "element_promoter_edges": {
                "path": str(edges_output),
                "sha256": sha256_file(edges_output),
            },
            "element_gene_candidates": {
                "path": str(genes_output),
                "sha256": sha256_file(genes_output),
            },
        },
    }
    write_json_atomic(metrics_output, metrics)
    return metrics


def aggregate_link_metrics(
    *,
    context_metric_paths: list[Path],
    source_manifest: Path,
    promoter_metrics: Path,
    contact_metrics: list[Path],
    powerlaw: Path,
    output_metrics: Path,
    output_provenance: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    contexts = []
    for path in context_metric_paths:
        with path.open(encoding="utf-8") as handle:
            contexts.append(json.load(handle))
    contexts.sort(key=lambda row: row["context"])
    metrics = {
        "schema_version": 1,
        "context_count": len(contexts),
        "observed_context_count": sum(
            row["contact_strategy"] == "observed" for row in contexts
        ),
        "powerlaw_context_count": sum(
            row["contact_strategy"] == "powerlaw" for row in contexts
        ),
        "element_promoter_edge_count": sum(
            int(row["element_promoter_edge_count"]) for row in contexts
        ),
        "element_gene_candidate_count": sum(
            int(row["element_gene_candidate_count"]) for row in contexts
        ),
        "contexts": {
            row["context"]: {
                key: row[key]
                for key in (
                    "contact_strategy",
                    "contact_assay",
                    "contact_match",
                    "contact_resolution_bp",
                    "element_node_count",
                    "promoter_node_count",
                    "active_promoter_count",
                    "element_promoter_edge_count",
                    "element_gene_candidate_count",
                )
            }
            for row in contexts
        },
    }
    write_json_atomic(output_metrics, metrics)
    inputs = [source_manifest, promoter_metrics, powerlaw, *contact_metrics, *context_metric_paths]
    provenance = {
        "schema_version": 1,
        "method": "context_contact_graph_v1",
        "inputs": {
            str(path): {"path": str(path), "sha256": sha256_file(path)}
            for path in inputs
        },
        "outputs": {
            "metrics": {
                "path": str(output_metrics),
                "sha256": sha256_file(output_metrics),
            }
        },
    }
    write_json_atomic(output_provenance, provenance)
    return metrics, provenance
