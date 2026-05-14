#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


# USER: set your base path here


CHECKV_DB = Path("/mnt/nvme/Taxonomia/Proteins/Checkv/DB_CLEAN/checkv-db-v1.5_skip71")
# USER: set your base path here
CHECKV_BIN = Path("/mnt/nvme/conda_envs/mvip/bin/checkv")
# USER: set your base path here
CHECKV_ENV_BIN = Path("/mnt/nvme/conda_envs/mvip/bin")
# USER: set your base path here
GENOMAD_ROOT = Path("/mnt/nvme2/Metagenome_Spain_Miuvigs")
# USER: set your base path here
FASTQ_ROOT = Path("/mnt/nvme2/Tati/FastQ")
# USER: set your base path here
OUTPUT_ROOT = Path("/mnt/nvme2/Metagenome_Spain_Miuvigs/Genomad_Brazil")
RCLONE_REMOTE_ROOT = "ubu:/home/vMags/Brazil"
# USER: set your base path here
PHABOX_DB = Path("/mnt/nvme/Taxonomia/phabox_db_v2_1")
# USER: set your base path here
PHABOX_ENV_PYTHON = Path("/mnt/nvme/conda_envs/phabox-env/bin/python")
# USER: set your base path here
PHABOX_ENV_BIN = Path("/mnt/nvme/conda_envs/phabox-env/bin")
PHABOX_RUNNER = OUTPUT_ROOT / "run_phabox_task.py"
PHABOX_INPUT_DIR = OUTPUT_ROOT / "Brazil_phabox_inputs"
PHABOX_OUTPUT_DIR = OUTPUT_ROOT / "Brazil_phabox_outputs"
PHAGE_QC_LOG = OUTPUT_ROOT / "Brazil_phage_qc.log"

PHAGE_MARKERS = [
    "caudoviricetes",
    "uroviricota",
    "crassvirales",
    "inoviridae",
    "tubulavirales",
]
NON_PHAGE_MARKERS = [
    "nucleocytoviricota",
    "mimiviridae",
    "poxviridae",
    "orthopoxvirus",
    "megaviricetes",
    "pokkesviricetes",
    "imitervirales",
    "chitovirales",
]
AMBIGUOUS_TAXONOMY_LABELS = {
    "",
    "-",
    "unknown",
    "unclassified",
    "viruses",
    "viruses;;;;;;;;;;;;;;;",
    "no orfs found",
    "no hits to database",
    "hits not found in taxonomy files",
    "no lineage larger than threshold.",
    "filtered",
}


def run(cmd: list[str], log_path: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    if log_path is None:
        subprocess.run(cmd, check=True, env=env)
        return
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(cmd) + "\n")
        handle.flush()
        subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.STDOUT, env=env)
        handle.write("\n")


def checkv_run_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = f"{CHECKV_ENV_BIN}:{env.get('PATH', '')}"
    return env


def capture(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return proc.stdout


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line)
    return {name: "".join(parts) for name, parts in seqs.items()}


def write_fasta(path: Path, seqs: dict[str, str], keep_ids: set[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for seq_id, seq in seqs.items():
            if seq_id not in keep_ids:
                continue
            handle.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def write_fasta_from_map(path: Path, seqs: dict[str, str], ordered_ids: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for seq_id in ordered_ids:
            seq = seqs.get(seq_id)
            if seq is None:
                continue
            handle.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def concat_fastas(inputs: list[Path], output: Path) -> None:
    with output.open("w", encoding="utf-8") as out_handle:
        for path in inputs:
            with path.open("r", encoding="utf-8") as in_handle:
                for line in in_handle:
                    out_handle.write(line)
                if not line.endswith("\n"):
                    out_handle.write("\n")


def write_empty_tsv(path: Path, fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()


def iter_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_header: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(parts)))
                current_header = line[1:]
                parts = []
            else:
                parts.append(line)
    if current_header is not None:
        records.append((current_header, "".join(parts)))
    return records


def write_fasta_records(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def append_prefixed_fasta(input_path: Path, output_handle, sample: str) -> None:
    with input_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                header = raw[1:].rstrip("\n")
                parts = header.split(maxsplit=1)
                seq_id = parts[0]
                suffix = f" {parts[1]}" if len(parts) > 1 else ""
                output_handle.write(f">{sample}__{seq_id}{suffix}\n")
            else:
                output_handle.write(raw)


def normalize_bucket_label(label: str) -> str:
    normalized = label.strip().lower()
    normalized = re.sub(r"[\s\-]+", "_", normalized)
    return normalized


def classify_taxonomy_label(lineage: str) -> str:
    normalized = normalize_bucket_label(lineage)
    if normalized in AMBIGUOUS_TAXONOMY_LABELS:
        return "ambiguous"
    if any(marker in normalized for marker in NON_PHAGE_MARKERS):
        return "non_phage"
    if any(marker in normalized for marker in PHAGE_MARKERS):
        return "phage"
    if normalized.startswith("viruses"):
        return "ambiguous"
    return "ambiguous"


def choose_taxonomy_class(genomad_status: str, phabox_status: str) -> str:
    if {genomad_status, phabox_status} == {"phage", "non_phage"}:
        return "taxonomy_conflict"
    if genomad_status == "non_phage" or phabox_status == "non_phage":
        return "non_phage_excluded"
    if genomad_status == "phage":
        return "phage_clear"
    if genomad_status == "ambiguous" and phabox_status in {"phage", "ambiguous"}:
        return "phage_ambiguous"
    if phabox_status == "phage":
        return "phage_ambiguous"
    return "phage_ambiguous"


def sample_paths(sample: str) -> dict[str, Path]:
    sample_root = GENOMAD_ROOT / sample / "Annotation" / "_sample_level" / "genomad_selected_mediumplus"
    prefix = f"{sample}.vMAGs.selected_mediumplus.merged"
    return {
        "virus_fna": sample_root / f"{prefix}_summary" / f"{prefix}_virus.fna",
        "virus_summary": sample_root / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv",
        "taxonomy": sample_root / f"{prefix}_annotate" / f"{prefix}_taxonomy.tsv",
    }


def discover_samples() -> list[str]:
    samples: list[str] = []
    for path in sorted(GENOMAD_ROOT.glob("P19109_*/Annotation/_sample_level/genomad_selected_mediumplus")):
        samples.append(path.parts[-4])
    return samples


def discover_remote_samples() -> list[str]:
    stdout = capture(
        [
            "conda",
            "run",
            "-n",
            "bioinfo",
            "bash",
            "-lc",
            f'rclone --config /home/rafael/.config/rclone/rclone.conf lsf "{RCLONE_REMOTE_ROOT}"',
        ]
    )
    samples: list[str] = []
    for line in stdout.splitlines():
        entry = line.strip().rstrip("/")
        if entry.startswith("P19109_"):
            samples.append(entry)
    return sorted(samples)


def discover_fastq_samples() -> list[str]:
    samples = {
        path.name.split("_S", 1)[0]
        for path in FASTQ_ROOT.glob("P19109_*_R[12]_001.fastq.gz")
    }
    return sorted(samples)


def resolve_fastqs(sample: str) -> tuple[Path, Path]:
    r1 = sorted(FASTQ_ROOT.glob(f"{sample}_S*_R1_001.fastq.gz"))
    r2 = sorted(FASTQ_ROOT.glob(f"{sample}_S*_R2_001.fastq.gz"))
    if len(r1) != 1 or len(r2) != 1:
        raise FileNotFoundError(f"{sample}: FASTQs not found uniquely in {FASTQ_ROOT}")
    return r1[0], r2[0]


def prepare_vrhyme_fastqs(sample: str) -> tuple[Path, Path]:
    source_r1, source_r2 = resolve_fastqs(sample)
    reads_dir = OUTPUT_ROOT / sample / "reads_for_vrhyme"
    reads_dir.mkdir(parents=True, exist_ok=True)
    alias_r1 = reads_dir / f"{sample}_R1.fastq.gz"
    alias_r2 = reads_dir / f"{sample}_R2.fastq.gz"
    if alias_r1.exists() or alias_r1.is_symlink():
        alias_r1.unlink()
    if alias_r2.exists() or alias_r2.is_symlink():
        alias_r2.unlink()
    alias_r1.symlink_to(source_r1)
    alias_r2.symlink_to(source_r2)
    return alias_r1, alias_r2


def remote_sample_paths(sample: str) -> dict[str, str]:
    return {
        "remote_dir": f"{RCLONE_REMOTE_ROOT}/{sample}",
        "virus_fna": f"{RCLONE_REMOTE_ROOT}/{sample}/{sample}_summary/{sample}_virus.fna",
        "virus_summary": f"{RCLONE_REMOTE_ROOT}/{sample}/{sample}_summary/{sample}_virus_summary.tsv",
        "taxonomy": f"{RCLONE_REMOTE_ROOT}/{sample}/{sample}_annotate/{sample}_taxonomy.tsv",
    }


def ensure_remote_sample(sample: str) -> dict[str, Path]:
    outdir = OUTPUT_ROOT / sample
    download_dir = outdir / "downloaded_genomad"
    download_dir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "pipeline.log"

    remote = remote_sample_paths(sample)
    local_paths = {
        "virus_fna": download_dir / f"{sample}_summary" / f"{sample}_virus.fna",
        "virus_summary": download_dir / f"{sample}_summary" / f"{sample}_virus_summary.tsv",
        "taxonomy": download_dir / f"{sample}_annotate" / f"{sample}_taxonomy.tsv",
    }
    if all(path.exists() for path in local_paths.values()):
        return local_paths

    remote_to_local = [
        (remote["virus_fna"], local_paths["virus_fna"]),
        (remote["virus_summary"], local_paths["virus_summary"]),
        (remote["taxonomy"], local_paths["taxonomy"]),
    ]
    for remote_path, local_path in remote_to_local:
        local_path.parent.mkdir(parents=True, exist_ok=True)
        run(
            [
                "conda",
                "run",
                "-n",
                "bioinfo",
                "rclone",
                "--config",
                "/home/rafael/.config/rclone/rclone.conf",
                "copyto",
                remote_path,
                str(local_path),
            ],
            log_path=log_path,
        )

    missing = [str(path) for path in local_paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"{sample}: missing downloaded geNomad files: {', '.join(missing)}")
    return local_paths


def build_filtered_tables(sample: str, contamination_cutoff: float, keep_low_quality: bool) -> tuple[Path, int]:
    paths = ensure_remote_sample(sample)
    outdir = OUTPUT_ROOT / sample
    outdir.mkdir(parents=True, exist_ok=True)

    log_path = outdir / "pipeline.log"
    checkv_out = outdir / "checkv_contigs"
    filtered_dir = outdir / "filtered_contigs"
    filtered_dir.mkdir(exist_ok=True)

    for key, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"{sample}: missing {key} -> {path}")

    if not CHECKV_DB.exists():
        raise FileNotFoundError(f"CheckV DB not found: {CHECKV_DB}")

    if not (checkv_out / "quality_summary.tsv").exists():
        attempts = [2, 1]
        last_exc: Exception | None = None
        for attempt_threads in attempts:
            if checkv_out.exists():
                shutil.rmtree(checkv_out, ignore_errors=True)
            try:
                run(
                    [
                        str(CHECKV_BIN),
                        "end_to_end",
                        str(paths["virus_fna"]),
                        str(checkv_out),
                        "-d",
                        str(CHECKV_DB),
                        "-t",
                        str(attempt_threads),
                    ],
                    log_path=log_path,
                    env=checkv_run_env(),
                )
                last_exc = None
                break
            except subprocess.CalledProcessError as exc:
                last_exc = exc
                with log_path.open("a", encoding="utf-8") as handle:
                    handle.write(f"Retrying contig CheckV for {sample} after failure with -t {attempt_threads}\n")
        if last_exc is not None:
            raise last_exc

    virus_rows = read_tsv(paths["virus_summary"])
    virus_by_id = {row["seq_name"]: row for row in virus_rows}

    taxonomy_rows = read_tsv(paths["taxonomy"])
    taxonomy_by_id = {row["seq_name"]: row for row in taxonomy_rows}

    quality_rows = read_tsv(checkv_out / "quality_summary.tsv")
    keep_ids: set[str] = set()
    ordered_keep_ids: list[str] = []
    kept_rows: list[dict[str, str]] = []
    removed_rows: list[dict[str, str]] = []

    for row in quality_rows:
        seq_id = row["contig_id"]
        geom = virus_by_id.get(seq_id, {})
        tax = taxonomy_by_id.get(seq_id, {})
        checkv_quality = row.get("checkv_quality", "")
        contamination = float(row.get("contamination", "0") or 0)

        quality_ok = checkv_quality != "Not-determined"
        if not keep_low_quality:
            quality_ok = checkv_quality in {"Medium-quality", "High-quality", "Complete"}

        contamination_ok = contamination <= contamination_cutoff

        merged = {
            "sample": sample,
            "contig_id": seq_id,
            "checkv_quality": checkv_quality,
            "completeness": row.get("completeness", ""),
            "contamination": row.get("contamination", ""),
            "provirus": row.get("provirus", ""),
            "viral_genes": row.get("viral_genes", ""),
            "host_genes": row.get("host_genes", ""),
            "virus_score": geom.get("virus_score", ""),
            "n_hallmarks": geom.get("n_hallmarks", ""),
            "marker_enrichment": geom.get("marker_enrichment", ""),
            "taxonomy": geom.get("taxonomy", tax.get("lineage", "")),
            "taxonomy_agreement": tax.get("agreement", ""),
            "kept": "yes" if (quality_ok and contamination_ok) else "no",
            "remove_reason": "",
        }

        reasons: list[str] = []
        if not quality_ok:
            reasons.append(f"checkv_quality={checkv_quality}")
        if not contamination_ok:
            reasons.append(f"contamination>{contamination_cutoff}")
        merged["remove_reason"] = ";".join(reasons)

        if merged["kept"] == "yes":
            keep_ids.add(seq_id)
            ordered_keep_ids.append(seq_id)
            kept_rows.append(merged)
        else:
            removed_rows.append(merged)

    fieldnames = [
        "sample",
        "contig_id",
        "checkv_quality",
        "completeness",
        "contamination",
        "provirus",
        "viral_genes",
        "host_genes",
        "virus_score",
        "n_hallmarks",
        "marker_enrichment",
        "taxonomy",
        "taxonomy_agreement",
        "kept",
        "remove_reason",
    ]
    write_tsv(filtered_dir / f"{sample}.kept_contigs.tsv", kept_rows, fieldnames)
    write_tsv(filtered_dir / f"{sample}.removed_contigs.tsv", removed_rows, fieldnames)

    original_seqs = read_fasta(paths["virus_fna"])
    checkv_virus_fna = checkv_out / "viruses.fna"
    checkv_seqs = read_fasta(checkv_virus_fna) if checkv_virus_fna.exists() else {}
    seqs = dict(original_seqs)

    replaced_proviruses = 0
    missing_provirus_replacements: list[str] = []
    for row in kept_rows:
        seq_id = row["contig_id"]
        if row.get("provirus") == "Yes":
            if seq_id in checkv_seqs:
                seqs[seq_id] = checkv_seqs[seq_id]
                replaced_proviruses += 1
            else:
                missing_provirus_replacements.append(seq_id)

    missing = keep_ids.difference(seqs)
    if missing:
        raise RuntimeError(f"{sample}: {len(missing)} kept contigs missing from FASTA")

    filtered_fasta = filtered_dir / f"{sample}.filtered_virus_contigs.fa"
    write_fasta_from_map(filtered_fasta, seqs, ordered_keep_ids)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            f"Using CheckV-trimmed sequences for kept proviruses: {replaced_proviruses} replaced; "
            f"{len(missing_provirus_replacements)} missing replacements\n"
        )
        if missing_provirus_replacements:
            handle.write("Missing CheckV provirus replacements: " + ", ".join(missing_provirus_replacements[:20]) + "\n")
    return filtered_fasta, len(keep_ids)


def run_vrhyme(sample: str, filtered_fasta: Path, threads: int) -> None:
    outdir = OUTPUT_ROOT / sample
    log_path = outdir / "pipeline.log"
    vrhyme_out = outdir / "vRhyme"
    fastq1, fastq2 = prepare_vrhyme_fastqs(sample)

    if (vrhyme_out / "vRhyme_best_bins.0.membership.tsv").exists() or any(vrhyme_out.glob("vRhyme_best_bins.*.membership.tsv")):
        print(f"{sample}: vRhyme output already exists, skipping", flush=True)
        return

    attempts: list[int] = []
    for attempt_threads in sorted({min(threads, 2), 1}, reverse=True):
        if attempt_threads not in attempts:
            attempts.append(attempt_threads)

    last_exc: Exception | None = None
    for attempt_threads in attempts:
        if vrhyme_out.exists():
            shutil.rmtree(vrhyme_out)
        try:
            run(
                [
                    "conda",
                    "run",
                    "-n",
                    "vrhyme310",
                    "vRhyme",
                    "-i",
                    str(filtered_fasta),
                    "-r",
                    str(fastq1),
                    str(fastq2),
                    "-o",
                    str(vrhyme_out),
                    "-t",
                    str(attempt_threads),
                    "--method",
                    "longest",
                ],
                log_path=log_path,
            )
            last_exc = None
            break
        except subprocess.CalledProcessError as exc:
            last_exc = exc
            with log_path.open("a", encoding="utf-8") as handle:
                handle.write(f"Retrying vRhyme for {sample} after failure with -t {attempt_threads}\n")
    if last_exc is not None:
        raise last_exc


def locate_best_iteration(vrhyme_out: Path) -> tuple[Path | None, Path | None]:
    membership_paths = sorted(
        vrhyme_out.glob("vRhyme_best_bins.*.membership.tsv"),
        key=lambda path: int(path.name.split(".")[1]),
    )
    if not membership_paths:
        return None, None

    membership = membership_paths[-1]
    iteration = membership.name.split(".")[1]
    summary = vrhyme_out / f"vRhyme_best_bins.{iteration}.summary.tsv"
    if not summary.exists():
        raise FileNotFoundError(f"Missing summary file for best iteration: {summary}")
    return membership, summary


def run_postbin_checkv(sample: str, threads: int) -> tuple[Path | None, Path | None]:
    outdir = OUTPUT_ROOT / sample
    log_path = outdir / "pipeline.log"
    vrhyme_out = outdir / "vRhyme"
    membership_tsv, summary_tsv = locate_best_iteration(vrhyme_out)
    if membership_tsv is None or summary_tsv is None:
        return None, None

    best_bins_fasta = vrhyme_out / "vRhyme_best_bins_fasta"
    if not best_bins_fasta.exists():
        raise FileNotFoundError(f"{sample}: missing vRhyme_best_bins_fasta directory")

    postbin_dir = outdir / "postbin_checkv"
    postbin_dir.mkdir(exist_ok=True)
    linked_dir = postbin_dir / "linked_bins"
    linked_files = sorted(linked_dir.glob("*.linked.fasta")) if linked_dir.exists() else []
    if not linked_files:
        if linked_dir.exists():
            shutil.rmtree(linked_dir)
        run(
            [
                "conda",
                "run",
                "-n",
                "vrhyme310",
                "link_bin_sequences.py",
                "-i",
                str(best_bins_fasta),
                "-o",
                str(linked_dir),
                "-e",
                "fasta",
                "-n",
                "1500",
                "-c",
                "N",
            ],
            log_path=log_path,
        )
        linked_files = sorted(linked_dir.glob("*.linked.fasta"))

    if not linked_files:
        raise RuntimeError(f"{sample}: no linked vRhyme bin FASTAs were generated")

    linked_multifasta = postbin_dir / f"{sample}.vRhyme_best_bins.linked.fasta"
    if not linked_multifasta.exists():
        concat_fastas(linked_files, linked_multifasta)

    prodigal_dir = postbin_dir / "prodigal_linked_bins"
    prodigal_gff = prodigal_dir / f"{sample}.linked_bins.gff"
    if not prodigal_gff.exists():
        prodigal_dir.mkdir(exist_ok=True)
        run(
            [
                "conda",
                "run",
                "-n",
                "vrhyme310",
                "prodigal",
                "-i",
                str(linked_multifasta),
                "-a",
                str(prodigal_dir / f"{sample}.linked_bins.faa"),
                "-d",
                str(prodigal_dir / f"{sample}.linked_bins.ffn"),
                "-f",
                "gff",
                "-o",
                str(prodigal_gff),
                "-m",
                "-p",
                "meta",
            ],
            log_path=log_path,
        )

    checkv_bins = postbin_dir / "checkv_bins"
    if not (checkv_bins / "quality_summary.tsv").exists():
        attempts: list[int] = []
        for attempt_threads in sorted({min(threads, 2), 1}, reverse=True):
            if attempt_threads not in attempts:
                attempts.append(attempt_threads)

        last_exc: Exception | None = None
        for attempt_threads in attempts:
            if checkv_bins.exists():
                shutil.rmtree(checkv_bins, ignore_errors=True)
            try:
                run(
                    [
                        str(CHECKV_BIN),
                        "end_to_end",
                        str(linked_multifasta),
                        str(checkv_bins),
                        "-d",
                        str(CHECKV_DB),
                        "-t",
                        str(attempt_threads),
                    ],
                    log_path=log_path,
                    env=checkv_run_env(),
                )
                last_exc = None
                break
            except subprocess.CalledProcessError as exc:
                last_exc = exc
                with log_path.open("a", encoding="utf-8") as handle:
                    handle.write(f"Retrying post-bin CheckV for {sample} after failure with -t {attempt_threads}\n")
        if last_exc is not None:
            raise last_exc

    return membership_tsv, summary_tsv


def build_bin_taxonomy_table(sample: str, membership_tsv: Path | None, summary_tsv: Path | None) -> Path:
    outdir = OUTPUT_ROOT / sample
    result_path = outdir / "bin_taxonomy.tsv"
    fieldnames = [
        "sample",
        "bin",
        "linked_bin_id",
        "members",
        "proteins",
        "redundancy",
        "vrhyme_low_contamination",
        "checkv_quality",
        "bin_completeness",
        "bin_contamination",
        "bin_provirus",
        "consensus_taxonomy",
        "consensus_support",
        "taxonomy_labels_seen",
        "member_contigs",
    ]

    if membership_tsv is None or summary_tsv is None:
        write_empty_tsv(result_path, fieldnames)
        return result_path

    kept_rows = read_tsv(outdir / "filtered_contigs" / f"{sample}.kept_contigs.tsv")
    kept_by_id = {row["contig_id"]: row for row in kept_rows}
    membership_rows = read_tsv(membership_tsv)
    summary_rows = read_tsv(summary_tsv)
    summary_by_bin = {row["bin"]: row for row in summary_rows}

    postbin_quality = outdir / "postbin_checkv" / "checkv_bins" / "quality_summary.tsv"
    checkv_by_id: dict[str, dict[str, str]] = {}
    if postbin_quality.exists():
        checkv_by_id = {row["contig_id"]: row for row in read_tsv(postbin_quality)}

    members_by_bin: dict[str, list[str]] = defaultdict(list)
    for row in membership_rows:
        members_by_bin[row["bin"]].append(row["scaffold"])

    rows: list[dict[str, str]] = []
    for bin_id in sorted(members_by_bin, key=lambda value: int(value)):
        member_ids = members_by_bin[bin_id]
        linked_bin_id = f"vRhyme_bin_{bin_id}"
        checkv_row = checkv_by_id.get(linked_bin_id, {})
        taxonomies = [kept_by_id.get(contig_id, {}).get("taxonomy", "") for contig_id in member_ids]
        nonempty_taxonomies = [taxonomy for taxonomy in taxonomies if taxonomy]
        taxonomy_counter = Counter(nonempty_taxonomies)

        consensus_taxonomy = ""
        consensus_support = ""
        taxonomy_labels_seen = ""
        if taxonomy_counter:
            consensus_taxonomy, consensus_count = taxonomy_counter.most_common(1)[0]
            consensus_support = f"{consensus_count / len(nonempty_taxonomies):.4f}"
            taxonomy_labels_seen = "; ".join(
                f"{taxonomy} ({count})" for taxonomy, count in taxonomy_counter.most_common()
            )

        summary_row = summary_by_bin.get(bin_id, {})
        redundancy = summary_row.get("redundancy", "")
        low_contamination = ""
        if redundancy != "":
            low_contamination = "yes" if int(redundancy) <= 1 else "no"

        rows.append(
            {
                "sample": sample,
                "bin": bin_id,
                "linked_bin_id": linked_bin_id,
                "members": summary_row.get("members", str(len(member_ids))),
                "proteins": summary_row.get("proteins", ""),
                "redundancy": redundancy,
                "vrhyme_low_contamination": low_contamination,
                "checkv_quality": checkv_row.get("checkv_quality", ""),
                "bin_completeness": checkv_row.get("completeness", ""),
                "bin_contamination": checkv_row.get("contamination", ""),
                "bin_provirus": checkv_row.get("provirus", ""),
                "consensus_taxonomy": consensus_taxonomy,
                "consensus_support": consensus_support,
                "taxonomy_labels_seen": taxonomy_labels_seen,
                "member_contigs": ";".join(member_ids),
            }
        )

    write_tsv(result_path, rows, fieldnames)
    return result_path


def write_master_bin_taxonomy() -> Path:
    master_path = OUTPUT_ROOT / "Brazil_vRhyme_bin_taxonomy.tsv"
    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None

    for sample_table in sorted(OUTPUT_ROOT.glob("P19109_*/bin_taxonomy.tsv")):
        sample_rows = read_tsv(sample_table)
        if fieldnames is None:
            fieldnames = sample_rows[0].keys() if sample_rows else [
                "sample",
                "bin",
                "linked_bin_id",
                "members",
                "proteins",
                "redundancy",
                "vrhyme_low_contamination",
                "checkv_quality",
                "bin_completeness",
                "bin_contamination",
                "bin_provirus",
                "consensus_taxonomy",
                "consensus_support",
                "taxonomy_labels_seen",
                "member_contigs",
            ]
        rows.extend(sample_rows)

    if fieldnames is None:
        fieldnames = [
            "sample",
            "bin",
            "linked_bin_id",
            "members",
            "proteins",
            "redundancy",
            "vrhyme_low_contamination",
            "checkv_quality",
            "bin_completeness",
            "bin_contamination",
            "bin_provirus",
            "consensus_taxonomy",
            "consensus_support",
            "taxonomy_labels_seen",
            "member_contigs",
        ]
    write_tsv(master_path, rows, list(fieldnames))
    return master_path


def build_brazil_phabox_inputs() -> tuple[Path, Path, Path]:
    PHABOX_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    linked_fasta = PHABOX_INPUT_DIR / "Brazil_vRhyme_bins.linked.prefixed.fasta"
    protein_faa = PHABOX_INPUT_DIR / "Brazil_vRhyme_bins.linked.prefixed.faa"
    mapping_tsv = PHABOX_INPUT_DIR / "Brazil_vRhyme_bin_accession_map.tsv"

    mapping_rows: list[dict[str, str]] = []
    with linked_fasta.open("w", encoding="utf-8") as fasta_out, protein_faa.open("w", encoding="utf-8") as faa_out:
        for sample_table in sorted(OUTPUT_ROOT.glob("P19109_*/bin_taxonomy.tsv")):
            sample = sample_table.parent.name
            sample_rows = read_tsv(sample_table)
            linked_input = OUTPUT_ROOT / sample / "postbin_checkv" / f"{sample}.vRhyme_best_bins.linked.fasta"
            protein_input = OUTPUT_ROOT / sample / "postbin_checkv" / "prodigal_linked_bins" / f"{sample}.linked_bins.faa"
            if not linked_input.exists():
                raise FileNotFoundError(f"{sample}: missing linked bin fasta for PhaBOX2 input")
            if not protein_input.exists():
                raise FileNotFoundError(f"{sample}: missing linked bin protein fasta for PhaBOX2 input")

            append_prefixed_fasta(linked_input, fasta_out, sample)
            append_prefixed_fasta(protein_input, faa_out, sample)

            for row in sample_rows:
                mapping_rows.append(
                    {
                        "sample": sample,
                        "bin": row["bin"],
                        "linked_bin_id": row["linked_bin_id"],
                        "phabox_accession": f"{sample}__{row['linked_bin_id']}",
                    }
                )

    write_tsv(mapping_tsv, mapping_rows, ["sample", "bin", "linked_bin_id", "phabox_accession"])
    return linked_fasta, protein_faa, mapping_tsv


def run_phabox_dataset_task(task: str, contigs: Path, proteins: Path, threads: int) -> Path:
    if task not in {"phagcn", "phatyp"}:
        raise ValueError(f"Unsupported PhaBOX2 task: {task}")
    if not PHABOX_DB.exists():
        raise FileNotFoundError(f"PhaBOX2 DB not found: {PHABOX_DB}")
    if not PHABOX_ENV_PYTHON.exists():
        raise FileNotFoundError(f"PhaBOX2 python not found: {PHABOX_ENV_PYTHON}")
    if not PHABOX_RUNNER.exists():
        raise FileNotFoundError(f"PhaBOX2 helper runner not found: {PHABOX_RUNNER}")

    task_out = PHABOX_OUTPUT_DIR / task
    prediction = task_out / "final_prediction" / f"{task}_prediction.tsv"
    if prediction.exists():
        return prediction

    task_out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = f"{PHABOX_ENV_BIN}:{env.get('PATH', '')}"
    run(
        [
            str(PHABOX_ENV_PYTHON),
            str(PHABOX_RUNNER),
            "--task",
            task,
            "--dbdir",
            str(PHABOX_DB),
            "--outpth",
            str(task_out),
            "--contigs",
            str(contigs),
            "--proteins",
            str(proteins),
            "--threads",
            str(threads),
            "--len",
            "2000",
        ],
        log_path=PHAGE_QC_LOG,
        env=env,
    )
    if not prediction.exists():
        raise FileNotFoundError(f"PhaBOX2 {task} did not produce {prediction}")
    return prediction


def parse_numeric(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def remap_lifestyle(label: str) -> str:
    normalized = normalize_bucket_label(label)
    if normalized == "virulent":
        return "lytic"
    if normalized == "temperate":
        return "temperate"
    if normalized in {"filtered", "-", "", "unknown"}:
        return "-"
    return label


def finalize_brazil_phage_qc(phabox_threads: int) -> dict[str, Path]:
    master_base = write_master_bin_taxonomy()
    linked_fasta, protein_faa, mapping_tsv = build_brazil_phabox_inputs()
    phagcn_tsv = run_phabox_dataset_task("phagcn", linked_fasta, protein_faa, phabox_threads)
    phatyp_tsv = run_phabox_dataset_task("phatyp", linked_fasta, protein_faa, phabox_threads)

    mapping_rows = read_tsv(mapping_tsv)
    phabox_accession_by_key = {
        (row["sample"], row["linked_bin_id"]): row["phabox_accession"] for row in mapping_rows
    }
    phagcn_by_accession = {row["Accession"]: row for row in read_tsv(phagcn_tsv)}
    phatyp_by_accession = {row["Accession"]: row for row in read_tsv(phatyp_tsv)}

    checkv_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for quality_path in sorted(OUTPUT_ROOT.glob("P19109_*/postbin_checkv/checkv_bins/quality_summary.tsv")):
        sample = quality_path.parts[-4]
        for row in read_tsv(quality_path):
            checkv_by_key[(sample, row["contig_id"])] = row

    base_rows = read_tsv(master_base)
    merged_rows: list[dict[str, str]] = []
    for row in base_rows:
        sample = row["sample"]
        linked_bin_id = row["linked_bin_id"]
        phabox_accession = phabox_accession_by_key.get((sample, linked_bin_id), f"{sample}__{linked_bin_id}")
        phagcn_row = phagcn_by_accession.get(phabox_accession, {})
        phatyp_row = phatyp_by_accession.get(phabox_accession, {})
        checkv_row = checkv_by_key.get((sample, linked_bin_id), {})

        redundancy_value = row.get("redundancy", "")
        redundancy_int = int(redundancy_value) if redundancy_value not in {"", "NA"} else 999
        checkv_quality = checkv_row.get("checkv_quality", row.get("checkv_quality", ""))
        miuvig_quality = checkv_row.get("miuvig_quality", "")
        qc_primary_pass = (
            redundancy_int <= 1 and checkv_quality in {"Medium-quality", "High-quality", "Complete"}
        )

        genomad_taxonomy = row.get("consensus_taxonomy", "")
        phabox_taxonomy = phagcn_row.get("Lineage", "")
        genomad_status = classify_taxonomy_label(genomad_taxonomy)
        phabox_status = classify_taxonomy_label(phabox_taxonomy)
        taxonomy_class = choose_taxonomy_class(genomad_status, phabox_status)
        final_keep_phage = "yes" if taxonomy_class in {"phage_clear", "phage_ambiguous"} else "no"

        if final_keep_phage == "yes" and qc_primary_pass:
            final_bucket = "primary"
        elif final_keep_phage == "yes":
            final_bucket = "supplementary_only"
        else:
            final_bucket = "excluded"

        merged = dict(row)
        merged.update(
            {
                "phabox_accession": phabox_accession,
                "miuvig_quality": miuvig_quality,
                "bin_length": checkv_row.get("contig_length", phagcn_row.get("Length", "")),
                "bin_gene_count": checkv_row.get("gene_count", ""),
                "phabox2_taxonomy": phabox_taxonomy,
                "phabox2_taxonomy_score": phagcn_row.get("PhaGCNScore", ""),
                "phabox2_genus": phagcn_row.get("Genus", ""),
                "phabox2_genus_cluster": phagcn_row.get("GenusCluster", ""),
                "phabox2_lifestyle": remap_lifestyle(phatyp_row.get("TYPE", "")),
                "phabox2_lifestyle_score": phatyp_row.get("PhaTYPScore", ""),
                "genomad_taxonomy_status": genomad_status,
                "phabox_taxonomy_status": phabox_status,
                "taxonomy_class": taxonomy_class,
                "qc_primary_pass": "yes" if qc_primary_pass else "no",
                "final_keep_phage": final_keep_phage,
                "final_bucket": final_bucket,
            }
        )
        merged_rows.append(merged)

    fieldnames = list(merged_rows[0].keys()) if merged_rows else [
        "sample",
        "bin",
        "linked_bin_id",
        "members",
        "proteins",
        "redundancy",
        "vrhyme_low_contamination",
        "checkv_quality",
        "bin_completeness",
        "bin_contamination",
        "bin_provirus",
        "consensus_taxonomy",
        "consensus_support",
        "taxonomy_labels_seen",
        "member_contigs",
        "phabox_accession",
        "miuvig_quality",
        "bin_length",
        "bin_gene_count",
        "phabox2_taxonomy",
        "phabox2_taxonomy_score",
        "phabox2_genus",
        "phabox2_genus_cluster",
        "phabox2_lifestyle",
        "phabox2_lifestyle_score",
        "genomad_taxonomy_status",
        "phabox_taxonomy_status",
        "taxonomy_class",
        "qc_primary_pass",
        "final_keep_phage",
        "final_bucket",
    ]

    master_qc = OUTPUT_ROOT / "Brazil_vRhyme_phage_qc_master.tsv"
    primary_tsv = OUTPUT_ROOT / "Brazil_vRhyme_phage_primary_mediumplus.tsv"
    supplementary_all_tsv = OUTPUT_ROOT / "Brazil_vRhyme_phage_supplementary_all.tsv"
    supplementary_only_tsv = OUTPUT_ROOT / "Brazil_vRhyme_phage_supplementary_only.tsv"
    exclusion_tsv = OUTPUT_ROOT / "Brazil_vRhyme_phage_exclusion_audit.tsv"

    write_tsv(master_qc, merged_rows, fieldnames)
    write_tsv(primary_tsv, [row for row in merged_rows if row["final_bucket"] == "primary"], fieldnames)
    write_tsv(
        supplementary_all_tsv,
        [row for row in merged_rows if row["final_keep_phage"] == "yes"],
        fieldnames,
    )
    write_tsv(
        supplementary_only_tsv,
        [row for row in merged_rows if row["final_bucket"] == "supplementary_only"],
        fieldnames,
    )
    write_tsv(exclusion_tsv, [row for row in merged_rows if row["final_bucket"] == "excluded"], fieldnames)

    return {
        "master_qc": master_qc,
        "primary": primary_tsv,
        "supplementary_all": supplementary_all_tsv,
        "supplementary_only": supplementary_only_tsv,
        "exclusion": exclusion_tsv,
        "phagcn": phagcn_tsv,
        "phatyp": phatyp_tsv,
        "linked_fasta": linked_fasta,
        "protein_faa": protein_faa,
        "mapping": mapping_tsv,
    }


def build_smoketest_subset(
    sample: str,
    max_contigs: int,
    min_completeness: float,
) -> Path:
    outdir = OUTPUT_ROOT / sample / "vamb_smoketest"
    outdir.mkdir(parents=True, exist_ok=True)
    kept_tsv = OUTPUT_ROOT / sample / "filtered_contigs" / f"{sample}.kept_contigs.tsv"
    filtered_fasta = OUTPUT_ROOT / sample / "filtered_contigs" / f"{sample}.filtered_virus_contigs.fa"
    keep_ids_path = outdir / "keepids.txt"
    subset_fasta = outdir / f"{sample}.subset.fa"

    kept_rows = read_tsv(kept_tsv)
    selected_ids: list[str] = []
    for row in kept_rows:
        try:
            completeness = float(row["completeness"] or 0)
        except ValueError:
            completeness = 0
        if completeness < min_completeness:
            continue
        selected_ids.append(row["contig_id"])
        if len(selected_ids) >= max_contigs:
            break
    if len(selected_ids) < 2:
        raise RuntimeError(f"{sample}: not enough contigs for Vamb smoketest")

    keep_ids_path.write_text("\n".join(selected_ids) + "\n", encoding="utf-8")
    seqs = read_fasta(filtered_fasta)
    write_fasta(subset_fasta, seqs, set(selected_ids))
    return subset_fasta


def run_vamb_smoketest(sample: str, subset_fasta: Path, threads: int) -> None:
    smoketest_dir = OUTPUT_ROOT / sample / "vamb_smoketest"
    log_path = OUTPUT_ROOT / sample / "pipeline.log"
    vrhyme_subset = smoketest_dir / "vrhyme_subset"
    vamb_subset = smoketest_dir / "vamb_subset"
    fastq1, fastq2 = prepare_vrhyme_fastqs(sample)
    bam = vrhyme_subset / "vRhyme_bam_files" / f"{fastq1.name}.sorted.bam"

    if not (vrhyme_subset / "vRhyme_best_bins.0.summary.tsv").exists():
        run(
            [
                "conda",
                "run",
                "-n",
                "vrhyme310",
                "vRhyme",
                "-i",
                str(subset_fasta),
                "-r",
                str(fastq1),
                str(fastq2),
                "-o",
                str(vrhyme_subset),
                "-t",
                str(threads),
                "--method",
                "longest",
            ],
            log_path=log_path,
        )
    if not bam.exists():
        raise FileNotFoundError(f"{sample}: missing smoketest BAM for Vamb")
    if not bam.with_suffix(".bam.bai").exists():
        run(["conda", "run", "-n", "vrhyme310", "samtools", "index", str(bam)], log_path=log_path)
    if not (vamb_subset / "vae_clusters_unsplit.tsv").exists():
        run(
            [
                "conda",
                "run",
                "-n",
                "vamb310",
                "vamb",
                "bin",
                "default",
                "--outdir",
                str(vamb_subset),
                "--fasta",
                str(subset_fasta),
                "--bamdir",
                str(bam.parent),
                "-p",
                str(threads),
                "-m",
                "2000",
                "--minfasta",
                "2000",
                "-o",
                "",
            ],
            log_path=log_path,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Download geNomad virus outputs, filter with CheckV, run vRhyme, and summarize bins.")
    parser.add_argument("--samples", nargs="*", help="Brazil sample IDs. Default: auto-discover directly from ubu:/home/vMags/Brazil and keep only samples with local FASTQs.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--phabox-threads", type=int, default=2)
    parser.add_argument("--contamination-cutoff", type=float, default=10.0)
    parser.add_argument(
        "--drop-low-quality",
        action="store_true",
        help="If set, keep only Medium-quality/High-quality/Complete. Default keeps Low-quality and above, drops only Not-determined.",
    )
    parser.add_argument("--checkv-only", action="store_true", help="Run CheckV and filtering but skip vRhyme.")
    parser.add_argument("--restart", action="store_true", help="Delete existing Refresh_Genomad/<sample> directory before running.")
    parser.add_argument("--vamb-smoketest", action="store_true", help="Run vRhyme and Vamb on a small filtered subset. Intended for a single sample.")
    parser.add_argument("--finalize-phage-qc", action="store_true", help="After sample processing, run dataset-level PhaBOX2 taxonomy/lifestyle and write phage QC tables.")
    parser.add_argument("--finalize-phage-qc-only", action="store_true", help="Skip sample processing and run only dataset-level PhaBOX2 taxonomy/lifestyle and phage QC tables.")
    parser.add_argument("--vamb-max-contigs", type=int, default=20)
    parser.add_argument("--vamb-min-completeness", type=float, default=10.0)
    args = parser.parse_args()

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.finalize_phage_qc_only:
        outputs = finalize_brazil_phage_qc(phabox_threads=args.phabox_threads)
        for label, path in outputs.items():
            print(f"{label}\t{path}")
        return 0

    if args.samples:
        samples = args.samples
    else:
        remote_samples = set(discover_remote_samples())
        fastq_samples = set(discover_fastq_samples())
        samples = sorted(remote_samples & fastq_samples)
        missing_fastqs = sorted(remote_samples - fastq_samples)
        if missing_fastqs:
            print(f"Skipping {len(missing_fastqs)} remote samples without local FASTQs: {', '.join(missing_fastqs)}", file=sys.stderr)
    if not samples:
        print("No samples discovered.", file=sys.stderr)
        return 1
    if args.vamb_smoketest and len(samples) != 1:
        print("--vamb-smoketest currently requires exactly one sample.", file=sys.stderr)
        return 1

    failures: list[str] = []

    for sample in samples:
        try:
            print(f"=== {sample} ===", flush=True)
            sample_outdir = OUTPUT_ROOT / sample
            if args.restart and sample_outdir.exists():
                shutil.rmtree(sample_outdir)
            filtered_fasta, n_kept = build_filtered_tables(
                sample=sample,
                contamination_cutoff=args.contamination_cutoff,
                keep_low_quality=not args.drop_low_quality,
            )
            print(f"{sample}: kept {n_kept} contigs after CheckV filtering", flush=True)
            if n_kept < 2:
                print(f"{sample}: fewer than 2 contigs survived; skipping vRhyme", flush=True)
                build_bin_taxonomy_table(sample=sample, membership_tsv=None, summary_tsv=None)
                continue
            if not args.checkv_only:
                run_vrhyme(sample=sample, filtered_fasta=filtered_fasta, threads=args.threads)
                membership_tsv, summary_tsv = run_postbin_checkv(sample=sample, threads=args.threads)
                build_bin_taxonomy_table(sample=sample, membership_tsv=membership_tsv, summary_tsv=summary_tsv)
            if args.vamb_smoketest:
                subset_fasta = build_smoketest_subset(
                    sample=sample,
                    max_contigs=args.vamb_max_contigs,
                    min_completeness=args.vamb_min_completeness,
                )
                run_vamb_smoketest(sample=sample, subset_fasta=subset_fasta, threads=args.threads)
        except Exception as exc:
            failures.append(f"{sample}: {exc}")
            print(f"ERROR {sample}: {exc}", file=sys.stderr, flush=True)

    if failures:
        print("Failures:", file=sys.stderr)
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    write_master_bin_taxonomy()
    if args.finalize_phage_qc:
        outputs = finalize_brazil_phage_qc(phabox_threads=args.phabox_threads)
        for label, path in outputs.items():
            print(f"{label}\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
