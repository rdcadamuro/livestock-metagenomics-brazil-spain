#!/usr/bin/env python3
from __future__ import annotations

import csv
import shutil
from collections import defaultdict
from pathlib import Path


# USER: set your base path here


ROOT = Path("/mnt/nvme2/Metagenome_Spain_Miuvigs/Refresh_Genomad")
EXCLUDED_TSV = ROOT / "Spain_vRhyme_final_taxonomy_excluded_audit.tsv"
PHAGE_TSV = ROOT / "Spain_vRhyme_final_taxonomy_unified.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def concat_fastas(inputs: list[Path], output: Path) -> None:
    with output.open("w", encoding="utf-8") as out_handle:
        for path in inputs:
            with path.open("r", encoding="utf-8") as in_handle:
                for line in in_handle:
                    out_handle.write(line)


def move_if_exists(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.unlink()
    shutil.move(str(source), str(destination))
    return True


def main() -> int:
    excluded_rows = read_tsv(EXCLUDED_TSV)
    phage_rows = read_tsv(PHAGE_TSV)

    excluded_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    phage_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in excluded_rows:
        excluded_by_sample[row["sample"]].append(row)
    for row in phage_rows:
        phage_by_sample[row["sample"]].append(row)

    for sample, rows in sorted(excluded_by_sample.items()):
        sample_dir = ROOT / sample
        not_phage_dir = sample_dir / "Not_Phages"
        not_phage_bins_dir = not_phage_dir / "vRhyme_best_bins_fasta"
        not_phage_linked_dir = not_phage_dir / "linked_bins"
        not_phage_dir.mkdir(parents=True, exist_ok=True)
        not_phage_bins_dir.mkdir(parents=True, exist_ok=True)
        not_phage_linked_dir.mkdir(parents=True, exist_ok=True)

        moved_records: list[dict[str, str]] = []
        for row in rows:
            bin_id = row["linked_bin_id"]
            fasta_dir = sample_dir / "vRhyme" / "vRhyme_best_bins_fasta"
            linked_dir = sample_dir / "postbin_checkv" / "linked_bins"

            moved_any = False
            for ext in ["fasta", "faa", "ffn"]:
                src = fasta_dir / f"{bin_id}.{ext}"
                dst = not_phage_bins_dir / f"{bin_id}.{ext}"
                moved_any = move_if_exists(src, dst) or moved_any

            linked_src = linked_dir / f"{bin_id}.linked.fasta"
            linked_dst = not_phage_linked_dir / f"{bin_id}.linked.fasta"
            moved_any = move_if_exists(linked_src, linked_dst) or moved_any

            manifest_row = dict(row)
            manifest_row["files_moved"] = "yes" if moved_any else "no"
            moved_records.append(manifest_row)

        excluded_fieldnames = list(moved_records[0].keys()) if moved_records else list(rows[0].keys()) + ["files_moved"]
        write_tsv(not_phage_dir / "not_phage_bins.tsv", moved_records, excluded_fieldnames)

        phage_sample_rows = phage_by_sample.get(sample, [])
        if phage_sample_rows:
            write_tsv(sample_dir / "phage_bins.tsv", phage_sample_rows, list(phage_sample_rows[0].keys()))

        linked_dir = sample_dir / "postbin_checkv" / "linked_bins"
        linked_multifasta = sample_dir / "postbin_checkv" / f"{sample}.vRhyme_best_bins.linked.fasta"
        linked_files = sorted(linked_dir.glob("vRhyme_bin_*.linked.fasta"))
        concat_fastas(linked_files, linked_multifasta)

    print(f"samples_with_not_phages\t{len(excluded_by_sample)}")
    print(f"excluded_bins\t{len(excluded_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
