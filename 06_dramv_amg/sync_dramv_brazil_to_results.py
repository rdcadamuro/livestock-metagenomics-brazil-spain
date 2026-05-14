#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import subprocess
from pathlib import Path
import pandas as pd

SAMPLES = [str(x) for x in range(1035, 1058)]
REMOTE_BASE = "ubu:home/Tati_Metagenome/DRAMv_VS2_phages_batch_2026-02-20"
# USER: set your base path here
LOCAL_BASE = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/DRamV")

FILES = [
    "amg_summary.tsv",
    "AMG_supporting_evidence.tsv",
    "AMG_strict_high_confidence.tsv",
    "AMG_potential_transfer.tsv",
    "DRAMv_AMG_final_report.txt",
    "quality_summary.tsv",
    "completeness.tsv",
    "contamination.tsv",
    "complete_genomes.tsv",
    "vMAG_stats.tsv",
    "product.html",
]
CORE = {
    "amg_summary.tsv",
    "AMG_supporting_evidence.tsv",
    "AMG_strict_high_confidence.tsv",
    "AMG_potential_transfer.tsv",
    "DRAMv_AMG_final_report.txt",
}


def run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)


def rclone_lsf(remote_dir: str) -> tuple[bool, list[str], str]:
    cmd = ["conda", "run", "-n", "bioinfo", "rclone", "lsf", remote_dir]
    r = run(cmd)
    if r.returncode != 0:
        return False, [], (r.stderr or r.stdout).strip()
    files = [x.strip() for x in r.stdout.splitlines() if x.strip()]
    return True, files, ""


def rclone_copyto(remote_file: str, local_file: Path) -> bool:
    local_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", "bioinfo", "rclone", "copyto", remote_file, str(local_file)]
    r = run(cmd)
    return r.returncode == 0


def main() -> None:
    LOCAL_BASE.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    for sample in SAMPLES:
        sample_tag = f"P19109_{sample}"
        remote_dir = f"{REMOTE_BASE}/{sample_tag}"
        local_dir = LOCAL_BASE / sample_tag
        local_dir.mkdir(parents=True, exist_ok=True)

        local_present = {f.name for f in local_dir.glob("*") if f.is_file()}
        core_local_count = len(CORE.intersection(local_present))
        if core_local_count == len(CORE):
            manifest_rows.append(
                {
                    "sample": sample,
                    "sample_tag": sample_tag,
                    "remote_dir": remote_dir,
                    "status": "already_local",
                    "source": "local",
                    "remote_ok": 1,
                    "remote_files_n": "NA",
                    "local_core_count": core_local_count,
                    "copied_n": 0,
                    "copied_files": "",
                    "missing_requested_files": "",
                    "missing_core_after_sync": "",
                    "note": "core files already present locally",
                }
            )
            continue

        ok, remote_files, err = rclone_lsf(remote_dir)
        if not ok or not remote_files:
            missing_core = sorted(list(CORE - local_present))
            manifest_rows.append(
                {
                    "sample": sample,
                    "sample_tag": sample_tag,
                    "remote_dir": remote_dir,
                    "status": "missing_remote",
                    "source": "none",
                    "remote_ok": 0,
                    "remote_files_n": 0,
                    "local_core_count": core_local_count,
                    "copied_n": 0,
                    "copied_files": "",
                    "missing_requested_files": ";".join(FILES),
                    "missing_core_after_sync": ";".join(missing_core),
                    "note": err or "remote sample folder not found",
                }
            )
            continue

        remote_set = set(remote_files)
        copied = []
        missing_requested = []
        for fname in FILES:
            if fname not in remote_set:
                missing_requested.append(fname)
                continue
            dst = local_dir / fname
            if dst.exists() and dst.stat().st_size > 0:
                continue
            src = f"{remote_dir}/{fname}"
            if rclone_copyto(src, dst):
                copied.append(fname)

        local_present = {f.name for f in local_dir.glob("*") if f.is_file()}
        core_local_count = len(CORE.intersection(local_present))
        missing_core = sorted(list(CORE - local_present))

        status = "found_synced" if core_local_count == len(CORE) else "missing_remote"
        source = "remote_synced" if copied else "local+remote_checked"
        note = ""
        if status != "found_synced":
            note = "incomplete core after sync"

        manifest_rows.append(
            {
                "sample": sample,
                "sample_tag": sample_tag,
                "remote_dir": remote_dir,
                "status": status,
                "source": source,
                "remote_ok": 1,
                "remote_files_n": len(remote_files),
                "local_core_count": core_local_count,
                "copied_n": len(copied),
                "copied_files": ";".join(copied),
                "missing_requested_files": ";".join(missing_requested),
                "missing_core_after_sync": ";".join(missing_core),
                "note": note,
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest = manifest.sort_values(by=["sample"])  # keep lexical sample ordering
    manifest_path = LOCAL_BASE / "DRAMv_sync_manifest_1035_1057.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    summary = (
        manifest.groupby(["status", "source"], dropna=False)
        .size()
        .reset_index(name="n_samples")
        .sort_values(["status", "source"])
    )
    summary_path = LOCAL_BASE / "DRAMv_sync_summary_1035_1057.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"[OK] manifest: {manifest_path}")
    print(f"[OK] summary: {summary_path}")


if __name__ == "__main__":
    main()
