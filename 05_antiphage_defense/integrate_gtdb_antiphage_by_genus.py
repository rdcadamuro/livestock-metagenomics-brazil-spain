#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SAMPLE_START = 1035
SAMPLE_END = 1057

RESULTS_DIR = Path(
    "/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/OneDrive_1_21-02-2026/"
    "Spacers_All/Results_Taxon_Spacers"
)
IMG_DIR = RESULTS_DIR / "IMG_CLEAN"

BIN_CONTEXT_TSV = RESULTS_DIR / "anti_phage_context_bin_table.tsv"
SYNC_MANIFEST_TSV = RESULTS_DIR / "gtdb_sync_manifest_1035_1057.tsv"

# USER: set your base path here

GTDB_SYNC_ROOT = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/GTDBK/Downloaded_ubu_GTDBK")
# USER: set your base path here
GTDB_FALLBACK_ROOT = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/GTDBK/ubu_home_GTDBK")

OUT_BIN_GTDB = RESULTS_DIR / "anti_phage_bin_with_gtdb.tsv"
OUT_GENUS_SYSTEM = RESULTS_DIR / "anti_phage_genus_system_common.tsv"
OUT_SYSTEM_COVER = RESULTS_DIR / "anti_phage_system_genus_coverage.tsv"
OUT_HEATMAP_TSV = RESULTS_DIR / "anti_phage_genus_system_top_heatmap.tsv"
OUT_REPORT = RESULTS_DIR / "anti_phage_gtdb_summary_report.txt"
FIG_HEATMAP_PNG = IMG_DIR / "Fig_Genus_vs_AntiPhageSystem_Heatmap.png"
FIG_HEATMAP_SVG = IMG_DIR / "Fig_Genus_vs_AntiPhageSystem_Heatmap.svg"

TOP_GENUS = 25
TOP_SYSTEMS = 20


def canon_sample(x: object) -> str:
    s = str(x).strip()
    if s.startswith("P19109_"):
        return s.split("_", 1)[1]
    return s


def canon_bin(x: object) -> str:
    s = str(x).strip()
    if s.startswith("bin."):
        return s
    if s.startswith("bin_"):
        return s.replace("bin_", "bin.")
    if s.startswith("vRhyme_bin_"):
        return f"bin.{s.replace('vRhyme_bin_', '')}"
    return s


def parse_rank(classification: str, rank_prefix: str) -> str:
    if not isinstance(classification, str) or not classification:
        return "Unclassified_bacteria"
    for token in classification.split(";"):
        t = token.strip()
        if t.startswith(rank_prefix):
            value = t.split("__", 1)[1] if "__" in t else ""
            value = value.strip()
            return value if value else "Unclassified_bacteria"
    return "Unclassified_bacteria"


def load_manifest_or_build() -> pd.DataFrame:
    if SYNC_MANIFEST_TSV.exists() and SYNC_MANIFEST_TSV.stat().st_size > 0:
        m = pd.read_csv(SYNC_MANIFEST_TSV, sep="\t", dtype=str).fillna("")
        if {"sample", "selected_file", "status", "source"}.issubset(m.columns):
            return m

    rows: list[dict[str, str]] = []
    for sample_n in range(SAMPLE_START, SAMPLE_END + 1):
        sample = str(sample_n)
        p_sync = GTDB_SYNC_ROOT / f"P19109_{sample}" / "classify" / "gtdbtk.bac120.summary.tsv"
        p_fb = GTDB_FALLBACK_ROOT / f"P19109_{sample}" / "classify" / "gtdbtk.bac120.summary.tsv"
        if p_sync.exists() and p_sync.stat().st_size > 0:
            status, source, selected = "found", "remote_synced", p_sync
        elif p_fb.exists() and p_fb.stat().st_size > 0:
            status, source, selected = "found_fallback", "local_fallback", p_fb
        else:
            status, source, selected = "missing_classify", "none", Path("")
        rows.append(
            {
                "sample": sample,
                "status": status,
                "source": source,
                "selected_file": str(selected),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    bin_df = pd.read_csv(BIN_CONTEXT_TSV, sep="\t", dtype=str).fillna("")
    bin_df["sample"] = bin_df["sample"].map(canon_sample)
    bin_df["bin"] = bin_df["bin"].map(canon_bin)

    metadata_cols = {
        "sample",
        "bin",
        "anti_phage_padloc_total_bin",
        "anti_phage_load_bin",
        "PlotGroup",
        "Host",
        "Stage",
        "Digester",
        "Context",
    }
    system_cols = [c for c in bin_df.columns if c not in metadata_cols]

    for c in system_cols + ["anti_phage_load_bin"]:
        if c in bin_df.columns:
            bin_df[c] = pd.to_numeric(bin_df[c], errors="coerce").fillna(0.0)

    manifest = load_manifest_or_build()
    if not manifest.empty:
        manifest["sample"] = manifest["sample"].map(canon_sample)
    manifest.to_csv(SYNC_MANIFEST_TSV, sep="\t", index=False)

    gtdb_rows = []
    for _, r in manifest.iterrows():
        sample = str(r["sample"])
        selected_raw = str(r.get("selected_file", "")).strip()
        selected = Path(selected_raw) if selected_raw else Path("")
        source = str(r.get("source", "none"))
        if not selected_raw or selected_raw == "." or (not selected.is_file()) or selected.stat().st_size == 0:
            continue
        g = pd.read_csv(selected, sep="\t", dtype=str, usecols=["user_genome", "classification"]).fillna("")
        g["sample"] = sample
        g["bin"] = g["user_genome"].map(canon_bin)
        g["gtdb_source"] = source
        gtdb_rows.append(g[["sample", "bin", "classification", "gtdb_source"]])

    gtdb_df = (
        pd.concat(gtdb_rows, ignore_index=True)
        if gtdb_rows
        else pd.DataFrame(columns=["sample", "bin", "classification", "gtdb_source"])
    )
    if not gtdb_df.empty:
        gtdb_df = gtdb_df.drop_duplicates(subset=["sample", "bin"], keep="first")

    merged = bin_df.merge(gtdb_df, on=["sample", "bin"], how="left")
    merged["classification"] = merged["classification"].fillna("")
    merged["gtdb_source"] = merged["gtdb_source"].fillna("none")

    merged["bact_domain"] = merged["classification"].map(lambda x: parse_rank(x, "d__"))
    merged["bact_phylum"] = merged["classification"].map(lambda x: parse_rank(x, "p__"))
    merged["bact_class"] = merged["classification"].map(lambda x: parse_rank(x, "c__"))
    merged["bact_order"] = merged["classification"].map(lambda x: parse_rank(x, "o__"))
    merged["bact_family"] = merged["classification"].map(lambda x: parse_rank(x, "f__"))
    merged["bact_genus"] = merged["classification"].map(lambda x: parse_rank(x, "g__"))
    merged["bact_species"] = merged["classification"].map(lambda x: parse_rank(x, "s__"))
    merged.to_csv(OUT_BIN_GTDB, sep="\t", index=False)

    mags = merged[["sample", "bin", "bact_genus"] + system_cols].copy()
    mags = mags.drop_duplicates(subset=["sample", "bin"], keep="first")
    for c in system_cols:
        mags[c] = (pd.to_numeric(mags[c], errors="coerce").fillna(0.0) > 0).astype(int)
    mags["bact_genus"] = mags["bact_genus"].replace("", "Unclassified_bacteria").fillna("Unclassified_bacteria")

    genus_tot = mags.groupby("bact_genus").size().rename("n_mags_total_genus").reset_index()

    gs_rows: list[dict[str, object]] = []
    for system in system_cols:
        present = mags[mags[system] > 0].copy()
        if present.empty:
            continue
        grp = (
            present.groupby("bact_genus")
            .agg(
                n_mags_with_system=("bin", "count"),
                n_samples_with_system=("sample", "nunique"),
            )
            .reset_index()
        )
        grp = grp.merge(genus_tot, on="bact_genus", how="left")
        grp["pct_mags_genus_with_system"] = (
            100.0 * grp["n_mags_with_system"] / grp["n_mags_total_genus"].replace(0, np.nan)
        ).fillna(0.0)
        grp["system"] = system
        gs_rows.append(grp)

    genus_system = (
        pd.concat(gs_rows, ignore_index=True)
        if gs_rows
        else pd.DataFrame(
            columns=[
                "system",
                "bact_genus",
                "n_mags_with_system",
                "n_samples_with_system",
                "n_mags_total_genus",
                "pct_mags_genus_with_system",
            ]
        )
    )
    if not genus_system.empty:
        genus_system = genus_system[
            [
                "system",
                "bact_genus",
                "n_mags_with_system",
                "n_mags_total_genus",
                "pct_mags_genus_with_system",
                "n_samples_with_system",
            ]
        ].sort_values(["n_mags_with_system", "pct_mags_genus_with_system"], ascending=[False, False])
    genus_system.to_csv(OUT_GENUS_SYSTEM, sep="\t", index=False)

    cov_rows = []
    n_mags_total = len(mags)
    for system in system_cols:
        present = mags[mags[system] > 0]
        cov_rows.append(
            {
                "system": system,
                "n_mags_with_system": int(len(present)),
                "pct_mags_with_system": (100.0 * len(present) / n_mags_total) if n_mags_total > 0 else 0.0,
                "n_genera_with_system": int(present["bact_genus"].nunique()),
                "n_samples_with_system": int(present["sample"].nunique()),
            }
        )
    system_cov = pd.DataFrame(cov_rows).sort_values("n_mags_with_system", ascending=False)
    system_cov.to_csv(OUT_SYSTEM_COVER, sep="\t", index=False)

    top_genera = genus_tot.sort_values("n_mags_total_genus", ascending=False).head(TOP_GENUS)["bact_genus"].tolist()
    top_systems = system_cov.head(TOP_SYSTEMS)["system"].tolist()
    mat = pd.DataFrame(0.0, index=top_genera, columns=top_systems)
    for genus in top_genera:
        subset = mags[mags["bact_genus"] == genus]
        denom = len(subset)
        if denom == 0:
            continue
        for system in top_systems:
            mat.loc[genus, system] = 100.0 * subset[system].sum() / denom

    mat.to_csv(OUT_HEATMAP_TSV, sep="\t", index=True, index_label="bact_genus")

    if mat.shape[0] > 0 and mat.shape[1] > 0:
        fig_w = max(10.0, 0.45 * mat.shape[1] + 3.8)
        fig_h = max(8.0, 0.28 * mat.shape[0] + 3.2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=65, ha="right", fontsize=8)
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_yticklabels(mat.index, fontsize=8)
        ax.set_xlabel("Anti-phage system")
        ax.set_ylabel("Bacterial genus")
        ax.set_title("Top genera x anti-phage systems (% MAGs in genus with system)")
        cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
        cbar.set_label("% of MAGs in genus")
        plt.tight_layout()
        fig.savefig(FIG_HEATMAP_PNG, dpi=450, bbox_inches="tight")
        fig.savefig(FIG_HEATMAP_SVG, format="svg", bbox_inches="tight")
        plt.close(fig)

    mapped = (merged["bact_genus"] != "Unclassified_bacteria").sum()
    total_bins = len(merged)
    mapped_pct = (100.0 * mapped / total_bins) if total_bins > 0 else 0.0
    missing_samples = manifest.loc[~manifest["status"].isin(["found", "found_fallback"]), "sample"].tolist()

    top_system_lines = []
    for _, r in system_cov.head(12).iterrows():
        top_system_lines.append(
            f"{r['system']}: n_mags={int(r['n_mags_with_system'])}, "
            f"n_genera={int(r['n_genera_with_system'])}, pct_mags={r['pct_mags_with_system']:.2f}%"
        )

    genus_div = []
    for genus, grp in genus_system.groupby("bact_genus"):
        genus_div.append((genus, int((grp["n_mags_with_system"] > 0).sum())))
    genus_div = sorted(genus_div, key=lambda x: x[1], reverse=True)[:12]

    lines = [
        "Anti-phage systems by bacterial genus (GTDB-Tk integration)",
        "",
        f"Samples expected: {SAMPLE_END - SAMPLE_START + 1}",
        f"Samples in manifest: {manifest['sample'].nunique()}",
        f"Samples missing classify: {len(missing_samples)} ({', '.join(missing_samples) if missing_samples else 'none'})",
        "",
        f"Total MAG bins: {total_bins}",
        f"Bins with bacterial genus assigned: {mapped} ({mapped_pct:.2f}%)",
        "",
        "Top systems shared across genera:",
        *[f"- {x}" for x in top_system_lines],
        "",
        "Top genera by anti-phage system diversity (#systems present):",
        *[f"- {g}: {n_systems}" for g, n_systems in genus_div],
        "",
        f"Output table: {OUT_GENUS_SYSTEM}",
        f"Output table: {OUT_SYSTEM_COVER}",
        f"Heatmap table: {OUT_HEATMAP_TSV}",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[OK] {OUT_BIN_GTDB}")
    print(f"[OK] {OUT_GENUS_SYSTEM}")
    print(f"[OK] {OUT_SYSTEM_COVER}")
    print(f"[OK] {OUT_HEATMAP_TSV}")
    print(f"[OK] {OUT_REPORT}")
    if FIG_HEATMAP_PNG.exists():
        print(f"[OK] {FIG_HEATMAP_PNG}")


if __name__ == "__main__":
    main()
