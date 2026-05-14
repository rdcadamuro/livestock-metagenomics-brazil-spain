#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# USER: set your base path here

BASE = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/Results/Spacers_AntiPhage_GTDB_Final_1035_1057")
TSV = BASE / "TSV"
OUT = BASE / "PAPER_READY"
FIG_MAIN = OUT / "FIG_MAIN"
FIG_SUPP = OUT / "FIG_SUPP"
TAB = OUT / "TABLES"

IN_BIN = TSV / "anti_phage_context_bin_table.tsv"
IN_BIN_GTDB = TSV / "anti_phage_bin_with_gtdb.tsv"
IN_STATS = TSV / "anti_phage_context_stats.tsv"
IN_PAIR = TSV / "anti_phage_biodigester_paired_medians.tsv"
IN_HET = TSV / "anti_phage_intrasample_heterogeneity.tsv"
IN_HEAT = TSV / "anti_phage_genus_system_top_heatmap.tsv"

COLOR = {
    "Broiler litter W/A": "#E69F00",
    "Broiler litter No/A": "#C44E00",
    "Biodigestor Swine Entrance No/A": "#009E73",
    "Biodigestor Swine In": "#56B4E9",
    "Biodigestor Swine Out": "#CC79A7",
    "Broiler litter": "#E69F00",
    "Biodigestor": "#009E73",
}
TOP_GENUS = 25
TOP_FAMILY = 25
TOP_SYSTEMS = 20


def bh_fdr(pvals: Iterable[float]) -> np.ndarray:
    p = np.array(list(pvals), dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = (ranked[i] * n) / rank
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def rank_biserial_from_u(u_stat: float, n1: int, n2: int) -> float:
    if n1 <= 0 or n2 <= 0:
        return np.nan
    return (2.0 * u_stat / (n1 * n2)) - 1.0


def standardize_labels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["sample"] = out["sample"].astype(str)
    out["PlotGroup_std"] = out["PlotGroup"].map(
        {
            "Chicken W/A": "Broiler litter W/A",
            "Chicken No/A": "Broiler litter No/A",
            "Biodigestor Swine Entrance No/A": "Biodigestor Swine Entrance No/A",
            "Biodigestor Swine In": "Biodigestor Swine In",
            "Biodigestor Swine Out": "Biodigestor Swine Out",
        }
    )
    return out


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "#333333",
            "grid.color": "#D9D9D9",
            "grid.alpha": 0.7,
            "grid.linestyle": "-",
            "svg.fonttype": "none",
            "figure.dpi": 140,
        }
    )


def is_unclassified_taxon(x: object) -> bool:
    s = str(x).strip()
    if s == "":
        return True
    low = s.lower()
    if low in {"na", "n/a", "none", "null", "nan", "unknown"}:
        return True
    return low.startswith("unclassified")


def infer_system_columns(df: pd.DataFrame) -> list[str]:
    meta = {
        "sample",
        "bin",
        "anti_phage_padloc_total_bin",
        "anti_phage_load_bin",
        "PlotGroup",
        "Host",
        "Stage",
        "Digester",
        "Context",
        "PlotGroup_std",
        "classification",
        "gtdb_source",
        "bact_domain",
        "bact_phylum",
        "bact_class",
        "bact_order",
        "bact_family",
        "bact_genus",
        "bact_species",
    }
    cols = [c for c in df.columns if c not in meta]
    return cols


def build_taxon_matrix(
    df: pd.DataFrame,
    taxon_col: str,
    system_cols: list[str],
    top_taxa: int,
    top_systems: int,
) -> pd.DataFrame:
    work = df.copy()
    for c in system_cols:
        work[c] = pd.to_numeric(work[c], errors="coerce").fillna(0.0)

    global_prev = (work[system_cols] > 0).sum(axis=0).sort_values(ascending=False)
    top_sys = global_prev.head(top_systems).index.tolist()

    valid = ~work[taxon_col].map(is_unclassified_taxon)
    w = work.loc[valid, [taxon_col] + top_sys].copy()
    taxa_counts = w.groupby(taxon_col).size().sort_values(ascending=False)
    taxa = taxa_counts.head(top_taxa).index.tolist()

    mat = pd.DataFrame(0.0, index=taxa, columns=top_sys)
    for t in taxa:
        sub = w[w[taxon_col] == t]
        if len(sub) == 0:
            continue
        mat.loc[t, :] = 100.0 * (sub[top_sys] > 0).sum(axis=0) / len(sub)
    mat.index.name = taxon_col
    return mat


def build_overlap_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    work = df.copy()
    valid = (~work["bact_genus"].map(is_unclassified_taxon)) & (~work["bact_family"].map(is_unclassified_taxon))
    w = work.loc[valid, ["bact_genus", "bact_family"]].drop_duplicates().copy()

    genus_to_family = w.groupby("bact_genus")["bact_family"].nunique().sort_values(ascending=False)
    family_to_genus = w.groupby("bact_family")["bact_genus"].nunique().sort_values(ascending=False)

    rows = []
    for tax, n in genus_to_family.items():
        rows.append({"relation": "genus_to_family", "taxon": tax, "linked_n": int(n)})
    for tax, n in family_to_genus.items():
        rows.append({"relation": "family_to_genus", "taxon": tax, "linked_n": int(n)})
    report = pd.DataFrame(rows).sort_values(["relation", "linked_n", "taxon"], ascending=[True, False, True])

    summary = {
        "n_genera_total": int(genus_to_family.shape[0]),
        "n_families_total": int(family_to_genus.shape[0]),
        "n_genera_with_multiple_families": int((genus_to_family > 1).sum()),
        "n_families_with_multiple_genera": int((family_to_genus > 1).sum()),
        "top_families_with_more_genera": family_to_genus.head(15),
    }
    return report, summary


def plot_main_context(df_bin: pd.DataFrame, df_stats: pd.DataFrame) -> dict:
    ctx = df_bin[df_bin["Context"].isin(["Broiler litter", "Biodigestor"])].copy()
    ctx["anti_phage_load_bin"] = pd.to_numeric(ctx["anti_phage_load_bin"], errors="coerce").fillna(0.0)

    order = ["Broiler litter", "Biodigestor"]
    a = ctx.loc[ctx["Context"] == "Broiler litter", "anti_phage_load_bin"].values
    b = ctx.loc[ctx["Context"] == "Biodigestor", "anti_phage_load_bin"].values
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = rank_biserial_from_u(float(mw.statistic), len(a), len(b))
    med_diff = float(np.median(b) - np.median(a))

    # q across primary load hypotheses
    tests = df_stats[df_stats["test"].isin(["Mann-Whitney", "Wilcoxon paired", "Kruskal-Wallis"])].copy()
    tests["p_value"] = pd.to_numeric(tests["p_value"], errors="coerce")
    tests = tests.dropna(subset=["p_value"])
    tests["q_bh"] = bh_fdr(tests["p_value"].values) if len(tests) > 0 else np.nan
    q_main = np.nan
    hit = tests[tests["analysis"] == "broiler_vs_biodigestor_load"]
    if not hit.empty:
        q_main = float(hit.iloc[0]["q_bh"])

    n_mags = ctx.groupby("Context")["bin"].count().to_dict()
    n_samples = ctx.groupby("Context")["sample"].nunique().to_dict()

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    palette = {k: COLOR[k] for k in order}
    sns.boxplot(
        data=ctx,
        x="Context",
        y="anti_phage_load_bin",
        order=order,
        palette=palette,
        ax=ax,
        fliersize=2,
        linewidth=1.0,
    )
    sns.stripplot(
        data=ctx,
        x="Context",
        y="anti_phage_load_bin",
        order=order,
        ax=ax,
        color="#1a1a1a",
        size=2,
        alpha=0.35,
        jitter=0.22,
    )
    ax.set_title("Main: anti-phage load per MAG by context (strict systems)", fontsize=11, pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("Anti-phage load per MAG (strict PADLOC + CRISPR_active)")

    y_max = float(ctx["anti_phage_load_bin"].max())
    y_line = y_max * 1.08 if y_max > 0 else 1.0
    h = max(0.3, y_max * 0.02)
    ax.plot([0, 0, 1, 1], [y_line, y_line + h, y_line + h, y_line], lw=1.1, c="#333333")
    txt = (
        f"p={mw.pvalue:.3g}  q={q_main:.3g}\n"
        f"rank-biserial={rbc:.3f}  medianΔ(Bio-Broiler)={med_diff:.2f}\n"
        f"n_MAG Broiler={n_mags.get('Broiler litter',0)}, Biodigestor={n_mags.get('Biodigestor',0)} | "
        f"n_sample Broiler={n_samples.get('Broiler litter',0)}, Biodigestor={n_samples.get('Biodigestor',0)}"
    )
    ax.text(0.5, y_line + h * 1.2, txt, ha="center", va="bottom", fontsize=8, color="#222222")
    ax.set_ylim(bottom=0, top=max(y_line + h * 4, y_max * 1.24 + 0.5))
    sns.despine(ax=ax)
    fig.tight_layout()

    fig.savefig(FIG_MAIN / "Fig_MAIN_Context_Boxplot_Strict.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_MAIN / "Fig_MAIN_Context_Boxplot_Strict.svg", format="svg", bbox_inches="tight")
    plt.close(fig)

    return {
        "test": "Mann-Whitney (Broiler litter vs Biodigestor)",
        "p_value": float(mw.pvalue),
        "q_value": float(q_main) if not math.isnan(q_main) else np.nan,
        "u_stat": float(mw.statistic),
        "rank_biserial": float(rbc),
        "median_diff_bio_minus_broiler": med_diff,
        "n_mag_broiler": int(n_mags.get("Broiler litter", 0)),
        "n_mag_biodigestor": int(n_mags.get("Biodigestor", 0)),
        "n_sample_broiler": int(n_samples.get("Broiler litter", 0)),
        "n_sample_biodigestor": int(n_samples.get("Biodigestor", 0)),
    }


def plot_supp_5group(df_bin: pd.DataFrame, df_stats: pd.DataFrame) -> None:
    order = [
        "Broiler litter W/A",
        "Broiler litter No/A",
        "Biodigestor Swine Entrance No/A",
        "Biodigestor Swine In",
        "Biodigestor Swine Out",
    ]
    dat = df_bin[df_bin["PlotGroup_std"].isin(order)].copy()
    dat["anti_phage_load_bin"] = pd.to_numeric(dat["anti_phage_load_bin"], errors="coerce").fillna(0.0)

    fig, ax = plt.subplots(figsize=(12.5, 6.2))
    palette = {k: COLOR[k] for k in order}
    sns.boxplot(
        data=dat,
        x="PlotGroup_std",
        y="anti_phage_load_bin",
        order=order,
        palette=palette,
        ax=ax,
        fliersize=2,
        linewidth=1.0,
    )
    sns.stripplot(
        data=dat,
        x="PlotGroup_std",
        y="anti_phage_load_bin",
        order=order,
        ax=ax,
        color="#1a1a1a",
        size=2,
        alpha=0.35,
        jitter=0.22,
    )
    ax.set_title("Supplementary: anti-phage load per MAG across context groups", fontsize=11)
    ax.set_xlabel("")
    ax.set_ylabel("Anti-phage load per MAG (strict PADLOC + CRISPR_active)")
    ax.tick_params(axis="x", rotation=18)

    kw = df_stats[(df_stats["analysis"] == "anti_phage_load_by_plotgroup") & (df_stats["test"] == "Kruskal-Wallis")]
    dunn_sig = df_stats[
        (df_stats["analysis"] == "anti_phage_load_by_plotgroup")
        & (df_stats["test"] == "Dunn_BH")
        & (pd.to_numeric(df_stats["p_value"], errors="coerce") < 0.05)
    ].copy()
    dunn_lines = []
    if not dunn_sig.empty:
        for _, r in dunn_sig.head(4).iterrows():
            dunn_lines.append(f"{r['detail'].split(' (raw_p=')[0]} | q={float(r['p_value']):.3g}")
    txt = []
    if not kw.empty:
        txt.append(f"Kruskal p={float(kw.iloc[0]['p_value']):.3g}")
    if dunn_lines:
        txt.append("Dunn-BH significant:")
        txt.extend([f"- {x}" for x in dunn_lines])
    if txt:
        ax.text(
            1.01,
            0.98,
            "\n".join(txt),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(facecolor="white", edgecolor="#CCCCCC", boxstyle="round,pad=0.35"),
        )
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "Fig_SUPP_Context_Boxplot_5Groups_Strict.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_Context_Boxplot_5Groups_Strict.svg", format="svg", bbox_inches="tight")
    plt.close(fig)


def plot_supp_paired(df_pair: pd.DataFrame, df_stats: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    if not df_pair.empty:
        for _, r in df_pair.iterrows():
            ax.plot([0, 1], [float(r["median_in"]), float(r["median_out"])], color="#666666", alpha=0.75, linewidth=1.1)
            ax.scatter([0], [float(r["median_in"])], color=COLOR["Biodigestor Swine In"], s=32, zorder=3)
            ax.scatter([1], [float(r["median_out"])], color=COLOR["Biodigestor Swine Out"], s=32, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Biodigestor Swine In", "Biodigestor Swine Out"])
        ax.set_ylabel("Sample median anti-phage load per MAG")
        ax.set_title(f"Supplementary: paired biodigestor contrast (n_pairs={len(df_pair)})")
        ww = df_stats[(df_stats["analysis"] == "biodigestor_paired_median_load") & (df_stats["test"] == "Wilcoxon paired")]
        mw = df_stats[(df_stats["analysis"] == "biodigestor_in_vs_out_load") & (df_stats["test"] == "Mann-Whitney")]
        lines = []
        if not mw.empty:
            lines.append(f"Mann-Whitney (bin-level): p={float(mw.iloc[0]['p_value']):.3g}")
        if not ww.empty:
            lines.append(f"Wilcoxon paired (sample medians): p={float(ww.iloc[0]['p_value']):.3g}")
        if lines:
            ax.text(
                0.02,
                0.98,
                "\n".join(lines),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=8,
                bbox=dict(facecolor="white", edgecolor="#CCCCCC", boxstyle="round,pad=0.30"),
            )
    else:
        ax.text(0.5, 0.5, "No paired biodigestor samples found", ha="center", va="center")
        ax.axis("off")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "Fig_SUPP_Biodigestor_Paired_Medians.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_Biodigestor_Paired_Medians.svg", format="svg", bbox_inches="tight")
    plt.close(fig)


def plot_supp_heterogeneity(df_het: pd.DataFrame) -> None:
    order = [
        "Broiler litter W/A",
        "Broiler litter No/A",
        "Biodigestor Swine Entrance No/A",
        "Biodigestor Swine In",
        "Biodigestor Swine Out",
    ]
    df = df_het.copy()
    df["mean_pairwise_bray"] = pd.to_numeric(df["mean_pairwise_bray"], errors="coerce")
    df["q_two_sided_bh"] = pd.to_numeric(df["q_two_sided_bh"], errors="coerce")
    df = df.dropna(subset=["mean_pairwise_bray"]).copy()
    df["PlotGroup_std"] = df["PlotGroup"].map(
        {
            "Chicken W/A": "Broiler litter W/A",
            "Chicken No/A": "Broiler litter No/A",
            "Biodigestor Swine Entrance No/A": "Biodigestor Swine Entrance No/A",
            "Biodigestor Swine In": "Biodigestor Swine In",
            "Biodigestor Swine Out": "Biodigestor Swine Out",
        }
    )
    df = df.sort_values("sample", key=lambda s: s.astype(int))
    colors = [COLOR.get(g, "#999999") for g in df["PlotGroup_std"]]

    fig, ax = plt.subplots(figsize=(12.5, 5.4))
    ax.bar(df["sample"].astype(str), df["mean_pairwise_bray"], color=colors, edgecolor="#333333", linewidth=0.25)
    for i, (_, r) in enumerate(df.iterrows()):
        if pd.notna(r["q_two_sided_bh"]) and float(r["q_two_sided_bh"]) < 0.05:
            ax.text(i, float(r["mean_pairwise_bray"]) + 0.005, "*", ha="center", va="bottom", fontsize=11)
    present = [g for g in order if (df["PlotGroup_std"] == g).any()]
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[g], ec="#333333", lw=0.25) for g in present]
    ax.legend(handles, present, frameon=False, fontsize=8, title="Context groups", title_fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.22))
    ax.set_title("Supplementary: intra-sample heterogeneity among MAGs (Bray-Curtis)")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Mean pairwise Bray-Curtis (within sample)")
    ax.tick_params(axis="x", rotation=70)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "Fig_SUPP_IntraSample_Heterogeneity_Mags.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_IntraSample_Heterogeneity_Mags.svg", format="svg", bbox_inches="tight")
    plt.close(fig)


def plot_heatmaps(df_heat: pd.DataFrame) -> None:
    mat = df_heat.copy()
    mat = mat.set_index("bact_genus")
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # no-dendrogram prevalence
    fig_w = max(10.5, 0.42 * mat.shape[1] + 4.0)
    fig_h = max(8.5, 0.28 * mat.shape[0] + 3.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=65, ha="right", fontsize=8)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xlabel("Anti-phage system")
    ax.set_ylabel("Bacterial genus")
    ax.set_title("Supplementary: prevalence heatmap (no dendrogram)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label("% MAGs in genus with system")
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Prevalence_NoDendrogram.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Prevalence_NoDendrogram.svg", format="svg", bbox_inches="tight")
    plt.close(fig)

    # clustered heatmap with dendrogram
    cg = sns.clustermap(
        mat,
        method="average",
        metric="euclidean",
        cmap="viridis",
        vmin=0,
        vmax=100,
        figsize=(fig_w, fig_h),
        row_cluster=True,
        col_cluster=True,
        linewidths=0.0,
        cbar_kws={"label": "% MAGs in genus with system"},
    )
    cg.ax_heatmap.set_xlabel("Anti-phage system")
    cg.ax_heatmap.set_ylabel("Bacterial genus")
    cg.ax_heatmap.set_title("Supplementary: clustered heatmap with dendrogram", pad=18)
    cg.fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Clustered_Dendrogram.png", dpi=450, bbox_inches="tight")
    cg.fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Clustered_Dendrogram.svg", format="svg", bbox_inches="tight")
    plt.close(cg.fig)


def plot_heatmaps_genus_family_panel(mat_genus: pd.DataFrame, mat_family: pd.DataFrame) -> None:
    # Individual genus heatmap
    fig_w_g = max(10.5, 0.42 * mat_genus.shape[1] + 4.0)
    fig_h_g = max(8.5, 0.28 * mat_genus.shape[0] + 3.2)
    fig, ax = plt.subplots(figsize=(fig_w_g, fig_h_g))
    im = ax.imshow(mat_genus.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(np.arange(mat_genus.shape[1]))
    ax.set_xticklabels(mat_genus.columns, rotation=65, ha="right", fontsize=8)
    ax.set_yticks(np.arange(mat_genus.shape[0]))
    ax.set_yticklabels(mat_genus.index, fontsize=8)
    ax.set_xlabel("Anti-phage system")
    ax.set_ylabel("Bacterial genus")
    ax.set_title("Supplementary: Genus-level anti-phage prevalence heatmap (Unclassified removed)")
    cb = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cb.set_label("% MAGs in genus with system")
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Genus_NoUnclassified.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Genus_NoUnclassified.svg", format="svg", bbox_inches="tight")
    plt.close(fig)

    # Individual family heatmap
    fig_w_f = max(10.5, 0.42 * mat_family.shape[1] + 4.0)
    fig_h_f = max(8.5, 0.28 * mat_family.shape[0] + 3.2)
    fig, ax = plt.subplots(figsize=(fig_w_f, fig_h_f))
    im = ax.imshow(mat_family.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(np.arange(mat_family.shape[1]))
    ax.set_xticklabels(mat_family.columns, rotation=65, ha="right", fontsize=8)
    ax.set_yticks(np.arange(mat_family.shape[0]))
    ax.set_yticklabels(mat_family.index, fontsize=8)
    ax.set_xlabel("Anti-phage system")
    ax.set_ylabel("Bacterial family")
    ax.set_title("Supplementary: Family-level anti-phage prevalence heatmap (Unclassified removed)")
    cb = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cb.set_label("% MAGs in family with system")
    fig.tight_layout()
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Family_NoUnclassified.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Family_NoUnclassified.svg", format="svg", bbox_inches="tight")
    plt.close(fig)

    # Combined 1x2 panel
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.5), constrained_layout=True)
    im1 = axes[0].imshow(mat_genus.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    axes[0].set_xticks(np.arange(mat_genus.shape[1]))
    axes[0].set_xticklabels(mat_genus.columns, rotation=65, ha="right", fontsize=7)
    axes[0].set_yticks(np.arange(mat_genus.shape[0]))
    axes[0].set_yticklabels(mat_genus.index, fontsize=7)
    axes[0].set_xlabel("Anti-phage system")
    axes[0].set_ylabel("Bacterial genus")
    axes[0].set_title("Genus")

    im2 = axes[1].imshow(mat_family.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    axes[1].set_xticks(np.arange(mat_family.shape[1]))
    axes[1].set_xticklabels(mat_family.columns, rotation=65, ha="right", fontsize=7)
    axes[1].set_yticks(np.arange(mat_family.shape[0]))
    axes[1].set_yticklabels(mat_family.index, fontsize=7)
    axes[1].set_xlabel("Anti-phage system")
    axes[1].set_ylabel("Bacterial family")
    axes[1].set_title("Family")

    cbar = fig.colorbar(im2, ax=axes.ravel().tolist(), fraction=0.02, pad=0.015)
    cbar.set_label("% MAGs with system")
    fig.suptitle("Supplementary panel: Genus vs Family anti-phage prevalence (Unclassified removed)", fontsize=12)
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Panel_Genus_Family_NoUnclassified.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_SUPP / "Fig_SUPP_Heatmap_Panel_Genus_Family_NoUnclassified.svg", format="svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG_MAIN.mkdir(parents=True, exist_ok=True)
    FIG_SUPP.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    setup_style()

    df_bin = pd.read_csv(IN_BIN, sep="\t", dtype=str).fillna("")
    df_bin = standardize_labels(df_bin)
    df_stats = pd.read_csv(IN_STATS, sep="\t", dtype=str).fillna("")
    df_pair = pd.read_csv(IN_PAIR, sep="\t", dtype=str).fillna("")
    df_het = pd.read_csv(IN_HET, sep="\t", dtype=str).fillna("")
    df_heat = pd.read_csv(IN_HEAT, sep="\t", dtype=str).fillna("")
    df_bin_gtdb = pd.read_csv(IN_BIN_GTDB, sep="\t", dtype=str).fillna("")

    ann = plot_main_context(df_bin, df_stats)
    plot_supp_5group(df_bin, df_stats)
    plot_supp_paired(df_pair, df_stats)
    plot_supp_heterogeneity(df_het)
    plot_heatmaps(df_heat)

    # New heatmaps: genus/family from final bin+GTDB table, removing only Unclassified
    df_bin_gtdb["sample"] = df_bin_gtdb["sample"].astype(str)
    df_bin_gtdb["bin"] = df_bin_gtdb["bin"].astype(str)
    system_cols = infer_system_columns(df_bin_gtdb)
    mat_genus = build_taxon_matrix(df_bin_gtdb, "bact_genus", system_cols, TOP_GENUS, TOP_SYSTEMS)
    mat_family = build_taxon_matrix(df_bin_gtdb, "bact_family", system_cols, TOP_FAMILY, TOP_SYSTEMS)

    mat_genus.to_csv(TAB / "heatmap_genus_matrix_no_unclassified.tsv", sep="\t", index=True, index_label="bact_genus")
    mat_family.to_csv(TAB / "heatmap_family_matrix_no_unclassified.tsv", sep="\t", index=True, index_label="bact_family")
    plot_heatmaps_genus_family_panel(mat_genus, mat_family)

    overlap_df, overlap_summary = build_overlap_tables(df_bin_gtdb)
    overlap_df.to_csv(TAB / "taxon_overlap_genus_family_report.tsv", sep="\t", index=False)
    overlap_lines = [
        "Taxon overlap summary (from anti_phage_bin_with_gtdb.tsv)",
        "",
        f"n_genera_total\t{overlap_summary['n_genera_total']}",
        f"n_families_total\t{overlap_summary['n_families_total']}",
        f"n_genera_with_multiple_families\t{overlap_summary['n_genera_with_multiple_families']}",
        f"n_families_with_multiple_genera\t{overlap_summary['n_families_with_multiple_genera']}",
        "",
        "Top families with more genera:",
    ]
    for fam, n in overlap_summary["top_families_with_more_genera"].items():
        overlap_lines.append(f"- {fam}: {int(n)}")
    (TAB / "taxon_overlap_summary.txt").write_text("\n".join(overlap_lines) + "\n", encoding="utf-8")

    ann_df = pd.DataFrame([ann])
    ann_df.to_csv(TAB / "main_panel_stats_annotation.tsv", sep="\t", index=False)

    lines = [
        "Paper-ready anti-phage figure suite",
        "",
        "Main figure:",
        "- Fig_MAIN_Context_Boxplot_Strict (Broiler litter vs Biodigestor; strict anti-phage)",
        "- Annotation includes p, q (BH on primary hypotheses), effect size, and n (MAGs/samples).",
        "",
        "Supplementary figures:",
        "- Fig_SUPP_Context_Boxplot_5Groups_Strict",
        "- Fig_SUPP_Biodigestor_Paired_Medians",
        "- Fig_SUPP_IntraSample_Heterogeneity_Mags",
        "- Fig_SUPP_Heatmap_Prevalence_NoDendrogram",
        "- Fig_SUPP_Heatmap_Clustered_Dendrogram",
        "- Fig_SUPP_Heatmap_Genus_NoUnclassified",
        "- Fig_SUPP_Heatmap_Family_NoUnclassified",
        "- Fig_SUPP_Heatmap_Panel_Genus_Family_NoUnclassified",
        "",
        f"Main annotation table: {TAB / 'main_panel_stats_annotation.tsv'}",
        f"Genus matrix: {TAB / 'heatmap_genus_matrix_no_unclassified.tsv'}",
        f"Family matrix: {TAB / 'heatmap_family_matrix_no_unclassified.tsv'}",
        f"Overlap report: {TAB / 'taxon_overlap_genus_family_report.tsv'}",
    ]
    (OUT / "paper_ready_readme.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[OK] {FIG_MAIN}")
    print(f"[OK] {FIG_SUPP}")
    print(f"[OK] {TAB / 'main_panel_stats_annotation.tsv'}")
    print(f"[OK] {TAB / 'heatmap_genus_matrix_no_unclassified.tsv'}")
    print(f"[OK] {TAB / 'heatmap_family_matrix_no_unclassified.tsv'}")
    print(f"[OK] {TAB / 'taxon_overlap_genus_family_report.tsv'}")
    print(f"[OK] {OUT / 'paper_ready_readme.txt'}")


if __name__ == "__main__":
    main()
