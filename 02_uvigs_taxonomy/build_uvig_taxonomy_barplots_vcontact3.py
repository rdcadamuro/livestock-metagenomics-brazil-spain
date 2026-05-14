#!/usr/bin/env python3
"""
Barplot de taxonomia uViG usando saída do vConTACT3.
Input: brazil_spain_taxonomy_merged_groups.tsv
Output: FIG/Brazil_uViG_taxonomy_barplots_vcontact3.{png,svg}
        FIG/Spain_uViG_taxonomy_barplots_vcontact3.{png,svg}
"""

from __future__ import annotations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams["svg.fonttype"] = "none"

# USER: set your base path here

BASE   = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil")
OUTDIR = BASE / "FIG"
OUTTSV = BASE / "TSV"
MERGED = BASE / "vcontact3_runs/all_uvigs_brazil_spain/results/exports/brazil_spain_taxonomy_merged_groups.tsv"

BRAZIL_ORDER = [
    "Broiler with AMR",
    "Broiler without AMR",
    "Biodigestor in NO/A",
    "Biodigestor in",
    "Biodigestor out",
]

SPAIN_ORDER = [
    "Andalucia",
    "Aragon",
    "Asturias",
    "Castilla Leon",
    "Castilla-La Mancha",
    "Catalonia",
    "Extremadura",
    "Galicia",
    "Murcia",
    "Valencian Community",
]

LEVELS = ["class", "family", "genus"]
VC_COLS = {
    "class":  "vcontact3_class_prediction",
    "family": "vcontact3_family_prediction",
    "genus":  "vcontact3_genus_prediction",
}
OTHER_COLOR = "#BDBDBD"


def is_named(val) -> bool:
    """Retorna True apenas se o valor é um nome ICTV real (sem novel_, unplaced, pipe)."""
    if pd.isna(val):
        return False
    s = str(val)
    return not any(x in s.lower() for x in ["unplaced", "novel_", "|"])


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(MERGED, sep="\t", low_memory=False)

    # Criar coluna limpa por nível — None se não determinado
    for level, col in VC_COLS.items():
        df[level] = df[col].where(df[col].apply(is_named), other=None)

    brazil = df[df["dataset"] == "Brazil"].copy()
    spain  = df[df["dataset"] == "Spain"].copy()

    # Renomear sample_group -> group para compatibilidade
    brazil = brazil.rename(columns={"sample_group": "group"})
    spain  = spain.rename(columns={"sample_group": "group"})

    return brazil, spain


def summarize_level(
    df: pd.DataFrame,
    group_order: list[str],
    level: str,
    coverage_cutoff: float | None = None,
) -> tuple[pd.DataFrame, str | None]:
    tmp = df[["group", level]].copy()
    tmp = tmp.dropna(subset=[level])
    counts = tmp.groupby(["group", level], observed=True).size().reset_index(name="count")
    note = None

    if coverage_cutoff is not None and not counts.empty:
        totals = counts.groupby(level)["count"].sum().sort_values(ascending=False)
        frac = totals.cumsum() / totals.sum()
        keep = totals.index[frac <= coverage_cutoff].tolist()
        if not keep and not totals.empty:
            keep = [totals.index[0]]
        if keep and frac.loc[keep[-1]] < coverage_cutoff and len(keep) < len(totals):
            keep.append(totals.index[len(keep)])
        counts[level] = counts[level].where(counts[level].isin(keep), "Other")
        counts = counts.groupby(["group", level], observed=True)["count"].sum().reset_index()
        n_other = int((~totals.index.isin(keep)).sum())
        if n_other > 0:
            note = f"Other groups {n_other} additional {level} taxa."

    pivot = counts.pivot(index="group", columns=level, values="count").fillna(0)
    pivot = pivot.reindex(group_order, fill_value=0)
    ordered_cols = sorted(
        pivot.columns,
        key=lambda col: (col == "Other", -pivot[col].sum(), str(col)),
    )
    pivot = pivot[ordered_cols]
    return pivot, note


def palette(labels: list[str]) -> dict[str, str]:
    colors = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
    mapping: dict[str, str] = {}
    color_idx = 0
    for label in labels:
        if label == "Other":
            mapping[label] = OTHER_COLOR
        else:
            mapping[label] = colors[color_idx % len(colors)]
            color_idx += 1
    return mapping


def resolution_note(df: pd.DataFrame, rank: str) -> tuple[int, int, float]:
    total = len(df)
    resolved = int(df[rank].notna().sum())
    pct = (resolved / total * 100.0) if total else 0.0
    return total, resolved, pct


def plot_country(
    country: str,
    df: pd.DataFrame,
    group_order: list[str],
    png_name: str,
    svg_name: str,
    cutoff_map: dict[str, float | None],
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    OUTTSV.mkdir(parents=True, exist_ok=True)

    # Only include levels that have at least one named taxon
    active_levels = [lv for lv in LEVELS if df[lv].notna().any()]
    panel_labels = {lv: label for lv, label in zip(active_levels, ["A", "B", "C"])}

    n_panels = len(active_levels)
    fig, axes = plt.subplots(n_panels, 1, figsize=(15, max(8, 1.0 * len(group_order) * n_panels + 4)))
    if n_panels == 1:
        axes = [axes]
    fig.subplots_adjust(left=0.23, right=0.72, top=0.94, bottom=0.14, hspace=0.38)
    fig.suptitle(f"{country} uViG taxonomy by group (vConTACT3)", fontsize=22, fontweight="bold", y=0.98)

    notes: list[str] = []
    tsv_rows: list[dict] = []

    for ax, level in zip(axes, active_levels):
        pivot, collapse_note = summarize_level(df, group_order, level, cutoff_map.get(level))
        tax_labels = list(pivot.columns)
        colors = palette(tax_labels)
        totals = pivot.sum(axis=1).replace(0, np.nan)
        pct = pivot.div(totals, axis=0).multiply(100).fillna(0)

        y = np.arange(len(group_order))
        left = np.zeros(len(group_order))

        for taxon in tax_labels:
            values = pct[taxon].values
            ax.barh(y, values, left=left, color=colors[taxon],
                    edgecolor="white", linewidth=0.5, label=taxon)
            for idx, val in enumerate(values):
                if val >= 7:
                    ax.text(left[idx] + val / 2, y[idx], f"{val:.0f}%",
                            ha="center", va="center", fontsize=9,
                            color="white", fontweight="bold")
            left += values

            for group, count in pivot[taxon].items():
                tsv_rows.append({
                    "country": country, "group": group,
                    "level": level, "taxon": taxon, "count": int(count),
                })

        ax.set_xlim(0, 100)
        ax.set_yticks(y)
        ax.set_yticklabels(group_order, fontsize=11)
        ax.invert_yaxis()
        ax.set_xlabel("Relative abundance within ICTV-named bins (%)", fontsize=11)
        panel_letter = panel_labels[level]
        ax.set_title(f"{panel_letter}  {level.title()}", loc="left", fontsize=16, fontweight="bold")
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)
        ax.set_axisbelow(True)

        total, resolved, resolved_pct = resolution_note(df, level)
        if level == "class":
            notes.append(f"* Class: {resolved_pct:.1f}% of uViGs with ICTV-named class ({resolved}/{total}).")
        elif level == "family":
            notes.append(f"** Family: {resolved_pct:.1f}% of uViGs with ICTV-named family ({resolved}/{total}).")
        elif level == "genus":
            notes.append(f"*** Genus: {resolved_pct:.1f}% of uViGs with ICTV-named genus ({resolved}/{total}).")
        if collapse_note:
            notes.append(f"    {collapse_note}")

        ncol = 1 if len(tax_labels) <= 15 else 2 if len(tax_labels) <= 32 else 3
        ax.legend(
            title=f"{level.title()} taxa",
            bbox_to_anchor=(1.01, 1.0),
            loc="upper left",
            frameon=False,
            fontsize=9,
            title_fontsize=10,
            ncol=ncol,
        )

    fig.text(0.23, 0.04, "\n".join(notes), ha="left", va="bottom", fontsize=10)
    fig.savefig(OUTDIR / png_name, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / svg_name, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTDIR / png_name}")

    pd.DataFrame(tsv_rows).to_csv(
        OUTTSV / png_name.replace(".png", "_source_data.tsv"), sep="\t", index=False
    )


def main() -> None:
    brazil, spain = load_data()

    plot_country(
        country="Brazil",
        df=brazil,
        group_order=BRAZIL_ORDER,
        png_name="Brazil_uViG_taxonomy_barplots_vcontact3.png",
        svg_name="Brazil_uViG_taxonomy_barplots_vcontact3.svg",
        cutoff_map={"class": None, "family": 0.90, "genus": 0.75},
    )

    plot_country(
        country="Spain",
        df=spain,
        group_order=SPAIN_ORDER,
        png_name="Spain_uViG_taxonomy_barplots_vcontact3.png",
        svg_name="Spain_uViG_taxonomy_barplots_vcontact3.svg",
        cutoff_map={"class": None, "family": None, "genus": None},
    )


if __name__ == "__main__":
    main()
