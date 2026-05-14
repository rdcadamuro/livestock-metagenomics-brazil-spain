#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams["svg.fonttype"] = "none"


# USER: set your base path here


BASE = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil")
OUTDIR = BASE / "FIG"
OUTTSV = BASE / "TSV"

BRAZIL_FILE = BASE / "Genomad_Brazil" / "Brazil_vRhyme_phage_supplementary_all.tsv"
SPAIN_FILE = BASE / "Genomad_Spain" / "Spain_vRhyme_final_taxonomy_unified.tsv"
SPAIN_GROUPS = BASE / "VOTUS_COMPARE_BRAZIL_SWINE_ONLY" / "tables" / "groups_swine_class.tsv"

BRAZIL_ORDER = [
    "Broiler with",
    "Broiler without",
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
OTHER_COLOR = "#BDBDBD"


def normalize_taxon(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    bad = {"na", "nan", "none", "unknown", "unclassified", "-", "singleton"}
    if text.lower() in bad:
        return None
    return text


def brazil_group(sample: str) -> str | None:
    match = re.search(r"(\d+)$", str(sample))
    if not match:
        return None
    snum = int(match.group(1)) % 100
    if 35 <= snum <= 41:
        return "Broiler with"
    if snum == 42:
        return "Broiler without"
    if snum == 43:
        return "Biodigestor in NO/A"
    if snum >= 44 and snum % 2 == 0:
        return "Biodigestor in"
    if snum >= 45 and snum % 2 == 1:
        return "Biodigestor out"
    return None


def extract_labeled_rank(text: object, rank: str) -> str | None:
    if pd.isna(text):
        return None
    match = re.search(rf"{rank}:([^;]+)", str(text))
    if not match:
        return None
    return normalize_taxon(match.group(1))


def extract_spain_rank(text: object, rank: str) -> str | None:
    if pd.isna(text):
        return None
    parts = [p.strip() for p in str(text).split(";") if str(p).strip()]
    index_map = {"class": 3, "family": 4, "genus": 5}
    idx = index_map[rank]
    if len(parts) <= idx:
        return None
    return normalize_taxon(parts[idx])


def load_brazil() -> pd.DataFrame:
    df = pd.read_csv(BRAZIL_FILE, sep="\t", low_memory=False).copy()
    df["group"] = df["sample"].map(brazil_group)
    df = df[df["group"].notna()].copy()
    for rank in LEVELS:
        df[rank] = df["phabox2_taxonomy"].map(lambda x, r=rank: extract_labeled_rank(x, r))
    return df


def load_spain() -> pd.DataFrame:
    df = pd.read_csv(SPAIN_FILE, sep="\t", low_memory=False).copy()
    groups = pd.read_csv(SPAIN_GROUPS, sep="\t")
    groups = groups[groups["Sample"].astype(str).str.startswith("SRR")].copy()
    groups["group"] = groups["Group"].astype(str).str.replace("Spain_", "", regex=False)
    df = df.merge(groups[["Sample", "group"]], left_on="sample", right_on="Sample", how="left")
    df = df[df["group"].notna()].copy()
    for rank in LEVELS:
        df[rank] = df["agreed_taxonomy_prefix"].map(lambda x, r=rank: extract_spain_rank(x, r))
    return df


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
        if keep and keep[-1] != totals.index[min(len(keep), len(totals)) - 1]:
            pass
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

    fig, axes = plt.subplots(3, 1, figsize=(15, max(11, 1.0 * len(group_order) + 10)))
    fig.subplots_adjust(left=0.23, right=0.72, top=0.94, bottom=0.14, hspace=0.38)
    fig.suptitle(f"{country} uViG taxonomy by group", fontsize=22, fontweight="bold", y=0.98)

    notes: list[str] = []
    tsv_rows: list[dict[str, object]] = []

    for ax, level in zip(axes, LEVELS):
        pivot, collapse_note = summarize_level(df, group_order, level, cutoff_map.get(level))
        tax_labels = list(pivot.columns)
        colors = palette(tax_labels)
        totals = pivot.sum(axis=1).replace(0, np.nan)
        pct = pivot.div(totals, axis=0).multiply(100).fillna(0)

        y = np.arange(len(group_order))
        left = np.zeros(len(group_order))

        for taxon in tax_labels:
            values = pct[taxon].values
            ax.barh(
                y,
                values,
                left=left,
                color=colors[taxon],
                edgecolor="white",
                linewidth=0.5,
                label=taxon,
            )
            for idx, val in enumerate(values):
                if val >= 7:
                    ax.text(left[idx] + val / 2, y[idx], f"{val:.0f}%", ha="center", va="center",
                            fontsize=9, color="white", fontweight="bold")
            left += values

            for group, count in pivot[taxon].items():
                tsv_rows.append({
                    "country": country,
                    "group": group,
                    "level": level,
                    "taxon": taxon,
                    "count": int(count),
                })

        ax.set_xlim(0, 100)
        ax.set_yticks(y)
        ax.set_yticklabels(group_order, fontsize=11)
        ax.invert_yaxis()
        ax.set_xlabel("Relative abundance within resolved bins (%)", fontsize=11)
        ax.set_title(level.title(), loc="left", fontsize=16, fontweight="bold")
        ax.grid(axis="x", alpha=0.2, linewidth=0.6)
        ax.set_axisbelow(True)

        total, resolved, resolved_pct = resolution_note(df, level)
        if level == "family":
            notes.append(f"* Family represents {resolved_pct:.1f}% of all uViGs ({resolved}/{total}).")
        elif level == "genus":
            notes.append(f"** Genus represents {resolved_pct:.1f}% of all uViGs ({resolved}/{total}).")
        if collapse_note:
            marker = "***" if level == "family" else "****" if level == "genus" else None
            if marker:
                notes.append(f"{marker} {collapse_note}")

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

    fig.text(0.23, 0.04, "\n".join(notes), ha="left", va="bottom", fontsize=11)
    fig.savefig(OUTDIR / png_name, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / svg_name, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(tsv_rows).to_csv(OUTTSV / png_name.replace(".png", "_source_data.tsv"), sep="\t", index=False)


def main() -> None:
    brazil = load_brazil()
    spain = load_spain()

    plot_country(
        country="Brazil",
        df=brazil,
        group_order=BRAZIL_ORDER,
        png_name="Brazil_uViG_taxonomy_barplots.png",
        svg_name="Brazil_uViG_taxonomy_barplots.svg",
        cutoff_map={"class": None, "family": 0.90, "genus": 0.75},
    )

    plot_country(
        country="Spain",
        df=spain,
        group_order=SPAIN_ORDER,
        png_name="Spain_uViG_taxonomy_barplots.png",
        svg_name="Spain_uViG_taxonomy_barplots.svg",
        cutoff_map={"class": None, "family": None, "genus": None},
    )


if __name__ == "__main__":
    main()
