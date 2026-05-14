#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

try:
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import pdist, squareform
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

# USER: set your base path here

BASE_OUT = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil')
TSV_OUT = BASE_OUT / 'TSV'
FIG_OUT = BASE_OUT / 'FIG'

# USER: set your base path here

GTDB_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/GTDBK/Downloaded_ubu_GTDBK')
# USER: set your base path here
GTDB_MASTER = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/GTDBK/MAGs_sample_bin_taxonomy_cultivability.tsv')
# USER: set your base path here
MAG_FUNC = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/ARGs/ALL_SAMPLES_UNIFIED_with_categoria_final.tsv')

OUT_FIG_PNG = FIG_OUT / 'Fig2_MAG_Taxonomy_Function_Brazil.png'
OUT_FIG_SVG = FIG_OUT / 'Fig2_MAG_Taxonomy_Function_Brazil.svg'
OUT_SOURCE = TSV_OUT / 'Fig2_source_data.tsv'
OUT_D_GENUS_LIST = TSV_OUT / 'Fig2_panelD_genera_by_context.tsv'
OUT_D_STACK_TSV = TSV_OUT / 'Fig2_panelD_stacked_genus_composition.tsv'
OUT_PLACEHOLDER_RESOLUTION_TSV = TSV_OUT / 'Fig2_placeholder_taxon_resolution.tsv'
OUT_FAMILY_PLACEHOLDER_RESOLUTION_TSV = TSV_OUT / 'Fig2_family_placeholder_resolution.tsv'
OUT_PANELB_CLUSTER_STATS_TSV = TSV_OUT / 'Fig2_panelB_clustering_stats.tsv'

GROUP_ORDER = [
    'Broiler litter W/A',
    'Broiler litter No/A',
    'Biodigestor Swine Entrance No/A',
    'Biodigestor Swine In',
    'Biodigestor Swine Out',
]
GROUP_COL = {
    'Broiler litter W/A': '#E69F00',
    'Broiler litter No/A': '#C44E00',
    'Biodigestor Swine Entrance No/A': '#009E73',
    'Biodigestor Swine In': '#56B4E9',
    'Biodigestor Swine Out': '#CC79A7',
    'NA': '#999999',
}
TOP_MAGS_HEATMAP = 45
TOP_GENUS_D = 35
TOP_SEGMENTS_PER_GROUP_D = 6
FALLBACK_RANK_ORDER = ['family', 'order', 'class', 'phylum']
N_PERM = 999
N_BOOT = 300
RNG_SEED = 42


def _as_text(v: object) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ''
    return str(v).strip()


def is_unclassified_taxon(v: object) -> bool:
    s = _as_text(v).lower()
    return s in {'', 'unclassified_bacteria', 'unclassified', 'unknown', 'nan', 'none'}


def is_gtdb_placeholder_taxon(name: object) -> bool:
    s = _as_text(name)
    if is_unclassified_taxon(s) or s == 'Remaining genera':
        return False
    if re.search(r'_[A-Z]+$', s):
        return True
    if re.search(r'\d', s):
        return True
    if re.fullmatch(r'[A-Z][A-Z0-9_-]+', s):
        return True
    if s.startswith(('CAG-', 'GCA-', 'GCF-', 'UBA', 'JAG')):
        return True
    return False


def normalize_rank_candidate(v: object) -> str:
    s = _as_text(v)
    if not s:
        return ''
    m = re.match(r'^(.+?)_[A-Z]$', s)
    if m:
        base = _as_text(m.group(1))
        if base and (not is_unclassified_taxon(base)) and (not is_gtdb_placeholder_taxon(base)):
            return base
    return s


def display_taxon(name: object, lineage: dict | None = None) -> str:
    s = _as_text(name)
    if not s:
        return 'Unclassified_bacteria'
    if not is_gtdb_placeholder_taxon(s):
        return s
    ln = lineage or {}
    for rk in FALLBACK_RANK_ORDER:
        cand = normalize_rank_candidate(ln.get(rk, ''))
        if cand and (not is_unclassified_taxon(cand)) and (not is_gtdb_placeholder_taxon(cand)):
            return cand
    return s


def build_genus_lineage_map(df: pd.DataFrame) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    needed = ['genus', 'family', 'order', 'class', 'phylum']
    if df.empty or any(c not in df.columns for c in needed):
        return out
    tmp = df[needed].copy()
    for c in needed:
        tmp[c] = tmp[c].map(_as_text)
    tmp = tmp[~tmp['genus'].map(is_unclassified_taxon)]
    if tmp.empty:
        return out
    for genus, gdf in tmp.groupby('genus', dropna=False):
        rec: dict[str, str] = {}
        for rk in FALLBACK_RANK_ORDER:
            vals = gdf[rk]
            vals = vals[~vals.map(is_unclassified_taxon)]
            if vals.empty:
                rec[rk] = ''
            else:
                rec[rk] = str(vals.value_counts().index[0])
        out[str(genus)] = rec
    return out


def build_family_lineage_map(df: pd.DataFrame) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    needed = ['family', 'order', 'class', 'phylum']
    if df.empty or any(c not in df.columns for c in needed):
        return out
    tmp = df[needed].copy()
    for c in needed:
        tmp[c] = tmp[c].map(_as_text)
    tmp = tmp[~tmp['family'].map(is_unclassified_taxon)]
    if tmp.empty:
        return out
    for fam, gdf in tmp.groupby('family', dropna=False):
        rec: dict[str, str] = {}
        for rk in ['order', 'class', 'phylum']:
            vals = gdf[rk]
            vals = vals[~vals.map(is_unclassified_taxon)]
            rec[rk] = str(vals.value_counts().index[0]) if not vals.empty else ''
        out[str(fam)] = rec
    return out


def resolve_placeholder_label(genus: object, lineage_map: dict[str, dict[str, str]]) -> tuple[str, str]:
    g = _as_text(genus)
    if not g:
        return 'Unclassified_bacteria', 'none'
    if not is_gtdb_placeholder_taxon(g):
        return g, 'genus'
    ln = lineage_map.get(g, {})
    for rk in FALLBACK_RANK_ORDER:
        cand = normalize_rank_candidate(ln.get(rk, ''))
        if cand and (not is_unclassified_taxon(cand)) and (not is_gtdb_placeholder_taxon(cand)):
            return cand, rk
    return g, 'genus_placeholder'


def resolve_family_label(family: object, family_map: dict[str, dict[str, str]]) -> tuple[str, str]:
    f = _as_text(family)
    if not f:
        return 'Unclassified_bacteria', 'none'
    if not is_gtdb_placeholder_taxon(f):
        return f, 'family'
    ln = family_map.get(f, {})
    for rk in ['order', 'class', 'phylum']:
        cand = normalize_rank_candidate(ln.get(rk, ''))
        if cand and (not is_unclassified_taxon(cand)) and (not is_gtdb_placeholder_taxon(cand)):
            return cand, rk
    return f, 'family_placeholder'


def p_to_sig(p: float) -> str:
    if not np.isfinite(p):
        return 'NA'
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    if p < 0.10:
        return '·'
    return 'ns'


def comb2(x: np.ndarray) -> np.ndarray:
    return x * (x - 1) / 2.0


def adjusted_rand_index(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    a = np.asarray(labels_a)
    b = np.asarray(labels_b)
    if a.shape[0] != b.shape[0] or a.shape[0] < 2:
        return np.nan
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    cont = np.zeros((len(ua), len(ub)), dtype=float)
    for i in range(a.shape[0]):
        cont[ia[i], ib[i]] += 1.0
    nij2 = comb2(cont).sum()
    ai = cont.sum(axis=1)
    bj = cont.sum(axis=0)
    ai2 = comb2(ai).sum()
    bj2 = comb2(bj).sum()
    n2 = comb2(np.array([a.shape[0]], dtype=float))[0]
    if n2 <= 0:
        return np.nan
    exp = (ai2 * bj2) / n2
    mx = 0.5 * (ai2 + bj2)
    den = mx - exp
    if den <= 0:
        return np.nan
    return float((nij2 - exp) / den)


def permanova_oneway(distance: np.ndarray, groups: np.ndarray, n_perm: int = N_PERM, seed: int = RNG_SEED) -> dict[str, float]:
    g = np.asarray(groups, dtype=object)
    D = np.asarray(distance, dtype=float)
    n = D.shape[0]
    lev, cnt = np.unique(g, return_counts=True)
    if len(lev) < 2 or np.any(cnt < 2) or n <= len(lev):
        return {'F': np.nan, 'R2': np.nan, 'p': np.nan, 'n': float(n), 'n_groups': float(len(lev))}

    D2 = D ** 2

    def calc_f(lbl: np.ndarray) -> tuple[float, float]:
        lv, ct = np.unique(lbl, return_counts=True)
        if len(lv) < 2 or np.any(ct < 2):
            return np.nan, np.nan
        ss_total = np.triu(D2, 1).sum() / n
        ss_within = 0.0
        for x, nx in zip(lv, ct):
            idx = np.where(lbl == x)[0]
            if nx <= 1:
                continue
            ss_within += np.triu(D2[np.ix_(idx, idx)], 1).sum() / nx
        ss_between = ss_total - ss_within
        dfb = len(lv) - 1
        dfw = n - len(lv)
        if dfb <= 0 or dfw <= 0 or ss_within <= 0 or ss_total <= 0:
            return np.nan, np.nan
        fval = (ss_between / dfb) / (ss_within / dfw)
        r2 = ss_between / ss_total
        return float(fval), float(r2)

    f_obs, r2_obs = calc_f(g)
    if not np.isfinite(f_obs):
        return {'F': np.nan, 'R2': np.nan, 'p': np.nan, 'n': float(n), 'n_groups': float(len(lev))}

    rng = np.random.default_rng(seed)
    ge = 1
    for _ in range(int(n_perm)):
        gp = rng.permutation(g)
        fp, _ = calc_f(gp)
        if np.isfinite(fp) and fp >= f_obs:
            ge += 1
    p = ge / (int(n_perm) + 1)
    return {'F': f_obs, 'R2': r2_obs, 'p': float(p), 'n': float(n), 'n_groups': float(len(lev))}


def permdisp_oneway(distance: np.ndarray, groups: np.ndarray, n_perm: int = N_PERM, seed: int = RNG_SEED) -> dict[str, float]:
    g = np.asarray(groups, dtype=object)
    D = np.asarray(distance, dtype=float)
    n = D.shape[0]
    lev, cnt = np.unique(g, return_counts=True)
    if len(lev) < 2 or np.any(cnt < 2) or n <= len(lev):
        return {'F': np.nan, 'p': np.nan, 'n': float(n), 'n_groups': float(len(lev))}

    # PCoA embedding from distance matrix
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    ord_idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[ord_idx]
    eigvecs = eigvecs[:, ord_idx]
    pos = eigvals > 1e-12
    if not np.any(pos):
        return {'F': np.nan, 'p': np.nan, 'n': float(n), 'n_groups': float(len(lev))}
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])

    def f_anova(lbl: np.ndarray) -> float:
        lv, ct = np.unique(lbl, return_counts=True)
        if len(lv) < 2 or np.any(ct < 2):
            return np.nan
        dists = np.zeros(coords.shape[0], dtype=float)
        for x in lv:
            idx = np.where(lbl == x)[0]
            cen = coords[idx].mean(axis=0)
            dists[idx] = np.sqrt(((coords[idx] - cen) ** 2).sum(axis=1))
        grand = dists.mean()
        ssb = 0.0
        ssw = 0.0
        for x in lv:
            idx = np.where(lbl == x)[0]
            ssb += len(idx) * (dists[idx].mean() - grand) ** 2
            ssw += ((dists[idx] - dists[idx].mean()) ** 2).sum()
        dfb = len(lv) - 1
        dfw = len(lbl) - len(lv)
        if dfb <= 0 or dfw <= 0 or ssw <= 0:
            return np.nan
        return float((ssb / dfb) / (ssw / dfw))

    f_obs = f_anova(g)
    if not np.isfinite(f_obs):
        return {'F': np.nan, 'p': np.nan, 'n': float(n), 'n_groups': float(len(lev))}

    rng = np.random.default_rng(seed + 11)
    ge = 1
    for _ in range(int(n_perm)):
        gp = rng.permutation(g)
        fp = f_anova(gp)
        if np.isfinite(fp) and fp >= f_obs:
            ge += 1
    p = ge / (int(n_perm) + 1)
    return {'F': f_obs, 'p': float(p), 'n': float(n), 'n_groups': float(len(lev))}


def cluster_bootstrap_stability(X: np.ndarray, groups: np.ndarray, n_boot: int = N_BOOT, seed: int = RNG_SEED) -> dict[str, float]:
    X = np.asarray(X, dtype=float)
    g = np.asarray(groups, dtype=object)
    n, p = X.shape
    lev = np.unique(g)
    k = max(2, min(len(lev), n - 1))
    if n < 4 or p < 2 or k < 2:
        return {'ari_mean': np.nan, 'ari_sd': np.nan, 'ari_q025': np.nan, 'ari_q975': np.nan, 'n_boot': float(n_boot)}

    X0 = np.log1p(X)
    z0 = linkage(X0, method='ward', metric='euclidean')
    l0 = fcluster(z0, t=k, criterion='maxclust')

    rng = np.random.default_rng(seed + 23)
    aris: list[float] = []
    for _ in range(int(n_boot)):
        cols = rng.integers(0, p, size=p)
        xb = X0[:, cols]
        zb = linkage(xb, method='ward', metric='euclidean')
        lb = fcluster(zb, t=k, criterion='maxclust')
        ari = adjusted_rand_index(l0, lb)
        if np.isfinite(ari):
            aris.append(float(ari))
    if not aris:
        return {'ari_mean': np.nan, 'ari_sd': np.nan, 'ari_q025': np.nan, 'ari_q975': np.nan, 'n_boot': float(n_boot)}
    arr = np.asarray(aris, dtype=float)
    return {
        'ari_mean': float(np.mean(arr)),
        'ari_sd': float(np.std(arr)),
        'ari_q025': float(np.quantile(arr, 0.025)),
        'ari_q975': float(np.quantile(arr, 0.975)),
        'n_boot': float(len(arr)),
    }


def compute_panelb_stats(mat: pd.DataFrame, groups: pd.Series) -> dict[str, float | str]:
    out: dict[str, float | str] = {
        'distance': 'braycurtis',
        'permanova_F': np.nan,
        'permanova_R2': np.nan,
        'permanova_p': np.nan,
        'permdisp_F': np.nan,
        'permdisp_p': np.nan,
        'bootstrap_ari_mean': np.nan,
        'bootstrap_ari_sd': np.nan,
        'bootstrap_ari_q025': np.nan,
        'bootstrap_ari_q975': np.nan,
        'n_samples': float(mat.shape[0]),
        'n_features': float(mat.shape[1]),
        'n_groups': np.nan,
    }
    if (not HAVE_SCIPY) or mat.empty:
        return out
    grp = groups.reindex(mat.index).fillna('NA').astype(str)
    keep = grp != 'NA'
    grp = grp[keep]
    X = mat.loc[keep].to_numpy(dtype=float)
    if X.shape[0] < 4:
        return out
    lev, cnt = np.unique(grp.to_numpy(), return_counts=True)
    valid_levels = lev[cnt >= 2]
    m = grp.isin(valid_levels)
    grp = grp[m]
    X = X[m.to_numpy()]
    if X.shape[0] < 4 or len(np.unique(grp)) < 2:
        return out
    D = squareform(pdist(X, metric='braycurtis'))
    pman = permanova_oneway(D, grp.to_numpy(), n_perm=N_PERM, seed=RNG_SEED)
    pdisp = permdisp_oneway(D, grp.to_numpy(), n_perm=N_PERM, seed=RNG_SEED)
    stab = cluster_bootstrap_stability(X, grp.to_numpy(), n_boot=N_BOOT, seed=RNG_SEED)
    out.update(
        {
            'permanova_F': pman['F'],
            'permanova_R2': pman['R2'],
            'permanova_p': pman['p'],
            'permdisp_F': pdisp['F'],
            'permdisp_p': pdisp['p'],
            'bootstrap_ari_mean': stab['ari_mean'],
            'bootstrap_ari_sd': stab['ari_sd'],
            'bootstrap_ari_q025': stab['ari_q025'],
            'bootstrap_ari_q975': stab['ari_q975'],
            'n_samples': float(X.shape[0]),
            'n_features': float(X.shape[1]),
            'n_groups': float(len(np.unique(grp))),
        }
    )
    return out


def panelb_stats_text(stats: dict[str, float | str]) -> str:
    p_perm = float(stats.get('permanova_p', np.nan))
    r2 = float(stats.get('permanova_R2', np.nan))
    p_disp = float(stats.get('permdisp_p', np.nan))
    ari = float(stats.get('bootstrap_ari_mean', np.nan))
    return (
        f"PERMANOVA (Bray-Curtis): R²={r2:.3f}, p={p_perm:.3g} [{p_to_sig(p_perm)}]\n"
        f"PERMDISP: p={p_disp:.3g} [{p_to_sig(p_disp)}]\n"
        f"Bootstrap stability (ARI): {ari:.3f}"
    )


def parse_tax(classification: str) -> dict:
    out = {
        'phylum': 'Unclassified_bacteria',
        'class': 'Unclassified_bacteria',
        'order': 'Unclassified_bacteria',
        'family': 'Unclassified_bacteria',
        'genus': 'Unclassified_bacteria',
    }
    if not isinstance(classification, str):
        return out
    for p in classification.split(';'):
        p = p.strip()
        if p.startswith('p__'):
            out['phylum'] = p[3:] or 'Unclassified_bacteria'
        elif p.startswith('c__'):
            out['class'] = p[3:] or 'Unclassified_bacteria'
        elif p.startswith('o__'):
            out['order'] = p[3:] or 'Unclassified_bacteria'
        elif p.startswith('f__'):
            out['family'] = p[3:] or 'Unclassified_bacteria'
        elif p.startswith('g__'):
            out['genus'] = p[3:] or 'Unclassified_bacteria'
    return out


def sample_to_group(sample: str) -> str:
    raw = str(sample)
    if raw.startswith('P19109_'):
        raw = raw.split('_', 1)[1]
    s = int(raw) % 100
    if 35 <= s <= 41:
        return 'Broiler litter W/A'
    if s == 42:
        return 'Broiler litter No/A'
    if s == 43:
        return 'Biodigestor Swine Entrance No/A'
    if s >= 44 and s % 2 == 0:
        return 'Biodigestor Swine In'
    if s >= 45 and s % 2 == 1:
        return 'Biodigestor Swine Out'
    return 'NA'


def load_mag_taxonomy(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample in sorted(master['sample'].unique(), key=lambda x: int(x)):
        fp = GTDB_BASE / f'P19109_{sample}' / 'classify' / 'gtdbtk.bac120.summary.tsv'
        if not fp.exists() or fp.stat().st_size == 0:
            continue
        df = pd.read_csv(fp, sep='\t', low_memory=False)
        if 'user_genome' not in df.columns or 'classification' not in df.columns:
            continue
        for _, r in df.iterrows():
            tx = parse_tax(r['classification'])
            rows.append(
                {
                    'sample': str(sample),
                    'bin': str(r['user_genome']),
                    'phylum': tx['phylum'],
                    'class': tx['class'],
                    'order': tx['order'],
                    'family': tx['family'],
                    'genus': tx['genus'],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    TSV_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(GTDB_MASTER, sep='\t', low_memory=False)
    master['sample'] = master['SampleName'].astype(str).str.replace('P19109_', '', regex=False)
    master['bin'] = master['Bin'].astype(str)

    tax = load_mag_taxonomy(master)
    all_bins = master[['sample', 'bin']].drop_duplicates().copy()
    if not tax.empty:
        all_bins = all_bins.merge(tax, on=['sample', 'bin'], how='left')
    for rk in ['phylum', 'class', 'order', 'family', 'genus']:
        if rk not in all_bins.columns:
            all_bins[rk] = 'Unclassified_bacteria'
        all_bins[rk] = all_bins[rk].fillna('Unclassified_bacteria')
    all_bins['PlotGroup'] = all_bins['sample'].map(sample_to_group)
    all_bins['sample_bin'] = all_bins['sample'].astype(str) + '|' + all_bins['bin'].astype(str)
    genus_lineage_map = build_genus_lineage_map(all_bins)
    family_lineage_map = build_family_lineage_map(all_bins)

    # Manual audit table for placeholder resolution (one row per placeholder genus).
    placeholder_counts = (
        all_bins[all_bins['genus'].map(is_gtdb_placeholder_taxon)]
        .groupby('genus', dropna=False)
        .size()
        .reset_index(name='n_bins')
        .sort_values('n_bins', ascending=False)
    )
    if not placeholder_counts.empty:
        audit_rows = []
        for _, r in placeholder_counts.iterrows():
            g = _as_text(r['genus'])
            ln = genus_lineage_map.get(g, {})
            chosen, chosen_rank = resolve_placeholder_label(g, genus_lineage_map)
            audit_rows.append(
                {
                    'placeholder_genus': g,
                    'chosen_label': chosen,
                    'chosen_rank': chosen_rank,
                    'family': _as_text(ln.get('family', '')),
                    'order': _as_text(ln.get('order', '')),
                    'class': _as_text(ln.get('class', '')),
                    'phylum': _as_text(ln.get('phylum', '')),
                    'n_bins': int(r['n_bins']),
                }
            )
        pd.DataFrame(audit_rows).to_csv(OUT_PLACEHOLDER_RESOLUTION_TSV, sep='\t', index=False)
    else:
        pd.DataFrame(
            columns=['placeholder_genus', 'chosen_label', 'chosen_rank', 'family', 'order', 'class', 'phylum', 'n_bins']
        ).to_csv(OUT_PLACEHOLDER_RESOLUTION_TSV, sep='\t', index=False)

    family_counts = (
        all_bins[all_bins['family'].map(is_gtdb_placeholder_taxon)]
        .groupby('family', dropna=False)
        .size()
        .reset_index(name='n_bins')
        .sort_values('n_bins', ascending=False)
    )
    if not family_counts.empty:
        fam_rows = []
        for _, r in family_counts.iterrows():
            f = _as_text(r['family'])
            chosen, chosen_rank = resolve_family_label(f, family_lineage_map)
            ln = family_lineage_map.get(f, {})
            fam_rows.append(
                {
                    'placeholder_family': f,
                    'chosen_label': chosen,
                    'chosen_rank': chosen_rank,
                    'order': _as_text(ln.get('order', '')),
                    'class': _as_text(ln.get('class', '')),
                    'phylum': _as_text(ln.get('phylum', '')),
                    'n_bins': int(r['n_bins']),
                }
            )
        pd.DataFrame(fam_rows).to_csv(OUT_FAMILY_PLACEHOLDER_RESOLUTION_TSV, sep='\t', index=False)
    else:
        pd.DataFrame(
            columns=['placeholder_family', 'chosen_label', 'chosen_rank', 'order', 'class', 'phylum', 'n_bins']
        ).to_csv(OUT_FAMILY_PLACEHOLDER_RESOLUTION_TSV, sep='\t', index=False)

    func = pd.read_csv(MAG_FUNC, sep='\t', low_memory=False)
    func['sample'] = func['Sample'].astype(str).str.replace('P19109_', '', regex=False)
    func['bin'] = func['MAG_norm'].astype(str)
    func['PlotGroup'] = func['sample'].map(sample_to_group)
    func['category'] = func['CATEGORIA_FINAL'].fillna('Unclassified')
    func['sample_bin'] = func['sample'].astype(str) + '|' + func['bin'].astype(str)

    # A) top family/genus
    fam_top = all_bins['family'].value_counts().head(10)
    gen_top = all_bins['genus'].value_counts().head(10)

    # B) Heatmap MAG x all categories (no unification)
    mag_cat = (
        func.groupby(['sample_bin', 'category'], dropna=False)
        .size()
        .reset_index(name='n_hits')
    )
    mat = mag_cat.pivot(index='sample_bin', columns='category', values='n_hits').fillna(0)
    # keep all categories, sorted by total
    cat_order = mat.sum(axis=0).sort_values(ascending=False).index.tolist()
    mat = mat[cat_order]

    # row metadata
    row_meta = all_bins[['sample_bin', 'PlotGroup', 'genus', 'family']].drop_duplicates().set_index('sample_bin')
    row_meta = row_meta.reindex(mat.index)

    # keep most informative MAGs for readability
    keep_rows = mat.sum(axis=1).sort_values(ascending=False).head(TOP_MAGS_HEATMAP).index.tolist()
    mat = mat.loc[keep_rows]
    row_meta = row_meta.loc[keep_rows]
    panelb_stats = compute_panelb_stats(mat, row_meta['PlotGroup'])
    pd.DataFrame([panelb_stats]).to_csv(OUT_PANELB_CLUSTER_STATS_TSV, sep='\t', index=False)

    z = None
    if HAVE_SCIPY and len(mat) > 2:
        z = linkage(np.log1p(mat.values), method='ward', metric='euclidean')
        order = dendrogram(z, no_plot=True)['leaves']
        mat = mat.iloc[order]
        row_meta = row_meta.iloc[order]

    # C) Detailed categories by context (no domain collapsing)
    ctx = (
        func.groupby(['PlotGroup', 'category'], dropna=False)
        .size()
        .reset_index(name='n_hits')
    )
    ctx = ctx[ctx['PlotGroup'].isin(GROUP_ORDER)].copy()
    ctx_p = ctx.pivot(index='PlotGroup', columns='category', values='n_hits').fillna(0)
    # same category order as panel B
    ctx_p = ctx_p.reindex(index=GROUP_ORDER, columns=cat_order, fill_value=0)

    # D) Genus composition by sample type (stacked; preserves total richness count)
    gctx = (
        all_bins[all_bins['PlotGroup'].isin(GROUP_ORDER)]
        .groupby(['genus', 'PlotGroup'], dropna=False)
        .size()
        .reset_index(name='n_mag_bins')
    )
    gctx = gctx[gctx['genus'] != 'Unclassified_bacteria'].copy()
    top_genus = gctx.groupby('genus')['n_mag_bins'].sum().sort_values(ascending=False).head(TOP_GENUS_D).index.tolist()
    gmat = (
        gctx[gctx['genus'].isin(top_genus)]
        .pivot(index='genus', columns='PlotGroup', values='n_mag_bins')
        .reindex(index=top_genus, columns=GROUP_ORDER)
        .fillna(0)
    )
    genus_rich = (
        gctx[gctx['n_mag_bins'] > 0]
        .groupby('PlotGroup')['genus']
        .nunique()
        .reindex(GROUP_ORDER)
        .fillna(0)
        .astype(int)
    )
    top_genera_by_group = (
        gctx.sort_values(['PlotGroup', 'n_mag_bins'], ascending=[True, False])
        .groupby('PlotGroup')
        .head(12)
        .copy()
    )
    top_genera_by_group.to_csv(OUT_D_GENUS_LIST, sep='\t', index=False)

    d_rows = []
    for grp in GROUP_ORDER:
        sub = gctx[gctx['PlotGroup'] == grp].sort_values('n_mag_bins', ascending=False).copy()
        total = int(sub['n_mag_bins'].sum())
        taken = 0
        seg_rank = 1
        for _, r in sub.head(TOP_SEGMENTS_PER_GROUP_D).iterrows():
            n = int(r['n_mag_bins'])
            taken += n
            d_rows.append({
                'PlotGroup': grp,
                'segment_rank': seg_rank,
                'genus': str(r['genus']),
                'n_mag_bins': n,
                'pct_within_group': (100.0 * n / total) if total > 0 else 0.0,
                'group_total': total,
            })
            seg_rank += 1
        rem = total - taken
        if rem > 0:
            d_rows.append({
                'PlotGroup': grp,
                'segment_rank': seg_rank,
                'genus': 'Remaining genera',
                'n_mag_bins': int(rem),
                'pct_within_group': (100.0 * rem / total) if total > 0 else 0.0,
                'group_total': total,
            })
    d_stack = pd.DataFrame(d_rows)
    d_stack.to_csv(OUT_D_STACK_TSV, sep='\t', index=False)

    # Source data
    src_rows = []
    for k, v in fam_top.items():
        src_rows.append({'panel': 'A_family', 'label': k, 'value': int(v)})
    for k, v in gen_top.items():
        src_rows.append({'panel': 'A_genus', 'label': k, 'value': int(v)})
    for idx in mat.index:
        for col in mat.columns:
            src_rows.append({'panel': 'B_heatmap_mag_category', 'label': f'{idx}|{col}', 'value': float(mat.loc[idx, col])})
    src_rows.extend(
        [
            {'panel': 'B_stats', 'label': 'PERMANOVA_R2', 'value': float(panelb_stats.get('permanova_R2', np.nan))},
            {'panel': 'B_stats', 'label': 'PERMANOVA_p', 'value': float(panelb_stats.get('permanova_p', np.nan))},
            {'panel': 'B_stats', 'label': 'PERMDISP_p', 'value': float(panelb_stats.get('permdisp_p', np.nan))},
            {'panel': 'B_stats', 'label': 'Bootstrap_ARI_mean', 'value': float(panelb_stats.get('bootstrap_ari_mean', np.nan))},
        ]
    )
    for idx in ctx_p.index:
        for col in ctx_p.columns:
            src_rows.append({'panel': 'C_context_category', 'label': f'{idx}|{col}', 'value': float(ctx_p.loc[idx, col])})
    for grp, val in genus_rich.items():
        src_rows.append({'panel': 'D_genus_richness_by_group', 'label': grp, 'value': int(val)})
    for _, r in top_genera_by_group.iterrows():
        src_rows.append({'panel': 'D_top_genera_by_group', 'label': f"{r['PlotGroup']}|{r['genus']}", 'value': int(r['n_mag_bins'])})
    for _, r in d_stack.iterrows():
        src_rows.append({'panel': 'D_stacked_genus', 'label': f"{r['PlotGroup']}|{r['genus']}", 'value': int(r['n_mag_bins'])})
    pd.DataFrame(src_rows).to_csv(OUT_SOURCE, sep='\t', index=False)

    # Plot style
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'svg.fonttype': 'none',
    })

    # Wider canvas + larger inter-panel spacing to keep A/B/C/D clearly separated.
    fig = plt.figure(figsize=(30, 16))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 1.50],
        height_ratios=[1, 1],
        wspace=0.42,
        hspace=0.55,
    )

    # A
    axA = fig.add_subplot(gs[0, 0])
    a_df = pd.concat([
        pd.DataFrame({'taxon': fam_top.index, 'n': fam_top.values, 'type': 'family'}),
        pd.DataFrame({'taxon': gen_top.index, 'n': gen_top.values, 'type': 'genus'}),
    ], ignore_index=True)
    a_df['taxon_disp'] = a_df.apply(
        lambda r: (
            display_taxon(r['taxon'], genus_lineage_map.get(str(r['taxon']), {}))
            if r['type'] == 'genus'
            else resolve_family_label(r['taxon'], family_lineage_map)[0]
        ),
        axis=1,
    )
    a_df['label'] = a_df['type'].str.upper().str[:3] + ': ' + a_df['taxon_disp']
    a_df = a_df.sort_values('n', ascending=True)
    colors = a_df['type'].map({'family': '#4C78A8', 'genus': '#72B7B2'})
    axA.barh(a_df['label'], a_df['n'], color=colors)
    axA.set_title('A) MAG taxonomy distribution (top family/genus)')
    axA.set_xlabel('n MAG bins')
    for t in axA.get_yticklabels():
        t.set_fontstyle('italic')

    # B (dendrogram + context stripe + full-category heatmap)
    gsB = gs[0, 1].subgridspec(1, 3, width_ratios=[0.30, 0.06, 1.0], wspace=0.06)
    axBd = fig.add_subplot(gsB[0, 0])
    axBc = fig.add_subplot(gsB[0, 1])
    axBh = fig.add_subplot(gsB[0, 2])

    if HAVE_SCIPY and z is not None and len(mat) > 2:
        dendrogram(z, orientation='left', no_labels=True, color_threshold=None, above_threshold_color='#666666', ax=axBd)
    axBd.set_xticks([])
    axBd.set_yticks([])
    axBd.set_title('B) MAG functional profile\n(all categories, clustered)', fontsize=11)
    axBd.text(
        0.0,
        -0.08,
        panelb_stats_text(panelb_stats),
        transform=axBd.transAxes,
        ha='left',
        va='top',
        fontsize=7,
        color='#333333',
    )

    ctx_vals = row_meta['PlotGroup'].fillna('NA').map(GROUP_COL).tolist()
    rgb = []
    for h in ctx_vals:
        h = h.lstrip('#')
        rgb.append([int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)])
    axBc.imshow(np.array(rgb).reshape(len(rgb), 1, 3), aspect='auto')
    axBc.set_xticks([])
    axBc.set_yticks([])

    hm_vals = np.log1p(mat.values)
    im = axBh.imshow(hm_vals, aspect='auto', interpolation='nearest')
    axBh.set_xticks(range(len(mat.columns)))
    axBh.set_xticklabels(mat.columns, rotation=70, ha='right', fontsize=7)
    axBh.set_yticks(range(len(mat.index)))
    ylabels = [
        display_taxon(row_meta.iloc[i]['genus'], genus_lineage_map.get(str(row_meta.iloc[i]['genus']), {}))
        for i in range(len(mat.index))
    ]
    axBh.set_yticklabels(ylabels, fontsize=4)
    for t in axBh.get_yticklabels():
        t.set_fontstyle('italic')
    axBh.set_xlabel('Functional category (full detail)')
    axBh.set_ylabel('MAG (clustered)')

    cbar = fig.colorbar(im, ax=axBh, fraction=0.024, pad=0.012)
    cbar.set_label('log(1 + n hits)')

    handles = [Patch(facecolor=GROUP_COL[g], label=g) for g in GROUP_ORDER]
    axBh.legend(handles=handles, title='Sample type', loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8)

    # C (restore detailed stacked by all categories)
    axC = fig.add_subplot(gs[1, 0])
    x = np.arange(len(ctx_p.index))
    bottoms = np.zeros(len(ctx_p.index))
    cmap = plt.cm.get_cmap('tab20', max(len(ctx_p.columns), 1))
    for i, col in enumerate(ctx_p.columns):
        vals = ctx_p[col].values
        axC.bar(x, vals, bottom=bottoms, label=col, color=cmap(i))
        bottoms += vals
    axC.set_xticks(x)
    axC.set_xticklabels(ctx_p.index, rotation=20, ha='right')
    axC.set_ylabel('n hits')
    axC.set_title('C) Functional stratification by context (all categories)')
    axC.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1.0),
        fontsize=6,
        frameon=False,
        borderaxespad=0.0,
    )

    # D (stacked genus composition by sample type)
    axD = fig.add_subplot(gs[1, 1])
    x = np.arange(len(GROUP_ORDER), dtype=float)
    totals = [int(genus_rich.get(g, 0)) for g in GROUP_ORDER]
    # Stable color map across groups for better reading.
    all_seg_genera = [g for g in d_stack['genus'].drop_duplicates().tolist() if g != 'Remaining genera']
    cmap_d = plt.cm.get_cmap('tab20', max(len(all_seg_genera), 1))
    genus_color = {g: cmap_d(i) for i, g in enumerate(all_seg_genera)}
    genus_color['Remaining genera'] = '#DDDDDD'

    for i, grp in enumerate(GROUP_ORDER):
        sub = d_stack[d_stack['PlotGroup'] == grp].sort_values('segment_rank')
        bottom = 0.0
        for _, r in sub.iterrows():
            gen = r['genus']
            n = float(r['n_mag_bins'])
            pct = float(r['pct_within_group'])
            axD.bar(
                i,
                n,
                bottom=bottom,
                color=genus_color.get(gen, '#AAAAAA'),
                edgecolor='white',
                linewidth=0.5,
                width=0.8,
            )
            if pct >= 8.0:
                axD.text(
                    i,
                    bottom + n / 2.0,
                    f'{pct:.0f}%',
                    ha='center',
                    va='center',
                    fontsize=8,
                    color='black',
                    fontweight='bold',
                )
            bottom += n

    # total labels on top
    for i, total in enumerate(totals):
        axD.text(i, total + 1.2, str(total), ha='center', va='bottom', fontsize=10, fontweight='bold')

    # top genera labels above each bar with segment color
    labeled = []
    ymax = max(totals) if totals else 1
    for i, grp in enumerate(GROUP_ORDER):
        top3 = d_stack[(d_stack['PlotGroup'] == grp) & (d_stack['genus'] != 'Remaining genera')].sort_values('n_mag_bins', ascending=False).head(3)
        y0 = totals[i] + ymax * 0.06
        for j, (_, r) in enumerate(top3.iterrows()):
            g = str(r['genus'])
            labeled.append(g)
            g_disp = display_taxon(g, genus_lineage_map.get(g, {}))
            axD.text(
                i,
                y0 + j * ymax * 0.045,
                g_disp,
                ha='center',
                va='bottom',
                fontsize=7,
                color=genus_color.get(g, '#333333'),
                style='italic',
            )

    axD.set_xticks(x)
    axD.set_xticklabels(GROUP_ORDER, rotation=20, ha='right', fontsize=9)
    axD.set_ylabel('n unique bacterial genera')
    axD.set_xlabel('Sample type')
    axD.set_title('D) Genus composition within each sample type (stacked)')
    axD.set_ylim(0, ymax * 1.35)

    # compact legend: genera labeled above bars + remaining
    legend_items = list(dict.fromkeys(labeled))
    if 'Remaining genera' in d_stack['genus'].values:
        legend_items.append('Remaining genera')
    handles = [
        Patch(
            facecolor=genus_color[g],
            label=display_taxon(g, genus_lineage_map.get(g, {})) if g != 'Remaining genera' else g,
        )
        for g in legend_items[:20]
    ]
    axD.legend(
        handles=handles,
        title='Top genera segments',
        loc='upper center',
        bbox_to_anchor=(0.5, 1.22),
        ncol=4,
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )

    fig.suptitle('Figure 2 - MAG taxonomy and full functional structure (Brazil-only)', fontsize=15)
    fig.subplots_adjust(left=0.06, right=0.95, top=0.91, bottom=0.08)
    fig.savefig(OUT_FIG_PNG, dpi=450, bbox_inches='tight')
    fig.savefig(OUT_FIG_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    print(f'[OK] figure: {OUT_FIG_PNG}')
    print(f'[OK] figure: {OUT_FIG_SVG}')
    print(f'[OK] source: {OUT_SOURCE}')
    print(f'[OK] panelB stats: {OUT_PANELB_CLUSTER_STATS_TSV}')


if __name__ == '__main__':
    main()
