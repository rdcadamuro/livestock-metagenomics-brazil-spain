#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import itertools
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform

# Inputs
# USER: set your base path here
BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/OneDrive_1_21-02-2026')
PADLOC_TSV = BASE / 'PADLOC_MASTER_COUNTS.tsv'
CRISPR_TSV = BASE / 'CRISPR_MASTER_SUMMARY.tsv'
STRICT_CLASS_TSV = BASE / 'Spacers_All' / 'Results_Taxon_Spacers' / 'anti_phage_system_classification.tsv'

# Outputs
OUT = BASE / 'Spacers_All' / 'Results_Taxon_Spacers'
IMG = OUT / 'IMG_CLEAN'

OUT_BIN = OUT / 'anti_phage_context_bin_table.tsv'
OUT_INTRA = OUT / 'anti_phage_intrasample_heterogeneity.tsv'
OUT_CONTEXT = OUT / 'anti_phage_context_stats.tsv'
OUT_PERM = OUT / 'anti_phage_context_permanova.tsv'
OUT_PAIR = OUT / 'anti_phage_biodigester_paired_medians.tsv'
OUT_REPORT = OUT / 'anti_phage_intrasample_context_report.txt'

FIG_GROUP_PNG = IMG / 'Fig_Context_Load_ByGroup_CompareScriptColors.png'
FIG_GROUP_SVG = IMG / 'Fig_Context_Load_ByGroup_CompareScriptColors.svg'
FIG_CONTEXT2_PNG = IMG / 'Fig_Context_Load_Broiler_vs_Biodigestor_Strict.png'
FIG_CONTEXT2_SVG = IMG / 'Fig_Context_Load_Broiler_vs_Biodigestor_Strict.svg'
FIG_PAIR_PNG = IMG / 'Fig_Biodigestor_InOut_Paired_Medians.png'
FIG_PAIR_SVG = IMG / 'Fig_Biodigestor_InOut_Paired_Medians.svg'
FIG_INTRA_PNG = IMG / 'Fig_IntraSample_Heterogeneity_Mags.png'
FIG_INTRA_SVG = IMG / 'Fig_IntraSample_Heterogeneity_Mags.svg'

DPI = 450
N_PERM = 9999
N_NULL_INTRA = 5000
RNG_SEED = 42

# Palette and groups copied from reference script
GROUP_ORDER = [
    'Chicken W/A',
    'Chicken No/A',
    'Biodigestor Swine Entrance No/A',
    'Biodigestor Swine In',
    'Biodigestor Swine Out',
]
GROUP_COL = {
    'Chicken W/A': '#E69F00',
    'Chicken No/A': '#C44E00',
    'Biodigestor Swine Entrance No/A': '#009E73',
    'Biodigestor Swine In': '#56B4E9',
    'Biodigestor Swine Out': '#CC79A7',
    'NA': '#999999',
}
GROUP_LABEL = {
    'Chicken W/A': 'Broiler litter W/A',
    'Chicken No/A': 'Broiler litter No/A',
    'Biodigestor Swine Entrance No/A': 'Biodigestor Swine Entrance No/A',
    'Biodigestor Swine In': 'Biodigestor Swine In',
    'Biodigestor Swine Out': 'Biodigestor Swine Out',
}
CONTEXT_ORDER = ['Broiler litter', 'Biodigestor']
CONTEXT_COL = {
    'Broiler litter': '#E69F00',
    'Biodigestor': '#009E73',
}


def canon_sample(x: object) -> str:
    s = str(x).strip()
    if s.startswith('P19109_'):
        return s.split('_', 1)[1]
    return s


def canon_bin(x: object) -> str:
    s = str(x).strip()
    if s.startswith('bin.'):
        return s
    if s.startswith('bin_'):
        return s.replace('bin_', 'bin.')
    if s.startswith('vRhyme_bin_'):
        return f"bin.{s.replace('vRhyme_bin_', '')}"
    return s


def sample_snum(sample_id: str) -> Optional[int]:
    try:
        n = int(str(sample_id))
        return n % 100
    except Exception:
        return None


def map_group(sample_id: str) -> str:
    s = sample_snum(sample_id)
    if s is None:
        return 'NA'
    if 35 <= s <= 41:
        return 'Chicken W/A'
    if s == 42:
        return 'Chicken No/A'
    if s == 43:
        return 'Biodigestor Swine Entrance No/A'
    if s >= 44 and s % 2 == 0:
        return 'Biodigestor Swine In'
    if s >= 45 and s % 2 == 1:
        return 'Biodigestor Swine Out'
    return 'NA'


def map_host(sample_id: str) -> str:
    s = sample_snum(sample_id)
    if s is None:
        return 'NA'
    return 'CHICKEN' if s <= 42 else 'SWINE'


def map_stage(sample_id: str) -> str:
    s = sample_snum(sample_id)
    if s is None:
        return 'NA'
    if s == 43:
        return 'IN'
    if s >= 44 and s % 2 == 0:
        return 'IN'
    if s >= 45 and s % 2 == 1:
        return 'OUT'
    return 'NA'


def map_digester(sample_id: str) -> str:
    s = sample_snum(sample_id)
    if s is None:
        return 'NA'
    if s == 43:
        return 'BD_NOAB'
    if s >= 44 and s % 2 == 0:
        pair = (s - 44) // 2 + 1
        return f'BD{pair:02d}'
    if s >= 45 and s % 2 == 1:
        pair = (s - 45) // 2 + 1
        return f'BD{pair:02d}'
    return 'NA'


def map_context(sample_id: str) -> str:
    g = map_group(sample_id)
    if g == 'Chicken W/A':
        return 'Broiler litter'
    if g.startswith('Biodigestor'):
        return 'Biodigestor'
    return 'Other'


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


def braycurtis_safe(u: np.ndarray, v: np.ndarray) -> float:
    den = float(np.sum(u) + np.sum(v))
    if den <= 0:
        return 0.0
    num = float(np.sum(np.abs(u - v)))
    return num / den


def dunn_test(df: pd.DataFrame, value_col: str, group_col: str) -> pd.DataFrame:
    work = df[[value_col, group_col]].dropna().copy()
    work[value_col] = pd.to_numeric(work[value_col], errors='coerce')
    work = work.dropna(subset=[value_col])
    if work.empty:
        return pd.DataFrame(columns=['group1', 'group2', 'z', 'p_raw', 'p_bh'])

    groups = sorted(work[group_col].unique())
    n = len(work)
    if len(groups) < 2:
        return pd.DataFrame(columns=['group1', 'group2', 'z', 'p_raw', 'p_bh'])

    ranks = stats.rankdata(work[value_col].values)
    work = work.assign(_rank=ranks)
    tie_counts = pd.Series(work[value_col].values).value_counts().values
    tie_term = np.sum(tie_counts ** 3 - tie_counts)
    tie_corr = 1.0 - tie_term / (n ** 3 - n) if n > 1 else 1.0
    tie_corr = 1.0 if tie_corr <= 0 else tie_corr

    grp = work.groupby(group_col)['_rank'].agg(['mean', 'count']).rename(columns={'mean': 'rbar', 'count': 'ni'})
    c = n * (n + 1) / 12.0

    rows = []
    for g1, g2 in itertools.combinations(groups, 2):
        r1, r2 = grp.loc[g1, 'rbar'], grp.loc[g2, 'rbar']
        n1, n2 = grp.loc[g1, 'ni'], grp.loc[g2, 'ni']
        se = math.sqrt(c * tie_corr * (1.0 / n1 + 1.0 / n2))
        z = (r1 - r2) / se if se > 0 else 0.0
        p = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
        rows.append((g1, g2, z, p))

    out = pd.DataFrame(rows, columns=['group1', 'group2', 'z', 'p_raw'])
    out['p_bh'] = bh_fdr(out['p_raw'].values)
    return out


def _permute_labels_within_strata(labels: np.ndarray, strata: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    if strata is None:
        return rng.permutation(labels)
    out = labels.copy()
    for st in np.unique(strata):
        mask = strata == st
        out[mask] = rng.permutation(out[mask])
    return out


def permanova_oneway(dmat: pd.DataFrame, group: pd.Series, n_perm: int = 999, seed: int = 42, strata: pd.Series | None = None) -> dict:
    samples = dmat.index.astype(str)
    group = group.reindex(samples)
    if strata is not None:
        strata = strata.reindex(samples)

    ok = group.notna()
    if strata is not None:
        ok = ok & strata.notna()

    samples = samples[ok.to_numpy()]
    d = dmat.loc[samples, samples].to_numpy(dtype=float)
    g = group.loc[samples].astype(str).to_numpy()
    st = strata.loc[samples].astype(str).to_numpy() if strata is not None else None

    n = d.shape[0]
    uniq, counts = np.unique(g, return_counts=True)
    k = uniq.size
    if n < 3 or k < 2 or np.any(counts < 2):
        return {'n': n, 'k': k, 'F': np.nan, 'R2': np.nan, 'p_perm': np.nan}

    a = -0.5 * (d ** 2)
    j = np.eye(n) - np.ones((n, n)) / n
    b = j @ a @ j
    ss_total = float(np.trace(b))

    def ss_among(labels):
        ss = 0.0
        for lab in np.unique(labels):
            idx = np.where(labels == lab)[0]
            ng = idx.size
            bg = b[np.ix_(idx, idx)]
            ss += float(bg.sum()) / ng
        return ss

    ss_a = ss_among(g)
    ss_w = ss_total - ss_a
    df_a = float(k - 1)
    df_w = float(n - k)
    ms_a = ss_a / df_a
    ms_w = ss_w / df_w if df_w > 0 else np.nan
    f_obs = ms_a / ms_w if (ms_w is not None and ms_w > 0) else np.nan
    r2 = ss_a / ss_total if ss_total > 0 else np.nan
    if not np.isfinite(f_obs):
        return {'n': n, 'k': k, 'F': float(f_obs), 'R2': float(r2), 'p_perm': np.nan}

    rng = np.random.default_rng(int(seed))
    ge = 0
    for _ in range(int(n_perm)):
        gp = _permute_labels_within_strata(g, st, rng=rng)
        ss_ap = ss_among(gp)
        ss_wp = ss_total - ss_ap
        ms_ap = ss_ap / df_a
        ms_wp = ss_wp / df_w if df_w > 0 else np.nan
        fp = ms_ap / ms_wp if (ms_wp is not None and ms_wp > 0) else np.nan
        if np.isfinite(fp) and np.isfinite(f_obs) and fp >= f_obs:
            ge += 1

    p = (ge + 1) / (int(n_perm) + 1)
    return {'n': n, 'k': k, 'F': float(f_obs), 'R2': float(r2), 'p_perm': float(p)}


def permdisp_oneway(dmat: pd.DataFrame, group: pd.Series, n_perm: int = 999, seed: int = 42, strata: pd.Series | None = None) -> dict:
    # PCoA coordinates from distance matrix
    labels = dmat.index.astype(str)
    group = group.reindex(labels)
    if strata is not None:
        strata = strata.reindex(labels)

    ok = group.notna()
    if strata is not None:
        ok = ok & strata.notna()
    labels = labels[ok.to_numpy()]
    d = dmat.loc[labels, labels].to_numpy(dtype=float)
    g = group.loc[labels].astype(str).to_numpy()
    st = strata.loc[labels].astype(str).to_numpy() if strata is not None else None

    n = len(labels)
    if n < 3 or len(np.unique(g)) < 2:
        return {'n': n, 'k': len(np.unique(g)), 'F': np.nan, 'p_perm': np.nan}

    d2 = d ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ d2 @ J
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    vals = vals[idx]
    vecs = vecs[:, idx]
    pos = vals > 1e-12
    coords = vecs[:, pos] * np.sqrt(vals[pos]) if np.any(pos) else np.zeros((n, 1))

    uniq = np.unique(g)
    if np.any([(g == u).sum() < 2 for u in uniq]):
        return {'n': n, 'k': len(uniq), 'F': np.nan, 'p_perm': np.nan}

    def dist_to_centroid(labels_arr):
        out = np.zeros(n, dtype=float)
        for u in np.unique(labels_arr):
            ii = np.where(labels_arr == u)[0]
            c = coords[ii].mean(axis=0)
            out[ii] = np.sqrt(((coords[ii] - c) ** 2).sum(axis=1))
        return out

    def anova_f(values, labels_arr):
        ddf = pd.DataFrame({'v': values, 'g': labels_arr})
        gs = ddf.groupby('g')['v'].agg(['count', 'mean'])
        if len(gs) < 2:
            return np.nan
        n_tot = len(ddf)
        k = len(gs)
        if n_tot <= k:
            return np.nan
        grand = ddf['v'].mean()
        ss_b = float((gs['count'] * (gs['mean'] - grand) ** 2).sum())
        merged = ddf.join(gs['mean'], on='g', rsuffix='_g')
        ss_w = float(((merged['v'] - merged['mean']) ** 2).sum())
        df_b = k - 1
        df_w = n_tot - k
        if df_w <= 0:
            return np.nan
        if ss_w == 0:
            return 0.0 if ss_b == 0 else np.inf
        return (ss_b / df_b) / (ss_w / df_w)

    obs = anova_f(dist_to_centroid(g), g)
    if not np.isfinite(obs):
        return {'n': n, 'k': len(uniq), 'F': float(obs), 'p_perm': np.nan}
    rng = np.random.default_rng(int(seed))
    ge = 0
    for _ in range(int(n_perm)):
        gp = _permute_labels_within_strata(g, st, rng=rng)
        fp = anova_f(dist_to_centroid(gp), gp)
        if np.isfinite(fp) and fp >= obs:
            ge += 1
    p = (ge + 1) / (int(n_perm) + 1)
    return {'n': n, 'k': len(uniq), 'F': float(obs), 'p_perm': float(p)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    IMG.mkdir(parents=True, exist_ok=True)

    pad = pd.read_csv(PADLOC_TSV, sep='\t', dtype=str)
    cri = pd.read_csv(CRISPR_TSV, sep='\t', dtype=str)
    cls = pd.read_csv(STRICT_CLASS_TSV, sep='\t', dtype=str)

    strict_set = set(cls.loc[cls['classification'] == 'anti_phage_established', 'system'].astype(str))

    # normalize padloc
    pad = pad.copy()
    pad['sample'] = pad['sample'].map(canon_sample)
    pad['bin'] = pad['bin'].map(canon_bin)
    pad['count'] = pd.to_numeric(pad['count'], errors='coerce').fillna(0.0)
    pad = pad[pad['system'].isin(strict_set)].copy()

    # normalize crispr
    sample_col = 'Sample' if 'Sample' in cri.columns else 'sample'
    bin_col = 'BIN' if 'BIN' in cri.columns else 'bin'
    nsp_col = 'n_spacers' if 'n_spacers' in cri.columns else 'nspacers'
    cri = cri.copy()
    cri['sample'] = cri[sample_col].map(canon_sample)
    cri['bin'] = cri[bin_col].map(canon_bin)
    cri['n_spacers'] = pd.to_numeric(cri[nsp_col], errors='coerce').fillna(0.0)
    cri['CRISPR_active_bin'] = (cri['n_spacers'] > 0).astype(int)

    # bin x system matrix
    mat = pad.pivot_table(index=['sample', 'bin'], columns='system', values='count', aggfunc='sum', fill_value=0.0)
    cri_bin = cri.groupby(['sample', 'bin'], as_index=True)['CRISPR_active_bin'].max().astype(float)
    mat = mat.join(cri_bin, how='outer').fillna(0.0)

    # bin table
    bin_df = mat.reset_index().copy()
    feature_cols = [c for c in bin_df.columns if c not in ['sample', 'bin']]
    bin_df['anti_phage_padloc_total_bin'] = bin_df[[c for c in feature_cols if c != 'CRISPR_active_bin']].sum(axis=1)
    bin_df['anti_phage_load_bin'] = bin_df[feature_cols].sum(axis=1)
    bin_df['PlotGroup'] = bin_df['sample'].map(map_group)
    bin_df['Host'] = bin_df['sample'].map(map_host)
    bin_df['Stage'] = bin_df['sample'].map(map_stage)
    bin_df['Digester'] = bin_df['sample'].map(map_digester)
    bin_df['Context'] = bin_df['sample'].map(map_context)
    bin_df = bin_df.sort_values(['sample', 'bin']).reset_index(drop=True)
    bin_df.to_csv(OUT_BIN, sep='\t', index=False)

    # ------- Intra-sample heterogeneity among MAGs -------
    # relative profiles
    X = mat.copy()
    rs = X.sum(axis=1)
    Xrel = X.div(rs.replace(0, np.nan), axis=0).fillna(0.0)
    Xrel_nz = Xrel.loc[rs > 0].copy()

    # full distance matrix once
    labels = Xrel_nz.index.to_list()
    D = squareform(pdist(Xrel_nz.values, metric=braycurtis_safe))
    Ddf = pd.DataFrame(D, index=pd.Index(labels), columns=pd.Index(labels))

    rng = np.random.default_rng(RNG_SEED)
    all_idx = np.arange(len(labels))

    rows = []
    sample_to_indices = {}
    for i, (s, b) in enumerate(labels):
        sample_to_indices.setdefault(str(s), []).append(i)

    for s in sorted(bin_df['sample'].astype(str).unique(), key=int):
        idxs = sample_to_indices.get(str(s), [])
        k = len(idxs)
        if k < 3:
            rows.append({'sample': s, 'n_bins': k, 'mean_pairwise_bray': np.nan, 'p_high_heterogeneity': np.nan, 'p_low_heterogeneity': np.nan, 'p_two_sided': np.nan, 'zscore_vs_null': np.nan})
            continue

        sub = D[np.ix_(idxs, idxs)]
        iu = np.triu_indices(k, 1)
        obs = float(np.mean(sub[iu]))

        null_vals = np.empty(N_NULL_INTRA, dtype=float)
        for t in range(N_NULL_INTRA):
            pick = rng.choice(all_idx, size=k, replace=False)
            ss = D[np.ix_(pick, pick)]
            null_vals[t] = float(np.mean(ss[np.triu_indices(k, 1)]))

        p_hi = (np.sum(null_vals >= obs) + 1) / (N_NULL_INTRA + 1)
        p_lo = (np.sum(null_vals <= obs) + 1) / (N_NULL_INTRA + 1)
        p_two = min(1.0, 2 * min(p_hi, p_lo))
        z = (obs - float(np.mean(null_vals))) / (float(np.std(null_vals, ddof=1)) if float(np.std(null_vals, ddof=1)) > 0 else np.nan)

        rows.append({'sample': s, 'n_bins': k, 'mean_pairwise_bray': obs, 'p_high_heterogeneity': p_hi, 'p_low_heterogeneity': p_lo, 'p_two_sided': p_two, 'zscore_vs_null': z})

    intra = pd.DataFrame(rows).sort_values('sample')
    intra['q_two_sided_bh'] = np.nan
    m_ok = intra['p_two_sided'].notna()
    if m_ok.any():
        intra.loc[m_ok, 'q_two_sided_bh'] = bh_fdr(intra.loc[m_ok, 'p_two_sided'].values)
    intra['PlotGroup'] = intra['sample'].map(map_group)
    intra.to_csv(OUT_INTRA, sep='\t', index=False)

    # ------- Context stats -------
    stats_rows = []

    # global across 5 plot groups
    dat = bin_df[bin_df['PlotGroup'].isin(GROUP_ORDER)].copy()
    groups = [g['anti_phage_load_bin'].values for _, g in dat.groupby('PlotGroup')]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        kw = stats.kruskal(*groups)
        stats_rows.append({'analysis': 'anti_phage_load_by_plotgroup', 'test': 'Kruskal-Wallis', 'statistic': float(kw.statistic), 'p_value': float(kw.pvalue), 'detail': 'Groups from compare script'})
    dunn = dunn_test(dat, 'anti_phage_load_bin', 'PlotGroup')
    if not dunn.empty:
        for _, r in dunn.iterrows():
            stats_rows.append({'analysis': 'anti_phage_load_by_plotgroup', 'test': 'Dunn_BH', 'statistic': float(r['z']), 'p_value': float(r['p_bh']), 'detail': f"{r['group1']} vs {r['group2']} (raw_p={r['p_raw']})"})

    # context: Broiler litter vs Biodigestor (exclude Other)
    ctx = bin_df[bin_df['Context'].isin(['Broiler litter', 'Biodigestor'])].copy()
    a = ctx.loc[ctx['Context'] == 'Broiler litter', 'anti_phage_load_bin'].values
    b = ctx.loc[ctx['Context'] == 'Biodigestor', 'anti_phage_load_bin'].values
    if len(a) >= 2 and len(b) >= 2:
        mw = stats.mannwhitneyu(a, b, alternative='two-sided')
        stats_rows.append({'analysis': 'broiler_vs_biodigestor_load', 'test': 'Mann-Whitney', 'statistic': float(mw.statistic), 'p_value': float(mw.pvalue), 'detail': 'bin-level anti-phage load'})

    # biodigester IN vs OUT (bin-level)
    bio = bin_df[bin_df['PlotGroup'].isin(['Biodigestor Swine In', 'Biodigestor Swine Out'])].copy()
    ain = bio.loc[bio['PlotGroup'] == 'Biodigestor Swine In', 'anti_phage_load_bin'].values
    aout = bio.loc[bio['PlotGroup'] == 'Biodigestor Swine Out', 'anti_phage_load_bin'].values
    if len(ain) >= 2 and len(aout) >= 2:
        mw2 = stats.mannwhitneyu(ain, aout, alternative='two-sided')
        stats_rows.append({'analysis': 'biodigestor_in_vs_out_load', 'test': 'Mann-Whitney', 'statistic': float(mw2.statistic), 'p_value': float(mw2.pvalue), 'detail': 'bin-level anti-phage load'})

    # paired digesters: compare sample-level medians IN vs OUT
    smed = bin_df.groupby('sample', as_index=False)['anti_phage_load_bin'].median().rename(columns={'anti_phage_load_bin': 'median_load'})
    smed['Stage'] = smed['sample'].map(map_stage)
    smed['Digester'] = smed['sample'].map(map_digester)
    pair = smed[smed['Digester'].str.startswith('BD', na=False) & smed['Digester'].str.match(r'BD\d\d')].copy()
    in_tab = pair[pair['Stage'] == 'IN'][['Digester', 'sample', 'median_load']].rename(columns={'sample': 'sample_in', 'median_load': 'median_in'})
    out_tab = pair[pair['Stage'] == 'OUT'][['Digester', 'sample', 'median_load']].rename(columns={'sample': 'sample_out', 'median_load': 'median_out'})
    paired = in_tab.merge(out_tab, on='Digester', how='inner').sort_values('Digester')
    if not paired.empty:
        paired['delta_out_minus_in'] = paired['median_out'] - paired['median_in']
        paired.to_csv(OUT_PAIR, sep='\t', index=False)
        if len(paired) >= 2:
            wx = stats.wilcoxon(paired['median_in'], paired['median_out'], alternative='two-sided', zero_method='wilcox')
            stats_rows.append({'analysis': 'biodigestor_paired_median_load', 'test': 'Wilcoxon paired', 'statistic': float(wx.statistic), 'p_value': float(wx.pvalue), 'detail': f'n_pairs={len(paired)}'})
    else:
        paired = pd.DataFrame(columns=['Digester', 'sample_in', 'median_in', 'sample_out', 'median_out', 'delta_out_minus_in'])
        paired.to_csv(OUT_PAIR, sep='\t', index=False)

    # composition stats via PERMANOVA/PERMDISP
    Xcomp = Xrel_nz.copy()
    idx_df = pd.DataFrame(Xcomp.index.tolist(), columns=['sample', 'bin'])
    idx_df['sample'] = idx_df['sample'].astype(str)
    idx_df['PlotGroup'] = idx_df['sample'].map(map_group)
    idx_df['Context'] = idx_df['sample'].map(map_context)

    Dcomp = pd.DataFrame(squareform(pdist(Xcomp.values, metric=braycurtis_safe)), index=Xcomp.index, columns=Xcomp.index)
    # convert index for reindexing by labels as string composite
    idx_labels = pd.Index([f"{s}|{b}" for s, b in Xcomp.index])
    Dcomp.index = idx_labels
    Dcomp.columns = idx_labels

    md = idx_df.copy()
    md['id'] = [f"{s}|{b}" for s, b in Xcomp.index]
    md = md.set_index('id')

    perm_rows = []
    # 5 groups (compare-script groups)
    g = md['PlotGroup']
    keep = g.isin(GROUP_ORDER)
    dsub = Dcomp.loc[keep.index[keep], keep.index[keep]]
    gsub = g.loc[dsub.index]
    if gsub.nunique() >= 2:
        p1 = permanova_oneway(dsub, gsub, n_perm=N_PERM, seed=RNG_SEED)
        p2 = permdisp_oneway(dsub, gsub, n_perm=N_PERM, seed=RNG_SEED)
        perm_rows.append({'analysis': 'composition_by_plotgroup', 'test': 'PERMANOVA', 'F': p1['F'], 'R2': p1['R2'], 'p_value': p1['p_perm'], 'n': p1['n'], 'k': p1['k']})
        perm_rows.append({'analysis': 'composition_by_plotgroup', 'test': 'PERMDISP', 'F': p2['F'], 'R2': np.nan, 'p_value': p2['p_perm'], 'n': p2['n'], 'k': p2['k']})

    # broiler vs biodigestor composition
    gc = md['Context']
    keep2 = gc.isin(['Broiler litter', 'Biodigestor'])
    dsub2 = Dcomp.loc[keep2.index[keep2], keep2.index[keep2]]
    gsub2 = gc.loc[dsub2.index]
    if gsub2.nunique() >= 2:
        p1 = permanova_oneway(dsub2, gsub2, n_perm=N_PERM, seed=RNG_SEED)
        p2 = permdisp_oneway(dsub2, gsub2, n_perm=N_PERM, seed=RNG_SEED)
        perm_rows.append({'analysis': 'composition_broiler_vs_biodigestor', 'test': 'PERMANOVA', 'F': p1['F'], 'R2': p1['R2'], 'p_value': p1['p_perm'], 'n': p1['n'], 'k': p1['k']})
        perm_rows.append({'analysis': 'composition_broiler_vs_biodigestor', 'test': 'PERMDISP', 'F': p2['F'], 'R2': np.nan, 'p_value': p2['p_perm'], 'n': p2['n'], 'k': p2['k']})

    # biodigestor in vs out composition
    gb = md['PlotGroup']
    keep3 = gb.isin(['Biodigestor Swine In', 'Biodigestor Swine Out'])
    dsub3 = Dcomp.loc[keep3.index[keep3], keep3.index[keep3]]
    gsub3 = gb.loc[dsub3.index]
    if gsub3.nunique() >= 2:
        p1 = permanova_oneway(dsub3, gsub3, n_perm=N_PERM, seed=RNG_SEED)
        p2 = permdisp_oneway(dsub3, gsub3, n_perm=N_PERM, seed=RNG_SEED)
        perm_rows.append({'analysis': 'composition_biodigestor_in_vs_out', 'test': 'PERMANOVA', 'F': p1['F'], 'R2': p1['R2'], 'p_value': p1['p_perm'], 'n': p1['n'], 'k': p1['k']})
        perm_rows.append({'analysis': 'composition_biodigestor_in_vs_out', 'test': 'PERMDISP', 'F': p2['F'], 'R2': np.nan, 'p_value': p2['p_perm'], 'n': p2['n'], 'k': p2['k']})

    context_stats = pd.DataFrame(stats_rows)
    if not context_stats.empty:
        context_stats = context_stats.sort_values(['analysis', 'test']).reset_index(drop=True)
    context_stats.to_csv(OUT_CONTEXT, sep='\t', index=False)

    perm_stats = pd.DataFrame(perm_rows)
    if not perm_stats.empty:
        perm_stats = perm_stats.sort_values(['analysis', 'test']).reset_index(drop=True)
    perm_stats.to_csv(OUT_PERM, sep='\t', index=False)

    # ------- Plots -------
    sns.set_theme(style='whitegrid', context='paper')
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'svg.fonttype': 'none',
    })

    # 1) context boxplot (main, 2 groups)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    datc = bin_df[bin_df['Context'].isin(CONTEXT_ORDER)].copy()
    order_ctx = [g for g in CONTEXT_ORDER if (datc['Context'] == g).any()]
    palette_ctx = {g: CONTEXT_COL[g] for g in order_ctx}
    sns.boxplot(data=datc, x='Context', y='anti_phage_load_bin', order=order_ctx, palette=palette_ctx, ax=ax, fliersize=2)
    sns.stripplot(data=datc, x='Context', y='anti_phage_load_bin', order=order_ctx, ax=ax, color='#111111', size=2, alpha=0.45, jitter=0.22)
    ax.set_title('Anti-phage load per MAG by context (strict anti-phage systems)')
    ax.set_xlabel('')
    ax.set_ylabel('Anti-phage load per MAG (strict PADLOC + CRISPR_active)')
    plt.tight_layout()
    fig.savefig(FIG_CONTEXT2_PNG, dpi=DPI, bbox_inches='tight')
    fig.savefig(FIG_CONTEXT2_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    # 1b) group boxplot (5 groups, supplementary)
    fig, ax = plt.subplots(figsize=(12.5, 6.2))
    datp = bin_df[bin_df['PlotGroup'].isin(GROUP_ORDER)].copy()
    datp['PlotGroup_label'] = datp['PlotGroup'].map(lambda x: GROUP_LABEL.get(x, x))
    order = [g for g in GROUP_ORDER if (datp['PlotGroup'] == g).any()]
    order_label = [GROUP_LABEL[g] for g in order]
    palette = {GROUP_LABEL[g]: GROUP_COL[g] for g in order}
    sns.boxplot(data=datp, x='PlotGroup_label', y='anti_phage_load_bin', order=order_label, palette=palette, ax=ax, fliersize=2)
    sns.stripplot(data=datp, x='PlotGroup_label', y='anti_phage_load_bin', order=order_label, ax=ax, color='#111111', size=2, alpha=0.45, jitter=0.22)
    ax.set_title('Anti-phage load per MAG by context (compare-script groups)')
    ax.set_xlabel('')
    ax.set_ylabel('Anti-phage load per MAG (strict PADLOC + CRISPR_active)')
    ax.tick_params(axis='x', rotation=20)
    plt.tight_layout()
    fig.savefig(FIG_GROUP_PNG, dpi=DPI, bbox_inches='tight')
    fig.savefig(FIG_GROUP_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    # 2) paired biodigester medians
    fig, ax = plt.subplots(figsize=(8, 5.5))
    if not paired.empty:
        x = [0, 1]
        for _, r in paired.iterrows():
            ax.plot(x, [r['median_in'], r['median_out']], color='#666666', alpha=0.8, linewidth=1)
            ax.scatter([0], [r['median_in']], color=GROUP_COL['Biodigestor Swine In'], s=35, zorder=3)
            ax.scatter([1], [r['median_out']], color=GROUP_COL['Biodigestor Swine Out'], s=35, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Affl (In)', 'Effl (Out)'])
        ax.set_ylabel('Sample median anti-phage load per MAG')
        ax.set_title(f'Biodigester paired contrast (n_pairs={len(paired)})')
    else:
        ax.text(0.5, 0.5, 'No IN/OUT pairs found', ha='center', va='center')
        ax.axis('off')
    plt.tight_layout()
    fig.savefig(FIG_PAIR_PNG, dpi=DPI, bbox_inches='tight')
    fig.savefig(FIG_PAIR_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    # 3) intra-sample heterogeneity bar
    fig, ax = plt.subplots(figsize=(12.5, 5.5))
    plot_intra = intra.dropna(subset=['mean_pairwise_bray']).copy()
    plot_intra = plot_intra.sort_values('sample', key=lambda s: s.astype(int))
    colors = [GROUP_COL.get(g, '#999999') for g in plot_intra['PlotGroup']]
    ax.bar(plot_intra['sample'].astype(str), plot_intra['mean_pairwise_bray'].astype(float), color=colors, edgecolor='black', linewidth=0.2)
    for i, (_, r) in enumerate(plot_intra.iterrows()):
        if float(r['q_two_sided_bh']) < 0.05:
            ax.text(i, float(r['mean_pairwise_bray']) + 0.005, '*', ha='center', va='bottom', fontsize=11)
    present_groups = [g for g in GROUP_ORDER if (plot_intra['PlotGroup'] == g).any()]
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=GROUP_COL[g], ec='black', lw=0.2) for g in present_groups
    ]
    legend_labels = [GROUP_LABEL.get(g, g) for g in present_groups]
    ax.legend(
        legend_handles,
        legend_labels,
        title='Context groups',
        frameon=False,
        fontsize=8,
        title_fontsize=8,
        ncol=2,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.22),
    )
    ax.set_title('Intra-sample heterogeneity among MAGs (Bray-Curtis on strict anti-phage profile)')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Mean pairwise Bray-Curtis (within sample)')
    ax.tick_params(axis='x', rotation=70)
    plt.tight_layout()
    fig.savefig(FIG_INTRA_PNG, dpi=DPI, bbox_inches='tight')
    fig.savefig(FIG_INTRA_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    # report
    lines = []
    lines.append('Anti-phage context + intrasample MAG analysis')
    lines.append('')
    lines.append(f'Bins analyzed: {len(bin_df)}')
    lines.append(f'Bins with anti-phage signal (>0) used in profile tests: {len(Xrel_nz)}')
    lines.append(f'Samples analyzed: {bin_df["sample"].nunique()}')
    lines.append(f'Strict anti-phage systems kept: {len(strict_set)}')
    lines.append('')

    # summarize key p-values
    def pick_stat(analysis: str, test: str, col='p_value'):
        if context_stats.empty:
            return np.nan
        sub = context_stats[(context_stats['analysis'] == analysis) & (context_stats['test'] == test)]
        if sub.empty:
            return np.nan
        return float(sub.iloc[0][col])

    lines.append('Key tests (load-based):')
    p1 = pick_stat('anti_phage_load_by_plotgroup', 'Kruskal-Wallis')
    lines.append(f'- Across compare-script groups (Kruskal): p={p1}')
    p2 = pick_stat('broiler_vs_biodigestor_load', 'Mann-Whitney')
    lines.append(f'- Broiler litter vs Biodigestor (Mann-Whitney): p={p2}')
    p3 = pick_stat('biodigestor_in_vs_out_load', 'Mann-Whitney')
    lines.append(f'- Biodigestor Affl vs Effl (Mann-Whitney): p={p3}')
    p4 = pick_stat('biodigestor_paired_median_load', 'Wilcoxon paired')
    lines.append(f'- Biodigestor paired Affl/Effl by digester (Wilcoxon): p={p4}')

    if not perm_stats.empty:
        lines.append('')
        lines.append('Key tests (composition-based):')
        for _, r in perm_stats.iterrows():
            lines.append(f"- {r['analysis']} {r['test']}: F={r['F']}, p={r['p_value']}, R2={r['R2']}")

    sig_intra = intra[(pd.to_numeric(intra['q_two_sided_bh'], errors='coerce') < 0.05) & intra['mean_pairwise_bray'].notna()]
    lines.append('')
    lines.append(f'Intra-sample heterogeneity significant after BH (q<0.05): {len(sig_intra)} samples')
    if len(sig_intra) > 0:
        lines.append('- Samples: ' + ', '.join(sig_intra['sample'].astype(str).tolist()))

    OUT_REPORT.write_text('\n'.join(lines) + '\n', encoding='utf-8')

    print('[OK]', OUT_BIN)
    print('[OK]', OUT_INTRA)
    print('[OK]', OUT_CONTEXT)
    print('[OK]', OUT_PERM)
    print('[OK]', OUT_PAIR)
    print('[OK]', OUT_REPORT)
    print('[OK]', FIG_CONTEXT2_PNG)
    print('[OK]', FIG_GROUP_PNG)
    print('[OK]', FIG_PAIR_PNG)
    print('[OK]', FIG_INTRA_PNG)


if __name__ == '__main__':
    main()
