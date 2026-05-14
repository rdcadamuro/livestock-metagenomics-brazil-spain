#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


CHECKV_RANK = {
    "Not-determined": 0,
    "Low-quality": 1,
    "Medium-quality": 2,
    "High-quality": 3,
    "Complete": 4,
}

CONF_RANK = {"Low": 0, "Medium": 1, "High": 2}


def to_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def to_int(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default


def to_bool(x):
    return str(x).strip().lower() in {"true", "1", "yes", "y"}


def parse_cov_from_contig(contig_id: str):
    m = re.search(r"_cov_([0-9]+(?:\.[0-9]+)?)", contig_id)
    if m:
        return to_float(m.group(1), default=math.nan)
    return math.nan


def contig_to_bin(contig_id: str):
    m = re.match(r"^vRhyme_(\d+)__", contig_id)
    if m:
        return f"vRhyme_bin_{m.group(1)}"
    return "unknown_bin"


def fasta_stats(fasta_path: Path):
    ids = []
    lengths = {}
    cur = None
    cur_len = 0

    if not fasta_path.exists():
        return ids, lengths

    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = cur_len
                cur = line[1:].strip().split()[0]
                ids.append(cur)
                cur_len = 0
            else:
                cur_len += len(line.strip())
    if cur is not None:
        lengths[cur] = cur_len

    return ids, lengths


def read_checkv_bin(qfile: Path):
    if not qfile.exists():
        return {
            "n_contigs": 0,
            "viral_positive": 0,
            "viral_genes_sum": 0,
            "host_genes_sum": 0,
            "max_completeness": 0.0,
            "best_quality": "Not-determined",
            "no_viral_warning": 0,
        }

    n_contigs = 0
    viral_positive = 0
    viral_genes_sum = 0
    host_genes_sum = 0
    max_completeness = 0.0
    best_quality = "Not-determined"
    no_viral_warning = 0

    with qfile.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_contigs += 1
            vg = to_int(row.get("viral_genes", 0))
            hg = to_int(row.get("host_genes", 0))
            comp = to_float(row.get("completeness", 0.0))
            qual = row.get("checkv_quality", "Not-determined")
            warn = (row.get("warnings", "") or "").lower()

            if vg > 0:
                viral_positive += 1
            viral_genes_sum += vg
            host_genes_sum += hg
            if comp > max_completeness:
                max_completeness = comp
            if CHECKV_RANK.get(qual, 0) > CHECKV_RANK.get(best_quality, 0):
                best_quality = qual
            if "no viral genes detected" in warn:
                no_viral_warning += 1

    return {
        "n_contigs": n_contigs,
        "viral_positive": viral_positive,
        "viral_genes_sum": viral_genes_sum,
        "host_genes_sum": host_genes_sum,
        "max_completeness": round(max_completeness, 2),
        "best_quality": best_quality,
        "no_viral_warning": no_viral_warning,
    }


def read_taxonomy_bin(bin_tsv: Path):
    out = {}
    if not bin_tsv.exists():
        return out

    with bin_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            bid = row.get("bin", "")
            if not bid:
                continue
            out[bid] = {
                "n_tools_called": to_int(row.get("n_tools_called", 0)),
                "all_tools_agree": str(row.get("all_tools_agree", "False")) == "True",
                "n_unique_calls": to_int(row.get("n_unique_calls", 0)),
                "consensus_rule": row.get("consensus_rule", ""),
                "consensus_taxon_bin": row.get("consensus_taxon_bin", ""),
                "n_contigs_ani_proxy_pass": to_int(row.get("n_contigs_ani_proxy_pass", 0)),
                "frac_contigs_ani_proxy_pass": to_float(row.get("frac_contigs_ani_proxy_pass", 0.0), default=0.0),
                "representative_ani_proxy_pass": to_bool(row.get("representative_ani_proxy_pass", "False")),
                "ani_proxy_tax_conflict": to_bool(row.get("ani_proxy_tax_conflict", "False")),
            }
    return out


def read_pharokka_mash(annotation_root: Path):
    dist = defaultdict(list)
    for f in sorted(annotation_root.glob("vRhyme_bin_*/vRhyme_bin_*_top_hits_mash_inphared.tsv")):
        bid = f.parent.name
        with f.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                d = to_float(row.get("mash_distance", "nan"), default=math.nan)
                if not math.isnan(d):
                    dist[bid].append(d)

    out = {}
    for bid, vals in dist.items():
        vals = sorted(vals)
        out[bid] = vals[0] if vals else math.nan
    return out


def read_pharokka_functional(annotation_root: Path):
    out = {}
    structural_kw = [
        "capsid", "portal", "terminase", "tail", "baseplate", "sheath", "head", "neck",
        "fiber", "fibre", "spike", "tape measure",
    ]
    replication_kw = ["helicase", "primase", "polymerase", "ligase", "ssb", "dna-binding", "replisome"]
    lysogeny_kw = ["integrase", "excisionase", "repressor", "cro", "xis", "prophage"]
    host_like_kw = ["gyrb", "rpob", "rpoc", "reca", "dnak", "atpa", "atpb", "rpl", "rps"]

    for gff in sorted(annotation_root.glob("vRhyme_bin_*/vRhyme_bin_*.gff")):
        bid = gff.parent.name
        term = capsid = portal = tail = 0
        structural_hits = replication_hits = lysogeny_hits = host_like_hits = 0

        with gff.open() as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                feature = cols[2].lower()
                if feature not in {"cds", "gene"}:
                    continue
                low = cols[8].lower()

                if any(k in low for k in structural_kw):
                    structural_hits += 1
                if any(k in low for k in replication_kw):
                    replication_hits += 1
                if any(k in low for k in lysogeny_kw):
                    lysogeny_hits += 1
                if any(k in low for k in host_like_kw):
                    host_like_hits += 1

                if "terminase" in low:
                    term += 1
                if "capsid" in low:
                    capsid += 1
                if "portal" in low:
                    portal += 1
                if "tail" in low or "baseplate" in low or "sheath" in low or "fiber" in low or "fibre" in low:
                    tail += 1

        structural_core_present = (term >= 1) and ((capsid >= 1) or (portal >= 1)) and (tail >= 1)
        dual_structural_red_flag = (term >= 2) and (capsid >= 2) and (portal >= 2)
        host_like_red_flag = host_like_hits >= 3

        if dual_structural_red_flag or host_like_red_flag:
            status = "red_flag"
        elif structural_core_present:
            status = "coherent"
        elif structural_hits >= 2:
            status = "partial"
        elif (structural_hits + replication_hits + lysogeny_hits) == 0:
            status = "no_signal"
        else:
            status = "weak"

        out[bid] = {
            "functional_status": status,
            "structural_hits": structural_hits,
            "replication_hits": replication_hits,
            "lysogeny_hits": lysogeny_hits,
            "host_like_hits": host_like_hits,
            "terminase_hits": term,
            "capsid_hits": capsid,
            "portal_hits": portal,
            "tail_hits": tail,
            "structural_core_present": structural_core_present,
            "dual_structural_red_flag": dual_structural_red_flag,
            "host_like_red_flag": host_like_red_flag,
        }

    return out


def score_viral_purity(v):
    if v >= 0.90:
        return 35
    if v >= 0.75:
        return 25
    if v >= 0.60:
        return 15
    return 5


def score_cov_cv(cv):
    if math.isnan(cv):
        return 10
    if cv <= 0.25:
        return 20
    if cv <= 0.45:
        return 14
    if cv <= 0.70:
        return 8
    return 2


def score_taxonomy(n_tools, all_agree):
    if n_tools >= 2 and all_agree:
        return 20
    if n_tools >= 2:
        return 14
    if n_tools == 1:
        return 8
    return 0


def score_completeness(max_comp):
    if max_comp >= 50:
        return 15
    if max_comp >= 30:
        return 10
    if max_comp >= 10:
        return 6
    return 2


def score_host_ratio(host_ratio):
    if host_ratio <= 0.10:
        return 10
    if host_ratio <= 0.30:
        return 7
    if host_ratio <= 0.60:
        return 4
    return 1


def score_mash(best_mash):
    if math.isnan(best_mash):
        return 5
    if best_mash <= 0.10:
        return 10
    if best_mash <= 0.20:
        return 7
    if best_mash <= 0.35:
        return 4
    return 2


def final_confidence(total):
    if total >= 70:
        return "High"
    if total >= 45:
        return "Medium"
    return "Low"


def write_selected_fasta(keep_bins, bins_dir: Path, out_fasta: Path):
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w") as out:
        for bid in sorted(keep_bins):
            fa = bins_dir / f"{bid}.clean.fasta"
            if not fa.exists():
                continue
            with fa.open() as fh:
                out.write(fh.read())


def integrate(args):
    annotation_root = Path(args.annotation_root)
    sample_level = annotation_root / "_sample_level"
    checkv_root = Path(args.checkv_root)
    bins_dir = Path(args.bins_dir)
    sample_tag = args.sample_tag or f"P19109_{args.sample}"

    tax_bin_tsv = sample_level / f"taxonomy_4tools_bin_{args.sample}.tsv"
    tax_bin = read_taxonomy_bin(tax_bin_tsv)
    mash_by_bin = read_pharokka_mash(annotation_root)
    func_by_bin = read_pharokka_functional(annotation_root)

    bin_ids = set()
    bin_ids.update([p.name.replace(".clean.fasta", "") for p in bins_dir.glob("vRhyme_bin_*.clean.fasta")])
    bin_ids.update([p.name for p in checkv_root.glob("vRhyme_bin_*") if p.is_dir()])
    bin_ids.update(tax_bin.keys())

    rows = []
    for bid in sorted(bin_ids):
        clean_fa = bins_dir / f"{bid}.clean.fasta"
        qfile = checkv_root / bid / "quality_summary.tsv"

        ids, lengths = fasta_stats(clean_fa)
        clean_n = len(ids)
        clean_bp = sum(lengths.values())
        covs = [parse_cov_from_contig(cid) for cid in ids]
        covs = [x for x in covs if not math.isnan(x)]
        if covs:
            mean_cov = sum(covs) / len(covs)
            sd_cov = math.sqrt(sum((x - mean_cov) ** 2 for x in covs) / len(covs))
            cov_cv = sd_cov / mean_cov if mean_cov > 0 else math.nan
        else:
            cov_cv = math.nan

        q = read_checkv_bin(qfile)
        tax = tax_bin.get(bid, {
            "n_tools_called": 0,
            "all_tools_agree": False,
            "n_unique_calls": 0,
            "consensus_rule": "",
            "consensus_taxon_bin": "",
            "n_contigs_ani_proxy_pass": 0,
            "frac_contigs_ani_proxy_pass": 0.0,
            "representative_ani_proxy_pass": False,
            "ani_proxy_tax_conflict": False,
        })
        func = func_by_bin.get(
            bid,
            {
                "functional_status": "no_data",
                "structural_hits": 0,
                "replication_hits": 0,
                "lysogeny_hits": 0,
                "host_like_hits": 0,
                "terminase_hits": 0,
                "capsid_hits": 0,
                "portal_hits": 0,
                "tail_hits": 0,
                "structural_core_present": False,
                "dual_structural_red_flag": False,
                "host_like_red_flag": False,
            },
        )

        viral_purity = (q["viral_positive"] / q["n_contigs"]) if q["n_contigs"] > 0 else 0.0
        host_ratio = q["host_genes_sum"] / (q["viral_genes_sum"] + q["host_genes_sum"] + 1e-9)
        best_mash = mash_by_bin.get(bid, math.nan)

        s1 = score_viral_purity(viral_purity)
        s2 = score_cov_cv(cov_cv)
        s3 = score_taxonomy(tax["n_tools_called"], tax["all_tools_agree"])
        s4 = score_completeness(q["max_completeness"])
        s5 = score_host_ratio(host_ratio)
        s6 = score_mash(best_mash)
        total = s1 + s2 + s3 + s4 + s5 + s6
        conf = final_confidence(total)

        cov_ok = math.isnan(cov_cv) or cov_cv <= 0.45
        ani_ok = tax["representative_ani_proxy_pass"] or tax["frac_contigs_ani_proxy_pass"] >= 0.30
        tax_ok = (tax["n_tools_called"] >= 2 and not tax["ani_proxy_tax_conflict"]) or ani_ok
        red_flag = tax["ani_proxy_tax_conflict"] or func["dual_structural_red_flag"] or func["host_like_red_flag"]

        if red_flag:
            same_virus_call = "Fail_or_mixed"
        elif viral_purity >= 0.75 and host_ratio <= 0.30 and q["no_viral_warning"] == 0 and cov_ok and tax_ok:
            same_virus_call = "Pass_strong"
        elif viral_purity >= 0.60 and host_ratio <= 0.60 and (tax["n_tools_called"] >= 1 or ani_ok or (not math.isnan(best_mash) and best_mash <= 0.35)):
            same_virus_call = "Pass_moderate"
        else:
            same_virus_call = "Uncertain"

        reasons = []
        if tax["ani_proxy_tax_conflict"]:
            reasons.append("ani_proxy_tax_conflict")
        if func["dual_structural_red_flag"]:
            reasons.append("dual_structural_modules")
        if func["host_like_red_flag"]:
            reasons.append("host_like_gene_enrichment")
        if not cov_ok:
            reasons.append("coverage_dispersion_high")
        if q["no_viral_warning"] > 0:
            reasons.append("checkv_no_viral_warning")
        same_virus_notes = ";".join(reasons)

        rows.append({
            "sample": args.sample,
            "bin": bid,
            "clean_contigs": clean_n,
            "clean_bp": clean_bp,
            "checkv_contigs": q["n_contigs"],
            "checkv_best_quality": q["best_quality"],
            "checkv_max_completeness": q["max_completeness"],
            "checkv_no_viral_warning_n": q["no_viral_warning"],
            "viral_purity": round(viral_purity, 4),
            "host_gene_ratio": round(host_ratio, 4),
            "coverage_cv": "" if math.isnan(cov_cv) else round(cov_cv, 4),
            "taxonomy_tools_called": tax["n_tools_called"],
            "taxonomy_all_tools_agree": str(tax["all_tools_agree"]),
            "taxonomy_n_unique_calls": tax["n_unique_calls"],
            "taxonomy_consensus_rule": tax["consensus_rule"],
            "taxonomy_consensus_taxon": tax["consensus_taxon_bin"],
            "taxonomy_n_contigs_ani_proxy_pass": tax["n_contigs_ani_proxy_pass"],
            "taxonomy_frac_contigs_ani_proxy_pass": round(tax["frac_contigs_ani_proxy_pass"], 4),
            "taxonomy_representative_ani_proxy_pass": str(tax["representative_ani_proxy_pass"]),
            "taxonomy_ani_proxy_tax_conflict": str(tax["ani_proxy_tax_conflict"]),
            "pharokka_best_mash_distance": "" if math.isnan(best_mash) else round(best_mash, 4),
            "functional_status": func["functional_status"],
            "functional_structural_hits": func["structural_hits"],
            "functional_replication_hits": func["replication_hits"],
            "functional_lysogeny_hits": func["lysogeny_hits"],
            "functional_host_like_hits": func["host_like_hits"],
            "functional_terminase_hits": func["terminase_hits"],
            "functional_capsid_hits": func["capsid_hits"],
            "functional_portal_hits": func["portal_hits"],
            "functional_tail_hits": func["tail_hits"],
            "functional_structural_core_present": str(func["structural_core_present"]),
            "functional_dual_structural_red_flag": str(func["dual_structural_red_flag"]),
            "functional_host_like_red_flag": str(func["host_like_red_flag"]),
            "score_viral_purity": s1,
            "score_coverage_cv": s2,
            "score_taxonomy": s3,
            "score_completeness": s4,
            "score_host_ratio": s5,
            "score_mash": s6,
            "score_total": total,
            "final_confidence": conf,
            "same_virus_call": same_virus_call,
            "same_virus_notes": same_virus_notes,
        })

    keep_min_rank = CONF_RANK[args.keep_min_confidence]
    for r in rows:
        r["keep"] = str(CONF_RANK[r["final_confidence"]] >= keep_min_rank)

    sample_level.mkdir(parents=True, exist_ok=True)
    all_out = sample_level / f"vmag_confidence_{args.sample}.tsv"
    keep_out = sample_level / f"vmag_confidence_filtered_{args.sample}_{args.keep_min_confidence.lower()}plus.tsv"
    summary_out = sample_level / f"vmag_confidence_{args.sample}.summary.tsv"
    selected_merged = sample_level / f"{sample_tag}.vMAGs.selected_{args.keep_min_confidence.lower()}plus.merged.fasta"

    fields = [
        "sample", "bin", "clean_contigs", "clean_bp", "checkv_contigs", "checkv_best_quality",
        "checkv_max_completeness", "checkv_no_viral_warning_n", "viral_purity", "host_gene_ratio",
        "coverage_cv", "taxonomy_tools_called", "taxonomy_all_tools_agree", "taxonomy_n_unique_calls",
        "taxonomy_consensus_rule", "taxonomy_consensus_taxon", "taxonomy_n_contigs_ani_proxy_pass",
        "taxonomy_frac_contigs_ani_proxy_pass", "taxonomy_representative_ani_proxy_pass",
        "taxonomy_ani_proxy_tax_conflict", "pharokka_best_mash_distance", "functional_status",
        "functional_structural_hits", "functional_replication_hits", "functional_lysogeny_hits",
        "functional_host_like_hits", "functional_terminase_hits", "functional_capsid_hits",
        "functional_portal_hits", "functional_tail_hits", "functional_structural_core_present",
        "functional_dual_structural_red_flag", "functional_host_like_red_flag", "score_viral_purity", "score_coverage_cv",
        "score_taxonomy", "score_completeness", "score_host_ratio", "score_mash", "score_total",
        "final_confidence", "same_virus_call", "same_virus_notes", "keep"
    ]

    with all_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    kept = [r for r in rows if r["keep"] == "True"]
    with keep_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(kept)

    write_selected_fasta([r["bin"] for r in kept], bins_dir, selected_merged)

    conf_count = Counter([r["final_confidence"] for r in rows])
    with summary_out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["sample", args.sample])
        w.writerow(["keep_min_confidence", args.keep_min_confidence])
        w.writerow(["n_bins", len(rows)])
        w.writerow(["n_keep", len(kept)])
        w.writerow(["n_high", conf_count.get("High", 0)])
        w.writerow(["n_medium", conf_count.get("Medium", 0)])
        w.writerow(["n_low", conf_count.get("Low", 0)])
        w.writerow(["selected_merged_fasta", str(selected_merged)])

    print(all_out)
    print(keep_out)
    print(summary_out)
    print(selected_merged)


def main():
    p = argparse.ArgumentParser(description="Score final confidence of vMAG bins and filter by threshold.")
    p.add_argument("--sample", required=True, help="Sample numeric id, e.g. 1035")
    p.add_argument("--annotation-root", required=True, help="Sample annotation root dir (contains _sample_level and vRhyme_bin_* dirs)")
    p.add_argument("--checkv-root", required=True, help="CheckV bins root dir for sample")
    p.add_argument("--bins-dir", required=True, help="vRhyme bins fasta dir for sample")
    p.add_argument("--sample-tag", default="", help="Sample tag, default P19109_<sample>")
    p.add_argument("--keep-min-confidence", choices=["Low", "Medium", "High"], default="Medium")

    args = p.parse_args()
    integrate(args)


if __name__ == "__main__":
    main()
