#!/usr/bin/env python3
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


def clean(s: str) -> str:
    return (s or "").strip()


def norm_tax(s: str) -> str:
    x = clean(s)
    if not x:
        return ""
    low = x.lower()
    if low in {"-", "na", "nan", "none", "unknown", "unclassified", "no_inphared_mash_hit", "no predication"}:
        return ""
    if "no predication" in low:
        return ""
    if low.startswith("unclassified"):
        return ""
    return x


def contig_to_bin(contig_id: str) -> str:
    m = re.match(r"^vRhyme_(\d+)__", contig_id)
    if m:
        return f"vRhyme_bin_{m.group(1)}"
    return "unknown_bin"


def parse_lineage_deepest(lineage: str) -> str:
    if not lineage:
        return ""
    pairs = {}
    for part in lineage.split(";"):
        if ":" in part:
            k, v = part.split(":", 1)
            pairs[k.strip().lower()] = v.strip()
    for rank in ["genus", "family", "order", "class", "phylum", "kingdom", "realm", "clade", "superkingdom"]:
        v = norm_tax(pairs.get(rank, ""))
        if v:
            return v
    return ""


def majority_vote(values):
    vals = [v for v in values if norm_tax(v)]
    if not vals:
        return ""
    c = Counter(vals)
    return c.most_common(1)[0][0]


def build_consensus(tool_vals: dict):
    vals = [v for v in tool_vals.values() if norm_tax(v)]
    if not vals:
        return ("", "no_call", 0, 0, False)
    c = Counter(vals)
    top, n = c.most_common(1)[0]
    n_tools = len(vals)
    unique_n = len(c)
    all_agree = unique_n == 1 and n_tools >= 2

    if n >= 2:
        return (top, "vote_ge2", n_tools, unique_n, all_agree)

    for key in ["pharokka", "phabox", "phagenus", "blastn"]:
        v = norm_tax(tool_vals.get(key, ""))
        if v:
            return (v, f"single_{key}", n_tools, unique_n, all_agree)

    return ("", "no_call", n_tools, unique_n, all_agree)


def read_merged_contigs(merged_fasta: Path):
    contigs = []
    lengths = {}
    cur = None
    cur_len = 0
    if not merged_fasta.exists():
        return contigs, lengths
    with merged_fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = cur_len
                cur = line[1:].strip().split()[0]
                contigs.append(cur)
                cur_len = 0
            else:
                cur_len += len(line.strip())
    if cur is not None:
        lengths[cur] = cur_len
    return contigs, lengths


def read_phabox(phabox_tsv: Path):
    out = {}
    if not phabox_tsv.exists():
        return out
    with phabox_tsv.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            acc = clean(row.get("Accession", ""))
            if not acc:
                continue
            genus = norm_tax(row.get("Genus", ""))
            lineage = parse_lineage_deepest(row.get("Lineage", ""))
            out[acc] = genus or lineage
    return out


def read_phagenus(phagenus_csv: Path):
    out = {}
    if not phagenus_csv.exists():
        return out
    with phagenus_csv.open() as fh:
        first = fh.readline()
        fh.seek(0)
        # Some runs write a tab-delimited file with .csv extension.
        delim = "\t" if first.count("\t") > first.count(",") else ","
        r = csv.DictReader(fh, delimiter=delim)
        for row in r:
            cid = clean(row.get("contigs", "")) or clean(row.get("contig", "")) or clean(row.get("seq_id", ""))
            if not cid:
                continue
            pred = norm_tax(row.get("predict_label", "")) or norm_tax(row.get("predict", "")) or norm_tax(row.get("label", ""))
            out[cid] = pred
    return out


def read_blast(blast_tsv: Path, args):
    hits = defaultdict(list)
    if not blast_tsv.exists():
        return {}, {}

    with blast_tsv.open() as fh:
        r = csv.reader(fh, delimiter="\t")
        for row in r:
            if len(row) < 9:
                continue
            qseqid, sseqid, pident, length, evalue, bitscore, staxids, ssciname, qlen = row[:9]
            try:
                rec = {
                    "qseqid": qseqid,
                    "ssciname": clean(ssciname),
                    "pident": float(pident),
                    "length": float(length),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                    "qlen": float(qlen) if float(qlen) > 0 else 1.0,
                }
            except Exception:
                continue
            hits[qseqid].append(rec)

    calls = {}
    stats = {}
    for q, lst in hits.items():
        lst.sort(key=lambda x: (-x["bitscore"], x["evalue"]))
        filtered = []
        for h in lst:
            if h["evalue"] > args.blast_evalue_max:
                continue
            if h["pident"] < args.blast_pident_min:
                continue
            if h["length"] < args.blast_min_aln_bp:
                continue
            if (h["length"] / h["qlen"]) < args.blast_min_aln_frac_of_contig:
                continue
            filtered.append(h)

        top = filtered[:args.blast_topn]
        valid_n = len(top)
        call = ""
        maj = 0.0
        best_pident = 0.0
        best_qcov = 0.0
        ani_proxy_hits = 0
        ani_proxy_maj = 0.0
        ani_proxy_pass = False
        ani_proxy_taxon = ""

        if top:
            best_pident = max(h["pident"] for h in top)
            best_qcov = max((h["length"] / h["qlen"]) for h in top)

        if valid_n >= args.blast_min_valid_hits:
            c = Counter([h["ssciname"] for h in top if h["ssciname"]])
            if c:
                name, n = c.most_common(1)[0]
                maj = n / valid_n
                if maj >= args.blast_majority_frac:
                    call = name

        ani_top = []
        for h in top:
            qcov = h["length"] / h["qlen"]
            if h["pident"] >= args.blast_ani_proxy_pident_min and qcov >= args.blast_ani_proxy_qcov_min:
                ani_top.append(h)
        ani_proxy_hits = len(ani_top)
        if ani_proxy_hits > 0:
            ani_proxy_pass = True
            c_ani = Counter([h["ssciname"] for h in ani_top if h["ssciname"]])
            if c_ani:
                name, n = c_ani.most_common(1)[0]
                ani_proxy_maj = n / ani_proxy_hits
                ani_proxy_taxon = name

        calls[q] = call
        stats[q] = {
            "blast_valid_hits": valid_n,
            "blast_majority_frac": round(maj, 4),
            "blast_best_pident": round(best_pident, 4),
            "blast_best_qcov": round(best_qcov, 4),
            "blast_ani_proxy_hits": ani_proxy_hits,
            "blast_ani_proxy_majority_frac": round(ani_proxy_maj, 4),
            "blast_ani_proxy_pass": str(ani_proxy_pass),
            "blast_ani_proxy_taxon": ani_proxy_taxon,
        }

    return calls, stats


def read_pharokka(annotation_root: Path):
    out = {}
    md = {}

    files = sorted(annotation_root.glob("vRhyme_bin_*/vRhyme_bin_*_top_hits_mash_inphared.tsv"))
    for f in files:
        with f.open() as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                cid = clean(row.get("contig", ""))
                if not cid:
                    continue
                acc = norm_tax(row.get("Accession", ""))
                if not acc:
                    continue

                tax = ""
                for key in ["Lowest_Taxa", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm", "Description"]:
                    v = norm_tax(row.get(key, ""))
                    if v:
                        tax = v
                        break

                if tax:
                    out[cid] = tax

                try:
                    md[cid] = float(row.get("mash_distance", ""))
                except Exception:
                    pass

    return out, md


def integrate(args):
    annotation_root = Path(args.annotation_root)
    sample_level = annotation_root / "_sample_level"
    sample_tag = args.sample_tag or f"P19109_{args.sample}"

    merged = sample_level / f"{sample_tag}.vMAGs.clean.merged.fasta"
    phabox_tsv = sample_level / "phabox" / "final_prediction" / "phagcn_prediction.tsv"
    phagenus_csv = sample_level / "phagenus" / "prediction_output.csv"
    blast_tsv = sample_level / f"{sample_tag}.vmags.blastn.tsv"

    contigs, contig_lens = read_merged_contigs(merged)
    phabox = read_phabox(phabox_tsv)
    phagenus = read_phagenus(phagenus_csv)
    blast_calls, blast_stats = read_blast(blast_tsv, args)
    pharokka, pharokka_md = read_pharokka(annotation_root)

    all_contigs = set(contigs) | set(phabox) | set(phagenus) | set(blast_calls) | set(pharokka)

    contig_rows = []
    for cid in sorted(all_contigs):
        bid = contig_to_bin(cid)
        tool_vals = {
            "pharokka": norm_tax(pharokka.get(cid, "")),
            "phagenus": norm_tax(phagenus.get(cid, "")),
            "phabox": norm_tax(phabox.get(cid, "")),
            "blastn": norm_tax(blast_calls.get(cid, "")),
        }
        consensus, rule, n_tools, n_unique, all_agree = build_consensus(tool_vals)
        bs = blast_stats.get(
            cid,
            {
                "blast_valid_hits": 0,
                "blast_majority_frac": 0.0,
                "blast_best_pident": 0.0,
                "blast_best_qcov": 0.0,
                "blast_ani_proxy_hits": 0,
                "blast_ani_proxy_majority_frac": 0.0,
                "blast_ani_proxy_pass": "False",
                "blast_ani_proxy_taxon": "",
            },
        )

        contig_rows.append({
            "sample": args.sample,
            "bin": bid,
            "contig_id": cid,
            "pharokka_taxon": tool_vals["pharokka"],
            "pharokka_mash_distance": pharokka_md.get(cid, ""),
            "phagenus_taxon": tool_vals["phagenus"],
            "phabox_taxon": tool_vals["phabox"],
            "blastn_taxon": tool_vals["blastn"],
            "blast_valid_hits": bs["blast_valid_hits"],
            "blast_majority_frac": bs["blast_majority_frac"],
            "blast_best_pident": bs["blast_best_pident"],
            "blast_best_qcov": bs["blast_best_qcov"],
            "blast_ani_proxy_hits": bs["blast_ani_proxy_hits"],
            "blast_ani_proxy_majority_frac": bs["blast_ani_proxy_majority_frac"],
            "blast_ani_proxy_pass": bs["blast_ani_proxy_pass"],
            "blast_ani_proxy_taxon": bs["blast_ani_proxy_taxon"],
            "n_tools_called": n_tools,
            "n_unique_calls": n_unique,
            "all_tools_agree": str(all_agree),
            "consensus_taxon": consensus,
            "consensus_rule": rule,
        })

    by_bin = defaultdict(list)
    for r in contig_rows:
        by_bin[r["bin"]].append(r)

    bin_rows = []
    for bid in sorted(by_bin):
        rows = by_bin[bid]
        bvals = {
            "pharokka": majority_vote([x["pharokka_taxon"] for x in rows]),
            "phagenus": majority_vote([x["phagenus_taxon"] for x in rows]),
            "phabox": majority_vote([x["phabox_taxon"] for x in rows]),
            "blastn": majority_vote([x["blastn_taxon"] for x in rows]),
        }
        consensus, rule, n_tools, n_unique, all_agree = build_consensus(bvals)
        ani_rows = [x for x in rows if x.get("blast_ani_proxy_pass", "False") == "True"]
        n_ani = len(ani_rows)
        frac_ani = (n_ani / len(rows)) if rows else 0.0
        ani_taxa = set(norm_tax(x.get("blast_ani_proxy_taxon", "")) for x in ani_rows if norm_tax(x.get("blast_ani_proxy_taxon", "")))
        ani_proxy_tax_conflict = len(ani_taxa) > 1
        representative_ani_proxy_pass = any(
            (x.get("blast_ani_proxy_pass", "False") == "True")
            and (int(contig_lens.get(x.get("contig_id", ""), 0)) >= args.rep_contig_min_bp)
            for x in rows
        )

        bin_rows.append({
            "sample": args.sample,
            "bin": bid,
            "n_contigs": len(rows),
            "pharokka_taxon_bin": bvals["pharokka"],
            "phagenus_taxon_bin": bvals["phagenus"],
            "phabox_taxon_bin": bvals["phabox"],
            "blastn_taxon_bin": bvals["blastn"],
            "n_tools_called": n_tools,
            "n_unique_calls": n_unique,
            "all_tools_agree": str(all_agree),
            "consensus_taxon_bin": consensus,
            "consensus_rule": rule,
            "n_contigs_with_any_call": sum(1 for x in rows if int(x["n_tools_called"]) > 0),
            "n_contigs_ani_proxy_pass": n_ani,
            "frac_contigs_ani_proxy_pass": round(frac_ani, 4),
            "representative_ani_proxy_pass": str(representative_ani_proxy_pass),
            "ani_proxy_tax_conflict": str(ani_proxy_tax_conflict),
        })

    contig_out = sample_level / f"taxonomy_4tools_contig_{args.sample}.tsv"
    bin_out = sample_level / f"taxonomy_4tools_bin_{args.sample}.tsv"
    summary_out = sample_level / f"taxonomy_4tools_summary_{args.sample}.tsv"

    contig_fields = [
        "sample", "bin", "contig_id", "pharokka_taxon", "pharokka_mash_distance", "phagenus_taxon",
        "phabox_taxon", "blastn_taxon", "blast_valid_hits", "blast_majority_frac", "blast_best_pident",
        "blast_best_qcov", "blast_ani_proxy_hits", "blast_ani_proxy_majority_frac", "blast_ani_proxy_pass",
        "blast_ani_proxy_taxon", "n_tools_called",
        "n_unique_calls", "all_tools_agree", "consensus_taxon", "consensus_rule"
    ]
    bin_fields = [
        "sample", "bin", "n_contigs", "pharokka_taxon_bin", "phagenus_taxon_bin", "phabox_taxon_bin",
        "blastn_taxon_bin", "n_tools_called", "n_unique_calls", "all_tools_agree", "consensus_taxon_bin",
        "consensus_rule", "n_contigs_with_any_call", "n_contigs_ani_proxy_pass", "frac_contigs_ani_proxy_pass",
        "representative_ani_proxy_pass", "ani_proxy_tax_conflict"
    ]

    sample_level.mkdir(parents=True, exist_ok=True)

    with contig_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=contig_fields, delimiter="\t")
        w.writeheader()
        w.writerows(contig_rows)

    with bin_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=bin_fields, delimiter="\t")
        w.writeheader()
        w.writerows(bin_rows)

    n_bins = len(bin_rows)
    n_any = sum(1 for x in bin_rows if int(x["n_tools_called"]) >= 1)
    n_ge2 = sum(1 for x in bin_rows if int(x["n_tools_called"]) >= 2)
    n_agree = sum(1 for x in bin_rows if x["all_tools_agree"] == "True")

    with summary_out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["sample", args.sample])
        w.writerow(["n_bins", n_bins])
        w.writerow(["n_bins_any_call", n_any])
        w.writerow(["n_bins_tools_ge2", n_ge2])
        w.writerow(["n_bins_all_tools_agree", n_agree])
        w.writerow(["blast_topn", args.blast_topn])
        w.writerow(["blast_max_target_seqs", args.blast_max_target_seqs])
        w.writerow(["blast_evalue_max", args.blast_evalue_max])
        w.writerow(["blast_pident_min", args.blast_pident_min])
        w.writerow(["blast_min_aln_bp", args.blast_min_aln_bp])
        w.writerow(["blast_min_aln_frac_of_contig", args.blast_min_aln_frac_of_contig])
        w.writerow(["blast_min_valid_hits", args.blast_min_valid_hits])
        w.writerow(["blast_majority_frac", args.blast_majority_frac])
        w.writerow(["blast_ani_proxy_pident_min", args.blast_ani_proxy_pident_min])
        w.writerow(["blast_ani_proxy_qcov_min", args.blast_ani_proxy_qcov_min])
        w.writerow(["rep_contig_min_bp", args.rep_contig_min_bp])

    print(contig_out)
    print(bin_out)
    print(summary_out)


def main():
    p = argparse.ArgumentParser(description="Build 4-tool taxonomy consensus (Pharokka, PhaGenus, PhaBOX, BLASTn).")
    p.add_argument("--sample", required=True, help="Sample numeric id, e.g. 1035")
    p.add_argument("--annotation-root", required=True, help="Sample annotation root dir containing _sample_level and vRhyme_bin_* dirs")
    p.add_argument("--sample-tag", default="", help="Sample tag, default P19109_<sample>")

    p.add_argument("--blast-topn", type=int, default=10)
    p.add_argument("--blast-max-target-seqs", type=int, default=50)
    p.add_argument("--blast-evalue-max", type=float, default=1e-20)
    p.add_argument("--blast-pident-min", type=float, default=85.0)
    p.add_argument("--blast-min-aln-bp", type=float, default=500)
    p.add_argument("--blast-min-aln-frac-of-contig", type=float, default=0.05)
    p.add_argument("--blast-min-valid-hits", type=int, default=3)
    p.add_argument("--blast-majority-frac", type=float, default=0.60)
    p.add_argument("--blast-ani-proxy-pident-min", type=float, default=95.0)
    p.add_argument("--blast-ani-proxy-qcov-min", type=float, default=0.85)
    p.add_argument("--rep-contig-min-bp", type=int, default=10000)

    args = p.parse_args()
    integrate(args)


if __name__ == "__main__":
    main()
