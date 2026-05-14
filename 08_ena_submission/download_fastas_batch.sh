#!/bin/bash
# Download all 729 MAG FASTAs from Google Drive (rclone remote: ubu:)
# and organize into fastas/Brazil/ and fastas/Spain/ with standard naming.
#
# Brazil/<mag_id>.fa  e  Spain/<mag_id>.fa
#
# Run from: /mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags/

set -e

OUTDIR="ENA_final_submission_729/fastas"
BR="$OUTDIR/Brazil"
SP="$OUTDIR/Spain"
TMP="$OUTDIR/_tmp"

mkdir -p "$BR" "$SP" "$TMP/brazil_rescued" "$TMP/spain_rescued" "$TMP/mags_bins2" "$TMP/spain_mags"

echo "========================================"
echo "PASSO 1/4 — Brazil rescued (113 MAGs)"
echo "========================================"
rclone copy "ubu:home/Tati_Metagenome_Mags_Filtered/Brazil/rescued" "$TMP/brazil_rescued/" \
    --progress --transfers=8
echo "  Renomeando..."
for f in "$TMP/brazil_rescued/"*.final.fa; do
    [ -f "$f" ] || continue
    base=$(basename "$f" .final.fa)
    cp "$f" "$BR/${base}.fa"
done
echo "  Concluido: $(ls $BR/*.fa 2>/dev/null | wc -l) em Brazil/"

echo ""
echo "========================================"
echo "PASSO 2/4 — Spain rescued (59 MAGs)"
echo "========================================"
rclone copy "ubu:home/Tati_Metagenome_Mags_Filtered/Spain/rescued" "$TMP/spain_rescued/" \
    --progress --transfers=8
echo "  Renomeando..."
for f in "$TMP/spain_rescued/"*.final.fa; do
    [ -f "$f" ] || continue
    base=$(basename "$f" .final.fa)
    cp "$f" "$SP/${base}.fa"
done
echo "  Concluido: $(ls $SP/*.fa 2>/dev/null | wc -l) em Spain/"

echo ""
echo "========================================"
echo "PASSO 3/4 — Brazil MAGS_BINS_2 metawrap (444 MAGs)"
echo "========================================"
rclone copy "ubu:home/MAGS_BINS_2" "$TMP/mags_bins2/" \
    --include "*/metawrap_50_10_bins/*.fa" \
    --progress --transfers=8
echo "  Renomeando..."
for sample_dir in "$TMP/mags_bins2/P19109_"*; do
    [ -d "$sample_dir" ] || continue
    sample=$(basename "$sample_dir")
    for f in "$sample_dir/metawrap_50_10_bins/"bin.*.fa; do
        [ -f "$f" ] || continue
        bnum=$(basename "$f" .fa | sed 's/bin\.//')
        mag_id="Brazil__${sample}__bin.${bnum}"
        # Only copy if this mag_id is in our submission list (skip excluded bins)
        dest="$BR/${mag_id}.fa"
        if [ ! -f "$dest" ]; then
            cp "$f" "$dest"
        fi
    done
done
echo "  Concluido: $(ls $BR/*.fa 2>/dev/null | wc -l) em Brazil/"

echo ""
echo "========================================"
echo "PASSO 4/4 — Spain Metagenome_Spain/MAGs dastool (113 MAGs)"
echo "========================================"
rclone copy "ubu:home/Metagenome_Spain/MAGs" "$TMP/spain_mags/" \
    --include "*/dastool/*_DASTOOL_bins/*.fa" \
    --progress --transfers=8
echo "  Renomeando..."
for srr_dir in "$TMP/spain_mags/SRR"*; do
    [ -d "$srr_dir" ] || continue
    srr=$(basename "$srr_dir")
    dastool_dir="$srr_dir/dastool/${srr}_DASTOOL_bins"
    [ -d "$dastool_dir" ] || dastool_dir=$(ls -d "$srr_dir"/dastool/*_DASTOOL_bins 2>/dev/null | head -1)
    [ -d "$dastool_dir" ] || continue
    for f in "$dastool_dir/"bin_*.fa; do
        [ -f "$f" ] || continue
        bname=$(basename "$f" .fa | sed 's/^bin_//')
        mag_id="Spain__${srr}__bin_${bname}"
        dest="$SP/${mag_id}.fa"
        if [ ! -f "$dest" ]; then
            cp "$f" "$dest"
        fi
    done
done
echo "  Concluido: $(ls $SP/*.fa 2>/dev/null | wc -l) em Spain/"

echo ""
echo "========================================"
echo "VERIFICAÇÃO FINAL"
echo "========================================"
brazil_count=$(ls "$BR/"*.fa 2>/dev/null | wc -l)
spain_count=$(ls "$SP/"*.fa 2>/dev/null | wc -l)
total=$((brazil_count + spain_count))
echo "  Brazil: $brazil_count FASTAs"
echo "  Spain:  $spain_count FASTAs"
echo "  Total:  $total / 729 esperados"

if [ "$total" -eq 729 ]; then
    echo "  SUCESSO — todos os 729 MAGs baixados!"
    echo "  Limpando temporários..."
    rm -rf "$TMP"
else
    echo "  ATENCAO — $((729 - total)) MAGs ainda faltando. Ver $TMP para debug."
fi
