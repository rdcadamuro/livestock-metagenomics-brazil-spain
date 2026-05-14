#!/bin/bash
# Submete os 729 MAGs via webin-cli
# ANTES de correr: actualizar campo SAMPLE nos manifests com SAMEA accessions
# Uso: bash submit_all_assemblies.sh 'senha_webin'

PASSWORD=${1:?'Fornecer senha: bash submit_all_assemblies.sh SENHA'}
WEBIN_JAR=/path/to/webin-cli.jar
USERNAME=Webin-XXXXX  # substituir pelo username ENA

MANIFEST_DIR=/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags/ENA_final_submission_729/submission_files/manifests

ok=0; fail=0
for manifest in "$MANIFEST_DIR"/*_manifest.txt; do
    mag=$(basename "$manifest" _manifest.txt)
    echo "Submetendo $mag..."
    java -jar "$WEBIN_JAR" \
        -context genome \
        -manifest "$manifest" \
        -username "$USERNAME" \
        -password "$PASSWORD" \
        -submit \
        2>&1 | tee -a submit_log_${mag}.txt
    if [ $? -eq 0 ]; then ((ok++)); else ((fail++)); fi
done
echo "Concluido: $ok OK, $fail erros"
