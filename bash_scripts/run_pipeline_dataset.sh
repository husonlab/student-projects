#!/bin/bash
#SBATCH --job-name=pangenome_matrix
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=8000m
#SBATCH --output=pangenome_%j.log
#
# Wie run_pipeline_with_members.sh, aber Input-Datei und Output-Ordner sind
# Parameter statt fest verdrahtet -- so kann derselbe Job fuer mehrere
# Datensaetze gestartet werden, ohne dass sie sich gegenseitig ueberschreiben.
#
# Jede Stufe laeuft unter /usr/bin/time -v, damit am Ende eine Tabelle mit
# Dauer + Speicherverbrauch PRO STUFE ausgegeben wird, nicht nur fuer den
# ganzen Job insgesamt.
#
# Aufruf: sbatch run_pipeline_dataset.sh <input_clusters.tsv> <output_dir_name>
# Beispiel:
#   sbatch run_pipeline_dataset.sh datasets/50_id-80_cov.tsv output_50id80cov
#   sbatch run_pipeline_dataset.sh datasets/50_id-80_cov-covmode1.tsv output_50id80cov_covmode1
# (beide koennen gleichzeitig laufen -- unabhaengige SLURM-Jobs)

set -e

cd /pfs/10/work/tu_zxozj42-lm_pangenome/magrix

if [ $# -ne 2 ]; then
    echo "Usage: sbatch run_pipeline_dataset.sh <input_clusters.tsv> <output_dir_name>"
    exit 1
fi
INPUT_TSV="$1"
OUT="$2"
mkdir -p "$OUT"

echo "Datensatz : $INPUT_TSV"
echo "Output    : $OUT"
# Hardware und Code-Stand fuer den Methodik-Teil festhalten
echo "Knoten    : $(hostname)"
LC_ALL=C lscpu | grep -E "Model name|^Socket|^CPU\(s\)" | sed 's/^/            /'
echo "Commit    : $(git rev-parse --short HEAD 2>/dev/null || echo unbekannt)"

# ── Stufen-Timing/Speicher-Tracking ─────────────────────────────────────────
declare -a STAGE_ORDER=()
declare -A STAGE_TIME
declare -A STAGE_MEM

# Fuehrt ein Kommando unter /usr/bin/time -v aus, zeigt dessen normale
# Ausgabe weiterhin im Log an, und merkt sich Dauer + Max-RAM fuer die
# Zusammenfassung am Ende. $1 = Label fuer die Stufe, Rest = Kommando.
run_stage() {
    local label="$1"
    shift
    STAGE_ORDER+=("$label")
    local timelog
    timelog=$(mktemp)
    /usr/bin/time -v "$@" 2>"$timelog"
    cat "$timelog"   # Stufen-eigene Fortschrittsausgabe bleibt im SLURM-Log sichtbar
    STAGE_TIME["$label"]=$(grep "Elapsed (wall clock)" "$timelog" | sed 's/.*: //')
    STAGE_MEM["$label"]=$(grep "Maximum resident set size" "$timelog" | sed 's/.*: //')
    rm -f "$timelog"
}

echo "Stufe 1: main pipeline"
date
run_stage "1_main" ./target/release/pangenome_matrices "$INPUT_TSV" "$OUT"

NUM_FAM=$(($(wc -l < "$OUT/family_dictionary.tsv") - 1))
NUM_MAGS=$(($(wc -l < "$OUT/mag_dictionary.tsv") - 1))
echo "num_families=$NUM_FAM num_mags=$NUM_MAGS"

echo "Stufe 2: CSC/CSR fuer paralog"
date
run_stage "2a_csc" ./target/release/build_csc "$OUT/mmseqs_long_paralog.tsv" $NUM_FAM "$OUT/paralog"
run_stage "2b_csr" ./target/release/build_csr "$OUT/mmseqs_long_paralog.tsv" $NUM_MAGS "$OUT/paralog"
rm -f "$OUT/mmseqs_long_paralog.tsv"

echo "Stufe 3: members (einzelne Gen-IDs pro Zelle) aus gene_list.tsv"
date
run_stage "3_members" ./target/release/build_members "$OUT/mmseqs_long_gene_list.tsv" "$OUT/paralog"
rm -f "$OUT/mmseqs_long_gene_list.tsv"

echo "Stufe 4: binaere Arrays in eine Container-Datei packen"
date
run_stage "4_pack" ./target/release/pack_container "$OUT/paralog"
rm -f "$OUT/paralog.csc.indptr" "$OUT/paralog.csc.indices" "$OUT/paralog.csc.data"
rm -f "$OUT/paralog.csr.indptr" "$OUT/paralog.csr.indices" "$OUT/paralog.csr.data"
rm -f "$OUT/paralog.members.family_ptr" "$OUT/paralog.members.gene_numbers"

echo "Fertig: $OUT"
date
ls -la "$OUT/"

echo ""
echo "=== Zusammenfassung pro Stufe ==="
printf "%-12s %-20s %-20s\n" "Stufe" "Dauer" "Max RAM (KB)"
for label in "${STAGE_ORDER[@]}"; do
    printf "%-12s %-20s %-20s\n" "$label" "${STAGE_TIME[$label]}" "${STAGE_MEM[$label]}"
done
