#!/bin/bash
#SBATCH --job-name=magrix_buffer_sweep
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=8000m
#SBATCH --output=buffer_sweep_%j.log
#
# Measures how the I/O buffer sizes and rows_per_run affect run time and peak
# memory of each pipeline stage. Every configuration is run REPEATS times;
# configurations are interleaved (repeat loop outside) so that slow phases of
# the shared file system don't hit one configuration only.
#
# Results go to <work_dir>/buffer_sweep.csv with the columns
#   experiment,config,repeat,stage,seconds,max_rss_kb,runs
# ("runs" is only filled for the rows_per_run experiment).
#
# The input is read once before the first measurement, so all runs see a warm
# page cache -- mention that when reporting the numbers.
#
# Usage: sbatch bash_scripts/buffer_sweep.sh <input_clusters.tsv> <work_dir> [repeats]
# Example (50M-row sample of the real data):
#   head -n 50000000 datasets/50_id-80_cov.tsv > datasets/sample_50M.tsv
#   sbatch bash_scripts/buffer_sweep.sh datasets/sample_50M.tsv sweep_work 3

set -e

cd /pfs/10/work/tu_zxozj42-lm_pangenome/magrix

if [ $# -lt 2 ]; then
    echo "Usage: sbatch buffer_sweep.sh <input_clusters.tsv> <work_dir> [repeats]"
    exit 1
fi
INPUT_TSV="$1"
WORK="$2"
REPEATS="${3:-3}"
BIN=./target/release
CSV="$WORK/buffer_sweep.csv"

mkdir -p "$WORK"
echo "experiment,config,repeat,stage,seconds,max_rss_kb,runs" > "$CSV"

# Runs one command under /usr/bin/time and appends one CSV line.
# $1 experiment, $2 config, $3 repeat, $4 stage, $5 runs (may be empty), rest: command
measure() {
    local experiment="$1" config="$2" repeat="$3" stage="$4" runs="$5"
    shift 5
    local t
    t=$(mktemp)
    /usr/bin/time -f "%e %M" -o "$t" "$@" > /dev/null 2> "$WORK/last_stderr.txt"
    read -r secs rss < "$t"
    rm -f "$t"
    echo "$experiment,$config,$repeat,$stage,$secs,$rss,$runs" >> "$CSV"
}

# Full pipeline (all four stages) into a fresh directory; deleted afterwards.
full_pipeline() {
    local config="$1" repeat="$2"
    local out="$WORK/run"
    rm -rf "$out"; mkdir -p "$out"
    measure io "$config" "$repeat" 1_main "" $BIN/pangenome_matrices "$INPUT_TSV" "$out"
    local nf nm
    nf=$(($(wc -l < "$out/family_dictionary.tsv") - 1))
    nm=$(($(wc -l < "$out/mag_dictionary.tsv") - 1))
    measure io "$config" "$repeat" 2a_csc "" $BIN/build_csc "$out/mmseqs_long_paralog.tsv" $nf "$out/paralog"
    measure io "$config" "$repeat" 2b_csr "" $BIN/build_csr "$out/mmseqs_long_paralog.tsv" $nm "$out/paralog"
    measure io "$config" "$repeat" 3_members "" $BIN/build_members "$out/mmseqs_long_gene_list.tsv" "$out/paralog"
    measure io "$config" "$repeat" 4_pack "" $BIN/pack_container "$out/paralog"
    rm -rf "$out"
}

echo "Warm-up: reading input once"
date
cat "$INPUT_TSV" > /dev/null

# ── Experiment 1: read/write buffer sizes (MiB), full pipeline ──────────────
# Sizes double from step to step (evenly spaced on a log2 axis). Each config
# changes one buffer and keeps the other at its default (read 64, write 32);
# the default combination itself is run once.
DEFAULT_READ=64
DEFAULT_WRITE=32
IO_CONFIGS=("default $DEFAULT_READ $DEFAULT_WRITE")
for rb in 1 2 4 8 16 32 128 256; do
    IO_CONFIGS+=("read_${rb}MiB $rb $DEFAULT_WRITE")
done
for wb in 1 2 4 8 16 64 128; do
    IO_CONFIGS+=("write_${wb}MiB $DEFAULT_READ $wb")
done

for r in $(seq 1 "$REPEATS"); do
    for line in "${IO_CONFIGS[@]}"; do
        read -r name rb wb <<< "$line"
        echo "[io] repeat $r: $name (read=$rb MiB, write=$wb MiB)"
        MAGRIX_READ_BUF_MB=$rb MAGRIX_WRITE_BUF_MB=$wb full_pipeline "$name" "$r"
    done
done

# ── Experiment 2: rows_per_run for build_csr only ───────────────────────────
echo "Preparing paralog table for the rows_per_run experiment"
date
PREP="$WORK/prep"
rm -rf "$PREP"; mkdir -p "$PREP"
$BIN/pangenome_matrices "$INPUT_TSV" "$PREP" 2> /dev/null
rm -f "$PREP/mmseqs_long_gene_list.tsv"
NM=$(($(wc -l < "$PREP/mag_dictionary.tsv") - 1))

for r in $(seq 1 "$REPEATS"); do
    for rpr in 1000000 2000000 4000000 8000000 16000000 32000000 64000000; do
        echo "[rows_per_run] repeat $r: $rpr"
        # number of runs is printed by build_csr; run once more for it would
        # double the cost, so read it from the stderr of the measured run
        measure rows_per_run "$rpr" "$r" 2b_csr "" $BIN/build_csr "$PREP/mmseqs_long_paralog.tsv" $NM "$PREP/paralog" "$rpr"
        runs=$(grep -o "in [0-9]* runs" "$WORK/last_stderr.txt" | grep -o "[0-9]*")
        sed -i '$ s/,$/,'"$runs"'/' "$CSV"
        rm -f "$PREP"/paralog.csr.*
    done
done
rm -rf "$PREP"

echo "Fertig"
date
echo "Ergebnisse: $CSV"
echo
echo "=== Mittelwerte (Sekunden / Peak-RAM in MB) pro Konfiguration und Stufe ==="
LC_ALL=C awk -F, 'NR > 1 { k = $1 "," $2 "," $4; s[k] += $5; m[k] += $6; n[k]++ }
         END { for (k in s) printf "%-40s %8.2f s %8.1f MB\n", k, s[k]/n[k], m[k]/n[k]/1024 }' "$CSV" | sort
