#!/bin/bash
#SBATCH --job-name=magrix_formats
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mem=16000m
#SBATCH --output=format_bench_%j.log
#
# Fairer Vergleich MAGrix-Container vs. HDF5 / Zarr / SQLite (bench_formats.py)
# auf Stichproben des echten Datensatzes: jeweils die ersten N Zeilen der
# MMseqs2-TSV werden mit `magrix build` gebaut und dann in alle Formate kopiert.
# Pro Groesse landet eine CSV in results/formats_<N>.csv.
#
# Aufruf: sbatch bash_scripts/run_format_benchmark.sh [input.tsv] [N ...]
# Beispiel:
#   sbatch bash_scripts/run_format_benchmark.sh
#   sbatch bash_scripts/run_format_benchmark.sh /pfad/zu/clusters.tsv 10000000

set -e

cd /pfs/10/work/tu_zxozj42-lm_pangenome/magrix
module load devel/miniforge
eval "$(conda shell.bash hook)"
conda activate pangenome

INPUT_TSV="${1:-/pfs/10/work/tu_zxozj42-lm_pangenome/Large-Matrix/datasets/50_id-80_cov.tsv}"
shift || true
SIZES=("$@")
[ ${#SIZES[@]} -eq 0 ] && SIZES=(10000000 50000000 200000000)

echo "Datensatz : $INPUT_TSV"
echo "Groessen  : ${SIZES[*]} Zeilen"
echo "Knoten    : $(hostname)"
LC_ALL=C lscpu | grep -E "Model name|^Socket|^CPU\(s\)" | sed 's/^/            /'
echo "Commit    : $(git rev-parse --short HEAD 2>/dev/null || echo unbekannt)"

mkdir -p results
for N in "${SIZES[@]}"; do
    W="format_bench/n$N"
    echo ""
    echo "=== Stichprobe: erste $N Zeilen ==="
    date
    rm -rf "$W"
    mkdir -p "$W"
    head -n "$N" "$INPUT_TSV" > "$W/sample.tsv"
    ./target/release/magrix build "$W/sample.tsv" "$W/data"
    rm -f "$W/sample.tsv"
    python3 bench_formats.py "$W/data" "$W/work" --queries 1000 --repeats 3 --csv "results/formats_$N.csv"
    rm -rf "$W"
done

echo ""
echo "Fertig."
date
