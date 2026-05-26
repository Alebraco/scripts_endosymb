#!/bin/bash
# Step 1 wrapper: marine free-livers feature extraction.
# Chain (all submitted/run from this single bsub job):
#   1. bsub bakta array, wait
#   2. organize_bakta_outputs.sh
#   3. bsub transposase array + collect_features (run_ncbi_features.sh), wait
#   4. make_marine_csv.py  -> files/marine_free_livers.csv  (GATE 1)
#
# Submit from login node: bsub < scripts_endosymb/prediction/run_step1_marine_features.sh
# Assumes Step 0 (download_marine_relatives.py --download) has already populated
# marine_free_livers/genomes/<Species>/<acc>.fna on a login node.

#BSUB -J marine_step1
#BSUB -W 36:00
#BSUB -q bobay
#BSUB -e marine_step1.err
#BSUB -o marine_step1.out
#BSUB -R "rusage[mem=8GB] span[hosts=1]"
#BSUB -n 1

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/annotation

WORKDIR="/rs1/researchers/l/ljbobay/asoneto/endosymb"
BASE="marine_free_livers"
OUTDIR="$BASE/bakta_results"
GENOME_LIST="$WORKDIR/bakta_genome_list_${BASE}.txt"

mkdir -p "$WORKDIR/bakta_logs"

# ---------- 1. Bakta array ----------
if [ ! -f "$GENOME_LIST" ]; then
    find "$WORKDIR/$BASE" -type f -name "*.fna" | sort > "$GENOME_LIST"
fi
NGENOMES=$(wc -l < "$GENOME_LIST")
if [ "$NGENOMES" -eq 0 ]; then
    echo "No .fna files found under $WORKDIR/$BASE. Did Step 0 download succeed?"
    exit 1
fi
NJOBS=$(( (NGENOMES + 99) / 100 ))
echo "Found $NGENOMES marine genomes; submitting $NJOBS bakta array task(s)."

BAKTA_JOB=$(bsub -J "marine_bakta[1-${NJOBS}]%30" \
    -n 8 -R "rusage[mem=20GB] span[hosts=1]" \
    -e "$WORKDIR/bakta_logs/%J_%I.err" -o "$WORKDIR/bakta_logs/%J_%I.out" \
    -W 200:00 -q bobay \
    -env "INPUT_DIR='$WORKDIR/$BASE',OUTDIR='$OUTDIR',GENOME_LIST='$GENOME_LIST',WORKDIR='$WORKDIR'" \
    bash scripts_endosymb/annotation/bakta_annotation.sh | grep -oP '(?<=<)\d+(?=>)')
echo "Submitted bakta array as job $BAKTA_JOB; waiting..."

while true; do
    job=$(bjobs "$BAKTA_JOB" 2>&1)
    if echo "$job" | grep -qE 'RUN|PEND'; then
        sleep 60
        continue
    fi
    if echo "$job" | grep -q 'is not found' || [ -z "$job" ]; then
        echo "Bakta job $BAKTA_JOB not found; assuming complete."
        break
    fi
    break
done
if bjobs -d "$BAKTA_JOB" 2>/dev/null | grep -q 'EXIT'; then
    echo "Bakta array job $BAKTA_JOB exited with failures. Check $WORKDIR/bakta_logs/."
    exit 1
fi
echo "Bakta done."

# ---------- 2. Organize bakta outputs ----------
bash scripts_endosymb/annotation/organize_bakta_outputs.sh "$BASE"

# ---------- 3. Transposase array + collect_features --infer ----------
# run_ncbi_features.sh has its own #BSUB headers; invoking via `bash` just runs
# its body, which submits the transposase array, waits, and runs collect_features.
bash scripts_endosymb/prediction/run_ncbi_features.sh "$BASE"

# ---------- 4. Build files/marine_free_livers.csv (GATE 1) ----------
python -u scripts_endosymb/prediction/make_marine_csv.py --path "$BASE"

echo "Step 1 complete. Review files/marine_free_livers.csv (GATE 1) before submitting Step 2."
