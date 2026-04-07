#!/bin/bash

# Worker script for downloading and extracting ATB Bakta annotations.
# Called as a SLURM array job by submit_atb_bakta.sh — do not run directly.
#
# Each task handles one bakta tar.xz batch:
#   1. Look up this task's tar filename and URL from the batch index.
#   2. Download the tar.xz (skipped if already present).
#   3. For each sample in this batch's list, extract SAMPLE.bakta.json from
#      the tar using a wildcard match (path inside the tar is unknown).
#      Falls back to SAMPLE.bakta.json.gz if the plain JSON is not found.
#   4. Delete the tar.xz after extraction to save disk space.
#
# Output: annotations are written to ANNOT_DIR/Species_name/SAMPLE.bakta.json
#
# SLURM_ARRAY_TASK_ID is set automatically by SLURM and corresponds to the
# 1-based line number in bkt_batch_index.tsv.

# Must match the OUTDIR set in submit_atb_bakta.sh
OUTDIR="atb_bacteria"

BATCH_INDEX="${OUTDIR}/bkt_batch_index.tsv"   # idx → tar_xz, url
BATCH_LISTS="${OUTDIR}/batch_lists"            # per-batch sample lists
DL_DIR="${OUTDIR}/downloads"                   # temporary tar.xz storage
ANNOT_DIR="${OUTDIR}/annotations"             # final annotation JSON output

mkdir -p "${DL_DIR}" "${ANNOT_DIR}"

# Read the row for this array task from the index file.
# Format: batch_idx <TAB> tar_xz_filename <TAB> download_url
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BATCH_INDEX}")
BATCH_IDX=$(echo "$LINE" | cut -f1)
TAR_XZ=$(echo "$LINE" | cut -f2)
URL=$(echo "$LINE" | cut -f3)

if [[ -z "$TAR_XZ" || -z "$URL" ]]; then
    echo "ERROR: No batch info for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

TARBALL="${DL_DIR}/${TAR_XZ}"

# Download the batch tar.xz if not already on disk.
# Write to a .tmp file first so an interrupted download never leaves a
# partial file that looks complete.
if [[ ! -f "${TARBALL}" ]]; then
    echo "Downloading: ${TAR_XZ}"
    wget -q -O "${TARBALL}.tmp" "${URL}"
    if [[ $? -eq 0 && -s "${TARBALL}.tmp" ]]; then
        mv "${TARBALL}.tmp" "${TARBALL}"
    else
        echo "ERROR: download failed for ${TAR_XZ}"
        rm -f "${TARBALL}.tmp"
        exit 1
    fi
else
    echo "Already downloaded: ${TAR_XZ}"
fi

# Per-batch list of samples to extract.
# Columns: sample, sylph_species
BATCH_LIST="${BATCH_LISTS}/bkt_batch_${BATCH_IDX}.tsv"
EXTRACTED=0
SKIPPED=0

while IFS=$'\t' read -r SAMPLE SPECIES; do
    # Build a safe output directory name from the species string
    SPECIES_DIR="${ANNOT_DIR}/$(echo "${SPECIES}" | tr ' ' '_' | tr '/' '_' | tr -d "'")"
    DEST="${SPECIES_DIR}/${SAMPLE}.bakta.json"

    # Skip if already extracted (allows safe re-runs)
    if [[ -f "${DEST}" ]]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    mkdir -p "${SPECIES_DIR}"

    # The exact path inside the tar is unknown, so use --wildcards to match
    # any directory prefix. Stream directly to DEST (-O = stdout).
    tar -xJf "${TARBALL}" --wildcards "*/${SAMPLE}.bakta.json" -O > "${DEST}" 2>/dev/null
    if [[ $? -eq 0 && -s "${DEST}" ]]; then
        EXTRACTED=$((EXTRACTED + 1))
    else
        # Some ATB releases store annotations gzip-compressed inside the tar.
        # Try the .json.gz variant and decompress after extraction.
        tar -xJf "${TARBALL}" --wildcards "*/${SAMPLE}.bakta.json.gz" -O > "${DEST}.gz" 2>/dev/null
        if [[ $? -eq 0 && -s "${DEST}.gz" ]]; then
            gzip -d "${DEST}.gz"
            EXTRACTED=$((EXTRACTED + 1))
        else
            echo "WARN: extraction failed for ${SAMPLE} from ${TAR_XZ}"
            rm -f "${DEST}" "${DEST}.gz"
        fi
    fi
done < "${BATCH_LIST}"

echo "Batch ${TAR_XZ}: extracted=${EXTRACTED} skipped=${SKIPPED}"

# Remove the tarball after extraction to free disk space.
# Re-running the job will re-download it if needed.
rm -f "${TARBALL}"
