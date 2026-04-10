#!/bin/bash
# Worker script for downloading and extracting ATB data (assemblies or bakta).
# Called as a SLURM array job by submit_atb.sh — do not run directly.
#
# Usage: run_atb_worker.sh <assemblies|bakta>
#
# Each task handles one tar.xz batch:
#   1. Look up this task's tar filename and URL from the batch index.
#   2. Download the tar.xz (skipped if already present).
#   3. Extract the target files for each sample in this batch.
#   4. Delete the tar.xz after extraction to save disk space.
#
# assemblies: extracts SAMPLE.fna to genomes/Species_name/
# bakta:      extracts SAMPLE.bakta.json to annotations/Species_name/


MODE="${MODE:?MODE environment variable must be set to 'assemblies' or 'bakta'}"

OUTDIR="atb_bacteria"
MANIFEST="${OUTDIR}/manifest.tsv"
DL_DIR="${OUTDIR}/downloads"

case "$MODE" in
    assemblies)
        BATCH_INDEX="${OUTDIR}/asm_batch_index.tsv"
        OUT_DIR="${OUTDIR}/genomes"
        ;;
    bakta)
        BATCH_INDEX="${OUTDIR}/bkt_batch_index.tsv"
        OUT_DIR="${OUTDIR}/annotations"
        ;;
    *)
        echo "ERROR: Unknown mode '${MODE}'. Use 'assemblies' or 'bakta'."
        exit 1
        ;;
esac

mkdir -p "${DL_DIR}" "${OUT_DIR}"

# Read the row for this array task from the index file.
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BATCH_INDEX}")
TAR_XZ=$(echo "$LINE" | cut -f1)
URL=$(echo "$LINE" | cut -f2)

if [[ -z "$TAR_XZ" || -z "$URL" ]]; then
    echo "ERROR: No batch info for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

TARBALL="${DL_DIR}/${TAR_XZ}"

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

EXTRACTED=0
SKIPPED=0

if [[ "$MODE" == "assemblies" ]]; then
    while IFS=$'\t' read -r SAMPLE SPECIES FILE_IN_TAR; do
        SPECIES_DIR="${OUT_DIR}/$(echo "${SPECIES}" | tr ' ' '_' | tr '/' '_' | tr -d "'")"
        DEST="${SPECIES_DIR}/${SAMPLE}.fna"

        if [[ -f "${DEST}" ]]; then
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        mkdir -p "${SPECIES_DIR}"

        tar -xf "${TARBALL}" "${FILE_IN_TAR}" -O > "${DEST}"
        if [[ $? -eq 0 && -s "${DEST}" ]]; then
            EXTRACTED=$((EXTRACTED + 1))
        else
            echo "WARN: extraction failed for ${SAMPLE} from ${TAR_XZ}"
            rm -f "${DEST}"
        fi
    done < <(awk -F'\t' -v tar="$TAR_XZ" 'NR>1 && $4==tar {print $1"\t"$2"\t"$3}' "$MANIFEST")

else 
    while IFS=$'\t' read -r SAMPLE SPECIES; do
        SPECIES_DIR="${OUT_DIR}/$(echo "${SPECIES}" | tr ' ' '_' | tr '/' '_' | tr -d "'")"
        DEST="${SPECIES_DIR}/${SAMPLE}.bakta.json"

        if [[ -f "${DEST}" ]]; then
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        mkdir -p "${SPECIES_DIR}"

        tar -xf "${TARBALL}" "${SAMPLE}.bakta.json" -O > "${DEST}"
        if [[ $? -eq 0 && -s "${DEST}" ]]; then
            EXTRACTED=$((EXTRACTED + 1))
        else
            echo "WARN: extraction failed for ${SAMPLE} from ${TAR_XZ}"
            rm -f "${DEST}"
        fi
    done < <(awk -F'\t' -v tar="$TAR_XZ" 'NR>1 && $6==tar {print $1"\t"$2}' "$MANIFEST")
fi

echo "Batch ${TAR_XZ}: extracted=${EXTRACTED} skipped=${SKIPPED}"

rm -f "${TARBALL}"
