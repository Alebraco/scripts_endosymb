#!/bin/bash
#BSUB -J "bakta_batch[1-1000]%500"
#BSUB -n 1
#BSUB -R "rusage[mem=10GB]"
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out
#BSUB -W 96:00

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/bakta
outdir="bakta_results"
genome_dir="candidate_genomes"

genomes=($(find $genome_dir -type f -name "*.fna" | sort))

start_index=$((($LSB_JOBINDEX - 1) * 10))
end_index=$(($start_index + 9))

for ((i=start_index; i<=end_index && i<${#genomes[@]}; i++)); do
    genome=${genomes[$i]}
    subpath=${genome#$genome_dir/}
    subpath=${subpath%.fna}
    outpath=$outdir/$subpath/
    prefix=$(basename $genome .fna)

    echo "Batch $LSB_JOBINDEX: Annotating $prefix"
    PYTHONWARNINGS="ignore" bakta --db /rs1/researchers/l/ljbobay/asoneto/endosymb/bakta_db \
        --output $outpath \
        --prefix $prefix \
        --threads 1 \
        $genome
done

conda deactivate