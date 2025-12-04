#!/bin/bash

#SBATCH --job-name=bakta
#SBATCH --mem=10GB
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -e %A_%a.err
#SBATCH -o %A_%a.out
#SBATCH -c 1
#SBATCH --array=1-942%320

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bakta
outdir="bakta_results"
genome_dir="merged_candidates"

# counter=0
# total=$(find $genome_dir -type f -name *.fna | wc -l)
genomes=($(find $genome_dir -type f -name *.fna | sort))

# Indexing the genomes per batch
start_index=$((($SLURM_ARRAY_TASK_ID - 1) * 10))
end_index=$(($start_index + 9))

for ((i=start_index; i<=end_index && i<${#genomes[@]}; i++)); do
    genome=${genomes[$i]}
    subpath=${genome#$genome_dir/}
    subpath=${subpath%.fna}
    outpath=$outdir/$subpath/
    prefix=$(basename $genome .fna)
    
    echo "Batch $SLURM_ARRAY_TASK_ID: Annotating $prefix"
    PYTHONWARNINGS="ignore" bakta --db ~/endosymb/bakta_db/db \
        --output $outpath \
        --prefix $prefix \
        --threads 1 \
        $genome
done


# for genome in $(find $genome_dir -type f -name *.fna); do
#     ((counter++))
#     subpath=${genome#$genome_dir/}
#     subpath=${subpath%.fna}
#     outpath=$outdir/$subpath/
    
#     prefix=$(basename "$genome")
#     prefix=${prefix%.fna}

#     echo "Annotating $prefix ($counter/$total)..."
#     PYTHONWARNINGS="ignore" bakta --db ~/endosymb/bakta_db/db \
#         --output $outpath \
#         --prefix $prefix \
#         --threads 32 \
#         $genome

# done
# echo "All $total genomes have been annotated."




# export outdir genome_dir
# find $genome_dir -type f -name "*.fna" | parallel -j 32 '
#     file={}
#     subpath=${file#$genome_dir/}
#     subpath=${subpath%.fna}
#     outpath="$outdir/$subpath"

#     echo "Annotating {/.}"
#     PYTHONWARNINGS="ignore" bakta --db ~/endosymb/bakta_db/db \
#             --output $outpath \
#             --prefix {/.} \
#             --threads 1 \
#             $file
 
# '
# echo "All annotations completed."

