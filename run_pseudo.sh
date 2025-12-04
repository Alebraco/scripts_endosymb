#!/bin/bash
#BSUB -W 200:00
#BSUB -J "pseudofinder[1-662]%30"
#BSUB -q bobay
#BSUB -e logs/%J_%I.err
#BSUB -o logs/%J_%I.out
#BSUB -R "rusage[mem=32GB] span[hosts=1]"
#BSUB -n 8

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/pseudofinder/
mkdir -p logs

blastdb='/gpfs_common/databases/ncbi/blast/nr/nr'
pseudopy='software/pseudofinder/pseudofinder.py'
all_files=($(find 'endosymb_only/gbff_files' 'relatives_only/gbff_files' -type f -name "*.gbff" | sort))
start_index=$((($LSB_JOBINDEX - 1) * 10))
end_index=$(($start_index + 9))

for ((i=start_index; i<=end_index && i<${#all_files[@]}; i++)); do
        gbff=${all_files[$i]}
        genome=$(basename $gbff .gbff)
        sp_name=$(basename $(dirname $gbff))
        group=$(basename $(dirname $(dirname $(dirname $gbff))))
        pseudo_dir="/share/metastrain/asoneto/endosymb/pseudofinder/${group}/pseudogenes"
        outprefix="${pseudo_dir}/${sp_name}/${genome}"
        mkdir -p $(dirname $outprefix)

        echo "Processing $sp_name - $genome"
        python3 $pseudopy annotate -g ${gbff} -db '/gpfs_common/databases/ncbi/blast/nr/nr' -op ${outprefix} -di --threads 8
done
