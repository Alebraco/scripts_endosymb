#!/bin/bash
#BSUB -W 10:00
#BSUB -J pseudofinder[1-662]
#BSUB -q bobay
#BSUB -e %J_%I.err
#BSUB -o %J_%I.out
#BSUB -R "rusage[mem=20GB] span[hosts=1]"

source ~/.bashrc
conda activate /usr/local/usrapps/metastrain/asoneto/pseudofinder/

blastdb='/gpfs_common/databases/ncbi/blast/nr/nr'

all_files=($(find 'endosymb_only/gbff_files' 'relatives_only/gbff_files' -type f -name "*.gbff" | sort))
start_index=$((($LSB_JOBINDEX - 1) * 10))
end_index=$(($start_index + 9))

for ((i=start_index; i<=end_index && i<${#all_files[@]}; i++)); do
        gbff=${all_files[$i]}
        genome=$(basename $gbff .gbff)
        sp_name=$(basename $(dirname $gbff))
        group=$(basename $(dirname $(dirname $(dirname $gbff))))
        pseudo_dir="${group}/pseudogenes"
        outprefix="${pseudo_dir}/${sp_name}/${genome}"
        mkdir -p $(dirname $outprefix)

        echo "Processing $sp_name - $genome"
        python3 pseudofinder.py annotate -g ${gbff} -db '/gpfs_common/databases/ncbi/blast/nr/nr' -op ${outprefix}
done

