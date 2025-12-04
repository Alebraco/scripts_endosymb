#!/bin/bash
#SBATCH --job-name=download_nt
#SBATCH --output=download_nt_%j.log
#SBATCH --error=download_nt_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

echo "Starting NT database download at $(date)"
LD_LIBRARY_PATH= update_blastdb.pl --decompress nt

echo "Download completed at $(date)"
echo "Testing database installation"
blastdbcmd -db nt -info
