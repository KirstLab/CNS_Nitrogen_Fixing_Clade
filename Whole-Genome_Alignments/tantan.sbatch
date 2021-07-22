#!/bin/sh
#SBATCH --job-name=genome_mask_tantan # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user>@ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=12gb                  # Total memory limit
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=genome_mask_tantan_%j.out     # Standard output and error log
#SBATCH --qos=<qos>

pwd; hostname; date

# Load required modules
module load tantan/13

for sample in `ls *.fasta`
do
base=$(basename $sample ".fasta")

tantan $sample > ${base}_masked.fa

done

# remove non-masked fasta files
rm *.fasta
date
