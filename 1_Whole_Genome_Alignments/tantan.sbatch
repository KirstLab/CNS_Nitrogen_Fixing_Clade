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

# Loads required modules
module load tantan/13

# Removes softmasking of all genomes
for sample in `ls *.fna`
do
base=$(basename $sample ".fna")

awk '{if(/^[^>]/)$0=toupper($0);print $0}' $sample > ${base}_temp.fna
mv ${base}_temp.fna ${base}.fna

done

# Masking of single repeats using tantan

for sample in `ls *.fna`
do
base=$(basename $sample ".fna")

tantan $sample > ${base}_masked.fna

done

date
