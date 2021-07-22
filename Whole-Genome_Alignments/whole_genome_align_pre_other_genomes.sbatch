#!/bin/sh
#SBATCH --job-name=whole_genome_align_step1.1 # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user>@ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=15gb                  # Total memory limit
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=whole_genome_align_step1.1_%j.out     # Standard output and error log
#SBATCH --qos=<qos>

# Load required modules
module load cnspipeline
module load samtools/1.9
module load bedtools/2.27.1
module load parallel

pwd; hostname; date

#ln -s ./../genomes/*.fna .

## Prepare the files to be used in the array to run the aligment of each sequence of the genome
THREADS=12
REF_GENOME_NAME=Medicago_truncatula_masked

for sample in `ls *_masked.fna`
do
    QUERY_GENOME=$sample
    python del-query_modified.py $QUERY_GENOME $THREADS $REF_GENOME_NAME
done

date
