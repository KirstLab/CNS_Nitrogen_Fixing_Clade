#!/bin/sh
#SBATCH --job-name=species # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<users>  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=24            # Number of CPU cores per task
#SBATCH --mem=40gb                  # Total memory limit
#SBATCH --time=96:00:00              # Time limit hrs:min:sec
#SBATCH --output=nspecies_%j.out     # Standard output and error log
#SBATCH --qos=kirst

# Load required modules
module load gcc/5.2.0
module load cnspipeline
module load samtools/1.9
module load bedtools/2.27.1
module load parallel

THREADS=24
REF_GENOME_NAME=Medicago_truncatula_masked

for sample in `ls ./*.fna`
do
base=$(basename $sample ".fna")

QUERY_GENOME=${base}.fna

python /home/wendellpereira/bin/CNSpipeline/del-query_chaining.py $QUERY_GENOME $THREADS $REF_GENOME_NAME

done
