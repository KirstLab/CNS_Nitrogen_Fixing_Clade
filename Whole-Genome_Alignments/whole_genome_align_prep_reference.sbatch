#!/bin/sh
#SBATCH --job-name=whole_genome_align_step1 # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user>@ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=10gb                  # Total memory limit
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=whole_genome_align_step1_%j.out     # Standard output and error log
#SBATCH --qos=<qos>

# Load required modules
module load cnspipeline
module load samtools/1.9
module load bedtools/2.27.1
module load parallel

# Run your program with correct path and command line options
## The files must be in the same dir. Use soft links if necessary.

## Prepare the reference genome for the aligments ##
REF_GENOME=Medicago_truncatula_masked.fasta
del-reference.py $REF_GENOME
