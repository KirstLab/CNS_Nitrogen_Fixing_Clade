#!/bin/bash
#SBATCH --job-name=roast_Nfix
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --time=48:00:00
#SBATCH --output=roast_Nfix_%j.out
#SBATCH --qos=<qos>

pwd; hostname; date
module load cnspipeline

bash ./roast_commands_Nfixing.sh

date
