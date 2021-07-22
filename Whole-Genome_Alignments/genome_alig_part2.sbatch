#!/bin/sh
#SBATCH --job-name=species # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user>@ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem-per-cpu=2gb           # Memory per processor
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --qos=<qos>
#SBATCH --output=species_%A-%a.out    # Standard output and error log
#SBATCH --array=1-12                 # Array range

# This is an example script that combines array tasks with
# bash loops to process many short runs. Array jobs are convenient
# for running lots of tasks, but if each task is short, they
# quickly become inefficient, taking more time to schedule than
# they spend doing any work and bogging down the scheduler for
# all users.

# Load required modules
module load cnspipeline
module load samtools
module load bedtools/2.27.1
module load parallel                # Array range

pwd; hostname; date

GENOME_NAME=<genome_species>
N_OF_TASKS=12

#Set the number of runs that each SLURM task should do
N_OF_SEQ=$(< ${GENOME_NAME}/list_of_chain.sh wc -l)
N_OF_SEQ_PER_TASK=$(echo "scale=2 ; $N_OF_SEQ / $N_OF_TASKS" | bc)
PER_TASK=$(echo $N_OF_SEQ_PER_TASK | awk '{print int($1+1)}') # round for a value bigger than the division result

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  sed -n "${run}p" ${GENOME_NAME}/list_of_chain.sh | bash
done

date
