#!/bin/bash
#SBATCH --job-name=MEdecay
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=NONE
#SBATCH --output=/home/mpib/ren/logs/MEdecay_%A_%a.out
#SBATCH --error=/home/mpib/ren/logs/MEdecay_%A_%a.err
#SBATCH --array=1-3

cd rxj-neurocode/HierarchicalCluster/clusterRep_modeling
## Load the relevant modules needed for the job
module load matlab/R2021b

# Choose expMode here: 1=mouse, 2=key
IMODE=1

echo "Running SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, IMODE=${IMODE}"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"

## Run your program or script
matlab -nodisplay -nosplash -nodesktop -r "ClusterRep_ChoiceModel_main_HPC(${SLURM_ARRAY_TASK_ID}, ${IMODE}); exit;"


