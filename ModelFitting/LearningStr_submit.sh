#!/bin/bash
#SBATCH --job-name=MEdecay_Att
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=NONE
#SBATCH --output=/home/mpib/ren/logs/MEdecay_Att_%A_%a.out
#SBATCH --error=/home/mpib/ren/logs/MEdecay_Att_%A_%a.err
#SBATCH --array=1-24

cd rxj-neurocode/HierarchicalCluster/clusterRep_modeling
## Load the relevant modules needed for the job
module load matlab/R2021b

# Choose experiment + mode here
IEXP=1   # 1=ImplicitExp, 2=ExplicitExp, 3=ImplicitRandExp
IMODE=1  # 1=mouse, 2=key

ISUB=${SLURM_ARRAY_TASK_ID}

echo "IEXP=${IEXP} IMODE=${IMODE} ISUB=${ISUB}"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"

## Run your program or script
matlab -nodisplay -nosplash -nodesktop -r "ClusterRep_ChoiceModel_main_HPC(${IEXP}, ${IMODE}, ${ISUB}); exit;"
