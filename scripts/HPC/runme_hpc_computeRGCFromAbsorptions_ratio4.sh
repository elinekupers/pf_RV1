#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time 20:00:00 # time (D-HH:MM)
#SBATCH --mem=250GB # memory pool for all cores
#SBATCH --job-name=linearRGCModel_Ratio4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99/pf_RV1/

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/pf_RV1')); addpath(genpath('/scratch/ek99/JWLOrientedGabor')); addpath(genpath('/scratch/ek99/isetbio')); [expName, subFolder, seed] = prepHPC($SLURM_ARRAY_TASK_ID); linearRGCModel(pfRV1rootPath, subFolder, expName, seed, 4); exit()"

exit

