#!/bin/bash
#SBATCH --clusters=arc
#SBATCH --partition=short
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=01:00:00
#SBATCH --job-name=AbaqusJob
#SBATCH --switch=1

module purge
module load Abaqus/2022
module load iimpi/2020a

. abaqus_arcmpi.sh

## abaqus fetch job=model1.inp
abaqus input=BrainFull.inp job=AbaqusJob user=UMAT_new.for cpus=${SLURM_NTASKS} interactive
