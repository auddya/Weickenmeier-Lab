#!/bin/bash
#SBATCH --clusters=arc
#SBATCH --partition=devel
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:10:00
#SBATCH --job-name=AbaqusJob2
#SBATCH --switch=1

module purge
module load Abaqus/2022
module load iimpi/2020a

. abaqus_arcmpi.sh

## abaqus fetch job=model1.inp
abaqus input=BrainFull3.inp job=AbaqusJob2 user=UMAT_new.for cpus=${SLURM_NTASKS} interactive
