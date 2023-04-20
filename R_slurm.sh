#!/bin/bash

#SBATCH --job-name=R_CompDist2axial
#SBATCH --output=log.%A_%a           # %J is the Slurm's JOBID
#SBATCH --error=log.%A_%a
#SBATCH --array=1-21              # array job with 21 array elements
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8       #number of parallel R workers +1 


##SBATCH --partition=TEST
##SBATCH --time=00:05:00           # HH:MM:SS
##SBATCH --mem=16G          

#SBATCH --partition=CPU-long       
#SBATCH --time=2-00:00:00           # D-HH:MM:SS
#SBATCH --mem=128G              

module load gcc

Rscript CompDist2VonMisesApril2022_cluster_23_03.R $SLURM_ARRAY_TASK_ID
