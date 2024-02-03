#!/bin/bash

## run_control_kernels_iowa.sh
##
## Bryan Daniels
## 2024/2/2 branched from run_control_kernels.sh
## 2.19.2019 branched from RunScan_BCD.sh
## 10.10.2018 branched from runWormFitting.sh
## 7.18.2018 branched from runWormFitting.sub
## 9.10.2015
## 3.11.2016
## 3.22.2016
## 4.4.2016
##

## 2.19.2019 Run this script using slurm arrays. E.g. for 87 networks,
##     sbatch --array=0-86 run_control_kernels.sh

#SBATCH -n 1 #40 #20       # number of cores
#SBATCH -t 1-00:00 # wall time (D-HH:MM)
#SBATCH -o slurm.control_kernels.%A_%a.out
#SBATCH -e slurm.control_kernels.%A_%a.err
#SBATCH -q normal

module load anaconda/py3

cd ~/Control-of-Synthetic-Boolean-Networks/code/

python ~/Control-of-Synthetic-Boolean-Networks/code/run_control_kernels_iowa.py $SLURM_ARRAY_TASK_ID


