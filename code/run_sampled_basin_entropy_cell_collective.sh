#!/bin/bash

## run_sampled_basin_entropy_cell_collective.sh
##
## Bryan Daniels
## 2023/12/8 branched from run_basin_entropy_scan.sh
## 2023/3/31 branched from runFittingAnalysis.sh
## 2021/9/27 branched from runLandauTestSimulations.sh
## 2021/8/4 branched from run_control_kernels.sh
## 2.19.2019 branched from RunScan_BCD.sh
## 10.10.2018 branched from runWormFitting.sh
## 7.18.2018 branched from runWormFitting.sub
## 9.10.2015
## 3.11.2016
## 3.22.2016
## 4.4.2016
##

## 2.19.2019 Run this script using slurm arrays.
## E.g. for 70 networks:
##     sbatch --array=0-69 run_sampled_basin_entropy_cell_collective.sh

#SBATCH -n 1 #40 #20       # number of cores
#SBATCH -t 3-00:00 # wall time (D-HH:MM)
#SBATCH -o slurm.sampled_basin_entropy_HCC.%A_%a.out
#SBATCH -e slurm.sampled_basin_entropy_HCC.%A_%a.err
#SBATCH -q normal

module load anaconda/py3

# you can install required dependencies on the cluster
# in the following way:
#
# > interactive
# > module load anaconda/py3
# > cd ~/InfEst
# > pip install -e .

#export MPLBACKEND=Agg

cd ~/Control-of-Synthetic-Boolean-Networks/code/

python ~/Control-of-Synthetic-Boolean-Networks/code/run_sampled_basin_entropy_cell_collective.py



