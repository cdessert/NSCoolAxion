#!/bin/bash
#SBATCH --job-name=NSCooling
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2000m
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr3_20
#SBATCH --qos=lr_normal
#SBATCH --mail-type=NONE
#PBATCH END

echo #cid
python load.py $cid
