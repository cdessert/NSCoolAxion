#!/bin/bash
#SBATCH --job-name=NSCooling
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2000m
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr3
#SBATCH --qos=lr_normal
#SBATCH --mail-type=NONE
#PBATCH END


#argProcess=65504 #nucleons
#argProcess=65536 #muon

argProcess=1 #m sync
for arglogBinit in 11 12 13 14 15 16;do
for argPairing in 0 1 2 3 4 5; do
for argEOS in 0 1 2 3 4; do
#for argPairing in 0; do
for argMass in 1.0 1.2 1.4 1.6 1.8 2.0; do
#arggapp=$arggann
arggann=0
arggapp=0
arggaee=0
#arggamm=0

./single_run.sh $argProcess $argEOS $argPairing $argMass $arglogBinit $arggann $arggapp $arggaee $arggamm &

done
done 
wait
done
done

