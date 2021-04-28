#!/bin/bash

#for i in $(<couplings_NN.dat); do
#sbatch --export=arggann=$i scan.sh
#done


#for i in $(<couplings_mm.dat); do
#sbatch --export=arggamm=$i scan.sh
#done

for i in $(<couplings_mm_sync.dat); do
sbatch --export=arggamm=$i scan.sh
done

#for i in 1e-10 7.8e-11 6.2e-11 4.8e-11 3.8e-11 3e-11 2.3e-11 1.8e-11 1.4e-11 1.1e-11 1e-11 8.9e-12 7e-12 5.5e-12 4.3e-12 3.4e-12 2.6e-12 2.1e-12 1.6e-12 1.3e-12 1e-12 -1e-10 -7.8e-11 -6.2e-11 -4.8e-11 -3.8e-11 -3e-11 -2.3e-11 -1.8e-11 -1.4e-11 -1.1e-11 -1e-11 -8.9e-12 -7e-12 -5.5e-12 -4.3e-12 -3.4e-12 -2.6e-12 -2.1e-12 -1.6e-12 -1.3e-12 -1e-12; do
#sbatch --export=arggamm=$i scan.sh
#done

#for i in 1e-9 7.8e-10 6.2e-10 4.8e-10 3.8e-10 3e-10 2.3e-10 1.8e-10 1.4e-10 1.1e-10 1e-10 -1e-9 -7.8e-10 -6.2e-10 -4.8e-10 -3.8e-10 -3e-10 -2.3e-10 -1.8e-10 -1.4e-10 -1.1e-10 -1e-10; do
#sbatch --export=arggamm=$i scan.sh
#done


#for i in -1e-7 -1e-8 -1e-9 -1e-10 -1e-11 -1e-12 -1e-13 0 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13; do
#sbatch --export=arggamm=$i scan.sh
#done

