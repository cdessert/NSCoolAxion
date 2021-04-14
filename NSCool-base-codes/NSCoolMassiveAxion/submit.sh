#!/bin/bash

for i in 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13; do
#for i in 1e-13; do
sbatch --export=arggann=$i scan.sh
done

