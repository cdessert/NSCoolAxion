#!/bin/bash
#SBATCH --job-name=NSCooling
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:05:00
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#PBATCH END

if [ -z "${TMPDIR+x}" ]; then
    echo "TMPDIR is not set"
elif [ -z "$TMPDIR" ]; then
    echo "TMPDIR is set but empty"
else

    for argProcess in 4259808; do
    for arglogBinit in 11;do
    for argPairing in 0; do
    for argEOS in 0; do

        for argMass in 1.0; do
        for arglogDeltaM in -20; do 

            arggaee=0
            arggamm=0
            argma=0
            arggapp=0
            arggann=0

            ./single_run.sh $argProcess $argEOS $argPairing $argMass $arglogBinit $arglogDeltaM $arggann $arggapp $arggaee $arggamm $argma &

        done
        done
        wait

        echo "Start collecting data:"
        ## Collect data here

        rm -rf $TMPDIR/*

    done
    done
    done
    done

fi
