#!/bin/sh

argProcess=$1
argEOS=$2
argPairing=$3
argMass=$4
arglogBinit=$5
arglogDeltaM=$6
arggann=$7
arggapp=$8
arggaee=${9}
arggamm=${10}
argma=${11}

echo "Print All:"
echo $argProcess $argEOS $argPairing $argMass $arglogBinit $arglogDeltaM $arggann $arggapp $arggaee $arggamm $argma
echo "----------------------"

folder=$TMPDIR/${argProcess}
mkdir $folder

folder=${folder}/${argEOS}_${argPairing}_${argMass}_${arglogBinit}_${arglogDeltaM}/
mkdir $folder

folder=${folder}/${arggann}_${arggapp}_${arggaee}_${arggamm}_${argma}/
mkdir  $folder

cp -r Model_1/* $folder/.

timeout -k 120 120 ./NSCool.out $argProcess $argEOS $argPairing $argMass $arglogBinit $arglogDeltaM $arggann $arggapp $arggaee $arggamm $argma

#python python/NSCool_Likelihood.py $folder
