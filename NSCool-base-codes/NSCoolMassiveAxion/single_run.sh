#!/bin/sh

argProcess=$1
argEOS=$2
argPairing=$3
argMass=$4
arggann=$5
arggapp=$6
arggaee=$7
arggamm=$8

echo "Print All:"
echo $argProcess $argEOS $argPairing $argMass $arggann $arggapp $arggaee $arggamm
echo "----------------------"

folder=Runs/${argProcess}
mkdir $folder
folder=${folder}/${argEOS}_${argPairing}_${argMass}
mkdir $folder
folder=${folder}/${arggann}_${arggapp}_${arggaee}_${arggamm}/
mkdir  $folder

rm $folder/Cool_Try.in  
rm $folder/I.dat  
rm $folder/Star_Try.dat  
rm $folder/Teff_Try.dat  
rm $folder/Temp_Try.dat
cp -r Model_1/* $folder/.

./NSCool.out $argProcess $argEOS $argPairing $argMass $arggann $arggapp $arggaee $arggamm
python python/NSCool_Likelihood.py $folder

