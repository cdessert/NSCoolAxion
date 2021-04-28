#!/bin/sh

argProcess=$1
argEOS=$2
argPairing=$3
argMass=$4
arglogBinit=$5
arggann=$6
arggapp=$7
arggaee=$8
arggamm=$9

echo "Print All:"
echo $argProcess $argEOS $argPairing $argMass $arglogBinit $arggann $arggapp $arggaee $arggamm
echo "----------------------"

folder=Runs/${argProcess}
mkdir $folder
folder=${folder}/${argEOS}_${argPairing}_${argMass}_${arglogBinit}/
mkdir $folder
folder=${folder}/${arggann}_${arggapp}_${arggaee}_${arggamm}/
mkdir  $folder

rm $folder/Cool_Try.in  
rm $folder/I.dat  
rm $folder/Star_Try.dat  
rm $folder/Teff_Try.dat  
rm $folder/Temp_Try.dat
cp -r Model_1/* $folder/.

./NSCool.out $argProcess $argEOS $argPairing $argMass $arglogBinit $arggann $arggapp $arggaee $arggamm
python python/NSCool_Likelihood.py $folder

