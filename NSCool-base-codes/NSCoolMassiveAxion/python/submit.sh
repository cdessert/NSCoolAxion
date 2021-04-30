#!/bin/sh

export PBATCH_DIRECTORY='/global/home/users/buschman/pbatch'
alias pbatch=$PBATCH_DIRECTORY'/pbatch.sh'

IDs=12

STR=
for i in `seq 0 64`;do	
	STR="$STR --export=cid=$i load.sh"
done

pbatch $STR

#for i in $(<couplings_mm_sync.dat); do
#pbatch --export=arggamm=$i scan.sh
#done

