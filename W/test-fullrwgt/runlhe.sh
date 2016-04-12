#!/bin/bash

PRG=../../lhef_analysis

for dir in testrun-wm-lhc-8TeV-ct10 testrun-wm-lhc-8TeV-MRST
do

cd $dir

for i in {1..64}
do
case $i in
?) ch=000$i ;;
??) ch=00$i ;;
???) ch=0$i ;;
????) ch=$i ;;
esac
#name: pwgevents-rwgt-0034.lhe
echo pwgevents-$ch.lhe | $PRG &
#case $i in
#9|17|25|33|41|49|57) wait;;
#esac


done

cd ../

done
