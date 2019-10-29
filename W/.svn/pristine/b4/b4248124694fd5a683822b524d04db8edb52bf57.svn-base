#!/bin/bash

> Timings.txt


PROG=../../pwhg_main

# First compile the pwhg_main executable in the ../ directory
#

# two stages of importance sampling grid calculation
for igrid in {1..2}
do

(echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt

cat powheg.input-save | sed "s/xgriditeration.*/xgriditeration $igrid/ ; s/parallelstage.*/parallelstage 1/" > powheg.input

for i in {1..64}
do
echo $i | $PROG > run-st1-xg$igrid-$i.log 2>&1 &
done
wait

done



# compute NLO and upper bounding envelope for underlying born comfigurations
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 2/ ' > powheg.input
(echo -n st2 ' ' ; date ) >> Timings.txt
for i in {1..64}
do
echo $i | $PROG > run-st2-$i.log 2>&1 &
done
wait


# compute upper bounding coefficients for radiation
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 3/' > powheg.input
(echo -n st3 ' ' ; date ) >> Timings.txt
for i in {1..64}
do
echo $i | $PROG > run-st3-$i.log 2>&1 &
done
wait



# generate events 
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 4/' > powheg.input
(echo -n st4 ' ' ; date ) >> Timings.txt
for i in {1..64}
do
echo $i | $PROG > run-st4-$i.log 2>&1 &
done
wait

(echo -n end ' ' ; date ) >> Timings.txt




# reweighting

cat powheg.input-save | sed 's/ 21100 .*/ 10800  ! CT10/' | grep -v 'lhrwgt' > powheg.input

cat <<EOF >> powheg.input

lhrwgt_group_name 'PDF reweighting'
lhrwgt_id 'ct10rw'
lhrwgt_descr 'MRST reweighted to ct10'

compute_rwgt 1
fullrwgt 1
fullrwgtmode 4

EOF

# generate reweighted events 
(echo -n rwgt ' ' ; date ) >> Timings.txt
for i in {1..64}
do
case $i in
?) ch=000$i ;;
??) ch=00$i ;;
???) ch=0$i ;;
????) ch=$i ;;
esac

(echo $i ; echo pwgevents-$ch.lhe) | $PROG > run-rwgt-$i.log 2>&1 &
done
wait

(echo -n end ' ' ; date ) >> Timings.txt
