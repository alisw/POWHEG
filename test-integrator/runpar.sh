#!/bin/bash

> Timings.txt

ncores=4
nprocesses=8
# First compile the pwhg_main executable in the ../ directory
#

# the following function limits the number of subprocesses
# to be not larger than the number of cores specified by the
# user

function limit_procs {
    while [ `jobs -p | wc -w` -gt $ncores ]
    do
	sleep 1
    done
}

PWHGMAIN=../pwhg_main

# two stages of importance sampling grid calculation
for igrid in 1 2
do

(echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt

cat powheg.input-save | sed "s/xgriditeration.*/xgriditeration $igrid/ ; s/parallelstage.*/parallelstage 1/" > powheg.input

for i in `seq $nprocesses`
do
    echo $i | $PWHGMAIN > run-st1-xg$igrid-$i.log 2>&1 &
    limit_procs
done
wait

done



# compute NLO and upper bounding envelope for underlying born comfigurations
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 2/ ' > powheg.input
(echo -n st2 ' ' ; date ) >> Timings.txt
for i in `seq $nprocesses`
do
    echo $i | $PWHGMAIN > run-st2-$i.log 2>&1 &
    limit_procs
done
wait


# compute upper bounding coefficients for radiation
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 3/' > powheg.input
(echo -n st3 ' ' ; date ) >> Timings.txt
for i in `seq $nprocesses`
do
    echo $i | $PWHGMAIN > run-st3-$i.log 2>&1 &
    limit_procs
done
wait



# generate events 
cat powheg.input-save | sed 's/parallelstage.*/parallelstage 4/' > powheg.input
(echo -n st4 ' ' ; date ) >> Timings.txt
for i in `seq $nprocesses`
do
    echo $i | $PWHGMAIN > run-st4-$i.log 2>&1 &
    limit_procs
done
wait

(echo -n end ' ' ; date ) >> Timings.txt




