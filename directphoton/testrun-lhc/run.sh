#!/bin/bash

# First compile the pwhg_main executable in the ../ directory

if [ ! $# -eq 7 ]
then
    echo "Usage: $0 0/1 0/1 0/1 0/1 [number of cores used]"
    echo "       where 0 or 1 indicates if stage 1 2 3 4 is to be started." # 1,2,3 grids; 4 events;
    exit 1
else
    stage1=$1
    stage2=$2
    stage3=$3
    stage4=$4
fi

> timings.txt

PRG=$PWD/../pwhg_main
INPUT=powheg.input-save
xgriditer=4
ncores=$7

# no of PDF error sets
ERRSETS=100

echo -n 'Running stages '
if [ $stage1 -eq 1 ]
then
    echo -n '1 '
fi
if [ $stage2 -eq 1 ]
then
    echo -n '2 '
fi
if [ $stage3 -eq 1 ]
then
    echo -n '3 '
fi
if [ $stage4 -eq 1 ]
then
    echo -n '4 '
fi
echo


if [ $stage1 -eq 1 ] 
then
    echo "stage 1: grids"
    # two stages of importance sampling grid calculation
    for igrid in `seq 1 $xgriditer`
    do
    	(echo -n st1 xg$igrid ' ' ; date ) >> timings.txt

    	cat $INPUT | sed "s/xgriditeration.*/xgriditeration $igrid/ ; s/parallelstage.*/parallelstage 1/" > powheg.input

    	for i in `seq 1 $ncores`
    	do
	    echo $i | $PRG > run-st1-xg$igrid-$i.log 2>&1 &    
    	done
    	wait

    done
fi

if [ $stage2 -eq 1 ]
then
    echo "stage 2: upper bounding function for inclusive cross section"
    # compute NLO and upper bounding envelope for underlying born configurations
    cat $INPUT | sed 's/parallelstage.*/parallelstage 2/' > powheg.input
    (echo -n st2 ' ' ; date ) >> timings.txt
    for i in `seq 1 $ncores`
    do
    	echo $i | $PRG > run-st2-$i.log 2>&1 &
    done
    wait
fi

if [ $stage3 -eq 1 ]
then
    echo "stage 3: upper bounding function for radiation"
    # compute upper bounding coefficients for radiation
    cat $INPUT | sed 's/parallelstage.*/parallelstage 3/' > powheg.input
    (echo -n st3 ' ' ; date ) >> timings.txt
    for i in `seq 1 $ncores`
    do
    	echo $i | $PRG > run-st3-$i.log 2>&1 &
    done
    wait
fi

if [ $stage4 -eq 1 ]
then
    echo "stage 4: event generation"
    # generate events 
    cat $INPUT | sed 's/parallelstage.*/parallelstage 4/' > powheg.input
    (echo -n st4 ' ' ; date ) >> timings.txt
    for i in `seq 1 $ncores`
    do
    	echo $i | $PRG > run-st4-$i.log 2>&1 &
    done
    wait
fi

(echo -n end ' ' ; date ) >> timings.txt

echo Finished.
