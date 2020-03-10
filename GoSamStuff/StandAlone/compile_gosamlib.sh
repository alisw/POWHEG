#!/bin/sh

cd $1/../GoSamlib/

FILES1=$(echo *.f *.f90 ' ' | sed 's/qlonshellcutoff.f// ; s/qlconstants.f// ; s/.f /.o /g ; s/.f90 /.o /g ; ')
FILES2=$(echo *.cc *.F90 ' ' | sed 's/.cc /.o /g ; s/.F90 /.o /g ; ')

cd $1

ar cru gosamlib.a $FILES1 $FILES2
