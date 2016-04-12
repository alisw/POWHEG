#!/bin/bash

RUNDIR=${PWD}

svn export ../MadGraphStuff MadTMP
#cp -ra ../MadGraphStuff_new MadTMP

if [ -e proc_card.dat ]
then
\cp proc_card.dat MadTMP/Cards/
else
echo proc_card.dat missing!
exit -1
fi

cd MadTMP

cat <<EOF >&2

Warning: as of revision 3087 the BuildMad code has been changed.  It
now modifies the write_proc_labels*.f files, so that they check for
all permutation of final state particles that match the internally
computed one. This generates more flexible code.  If you want the
older behaviour, edit the BuildMad.sh script and delete the lines
indicated after this comment.

Press <enter> to continue ...
EOF

read line

# Modify the sborn_proc and sborn_real so that all fs particle permutations are tried.

ed ./MadGraph_POWHEG/Template/SubProcesses/write_proc_labels.f <<EOF
42
i
      l=2
.
wq
EOF

ed ./MadGraph_POWHEG/Template/SubProcesses/write_proc_labels_real.f <<EOF
41
i
      l=2
.
wq
EOF

# End modification.



./NewProcess.sh $*


\cp -a MadGraph_POWHEG/my_proc/Source/DHELAS .
\cp -a MadGraph_POWHEG/my_proc/Source/MODEL .
\cp -a MadGraph_POWHEG/my_proc/SubProcesses Madlib

\rm Madlib/coupl.inc

cp MODEL/coupl.inc Madlib/

# commented wrt the original version  26/3/2013
#echo editing DHELAS/Makefile
#sed -i 's/FC[\t ]*=.*// ; s/^DEST[\t ]*=.*/DEST = ..\// ;  s/^LIBRARY[\t ]*=.*/LIBRARY = ..\/libdhelas3.a/'   DHELAS/Makefile
#echo editing MODEL/makefile
#sed -i 's/^F77[\t ]*=.*// ; s/-ffixed-line-length-132// ; s/^LIBDIR[\t ]*=.*/LIBDIR = ..\// ' MODEL/makefile
#echo editing Madlib/makefile
#sed -i 's/^F77[\t ]*=.*// ; s/^LIBDIR[\t ]*=.*/LIBDIR = ..\// ' Madlib/makefile


# -n: do not overwrite existing files

echo > ../MGfiles.list

for i in *
do
# -a : exists, file or directory
if ! [ -a ../$i ]
then
    mv $i ../
    echo $i >> ../MGfiles.list
fi
done

cd $RUNDIR
rm -fr MadTMP
