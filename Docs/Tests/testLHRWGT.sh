#!/bin/bash

# run the program the old fashion way

\rm Makefile

svn up Makefile

sed -i 's/PDF=lhapdf/PDF=native/
        s/COMPILER=ifort/COMPILER=gfortran/
        s/ANALYSIS=.*/ANALYSIS=none/' Makefile


make clean
make -j > make-orig.log 2>&1

mv pwhg_main  pwhg_main-orig

make clean

sed -i \
's,VPATH= \./,VPATH= ./:../StandardRW/,
\,INCLUDE1=$(PWD)/include, aINCLUDE1.5=\$(shell dirname \$(PWD))/StandardRW/include
s/-I\$(INCLUDE1)/-I\$(INCLUDE1) -I\$(INCLUDE1.5)/ ' Makefile

make -j > make-lhrwgt.log 2>&1
mv pwhg_main  pwhg_main-lhrwgt

rm -rf testLHRWGT-1 testLHRWGT-2

mkdir testLHRWGT-1
mkdir testLHRWGT-2

cp testrun-lhc/*  testLHRWGT-1/
cp testrun-lhc/*  testLHRWGT-2/

#************************************ Run with the same setup, no rwgt flags
cd testLHRWGT-1
sed -i 's/^ *numevts .*/numevts 50/ ; s/^ *ncall1 .*/ncall1 2000/
        s/^ *itmx1 .*/itmx1 5/ ; s/^ *ncall2 .*/ncall2 2000/
        s/^ *itmx2 .*/itmx2 5/ ; s/^ *nubound .*/nubound 5000/
        s/^ *foldcsi .*/foldcsi 1/;  s/^ *foldy .*/foldy 1/;  s/^ *foldphi .*/foldphi 1/ ' powheg.input
../pwhg_main-orig > run.log 2>&1 &

cp powheg.input ../testLHRWGT-2/
cd ../testLHRWGT-2/
../pwhg_main-lhrwgt > run.log 2>&1 &

wait

cd ..
echo "we have now ran the old and new code with no reweight info"
echo "check that the lhe files are the same"
echo diff -b testLHRWGT-1/pwgevents.lhe testLHRWGT-2/pwgevents.lhe
diff -b testLHRWGT-1/pwgevents.lhe testLHRWGT-2/pwgevents.lhe
echo "should generate no output"
echo also
echo grep rwgt testLHRWGT-1/pwgevents.lhe
grep rwgt testLHRWGT-1/pwgevents.lhe
echo "should generate output (storeinfo_rwgt defaults to 1)"
read line

#************************************ Run with the same setup, storeinfo_rwgt 1
cd testLHRWGT-1
cat <<EOF >> powheg.input
storeinfo_rwgt 1
EOF
\rm pwgevents.lhe
../pwhg_main-orig > run.log 2>&1 &

cd ../testLHRWGT-2/
\cp ../testLHRWGT-1/powheg.input .
\rm pwgevents.lhe
../pwhg_main-lhrwgt > run.log 2>&1 &

wait

cd ..
echo "we have now ran the old and new code with storeinfo_rwgt 1"
echo "check that the lhe files are the same"
echo diff -b testLHRWGT-1/pwgevents.lhe testLHRWGT-2/pwgevents.lhe
diff -b testLHRWGT-1/pwgevents.lhe testLHRWGT-2/pwgevents.lhe
echo "should generate no output"
read line

#************************************ Run with the same setup, compute_rwgt 1
cd testLHRWGT-1/

sed -i 's/^storeinfo_rwgt 1/compute_rwgt 1/
        s/[# ]*renscfact .*/renscfact 0.5d0/
        s/[# ]*facscfact .*/facscfact 0.5d0/' powheg.input

../pwhg_main-orig > run1.log 2>&1 &

cd ../testLHRWGT-2/

\cp ../testLHRWGT-1/powheg.input .
../pwhg_main-lhrwgt > run1.log 2>&1 &
wait

cd ..
echo "we have now ran the old and new code with compute_rwgt 1"
echo "check that the lhe files are the same"

echo diff -b testLHRWGT-1/pwgevents-rwgt.lhe testLHRWGT-2/pwgevents-rwgt.lhe
diff -b testLHRWGT-1/pwgevents-rwgt.lhe testLHRWGT-2/pwgevents-rwgt.lhe
echo "should generate no output, except for the missing"
echo "</LesHouchesEvents> in the old implementation"
read line
#************************************ Run new code, compute_rwgt 1,
#   lhrwgt_group_name 'scales' lhrwgt_id 'c' lhrwgt_descr 'central'
cd testLHRWGT-2/

sed -i 's/^compute_rwgt 1/storeinfo_rwgt 1/
        s/[# ]*renscfact .*/# renscfact 1/
        s/[# ]*facscfact .*/# facscfact 1/' powheg.input

cat <<EOF >> powheg.input
lhrwgt_group_name 'scales'
lhrwgt_id 'c'
lhrwgt_descr 'central'
EOF

\rm pwgevents.lhe
../pwhg_main-lhrwgt > run.log 2>&1 &
wait

cat <<'EOF'
We have run the program with:
storeinfo_rwgt 1
lhrwgt_group_name 'scales'
lhrwgt_id 'c'
lhrwgt_descr 'central'
EOF
cd ..
echo "check the lhe files are different"
echo diff -b testLHRWGT-1/pwgevents.lhe testLHRWGT-2/pwgevents.lhe
diff -b testLHRWGT-1/pwgevents.lhe testLHRWGT-2/pwgevents.lhe
echo "should generate output"
read line
#************************************
cd testLHRWGT-2/

sed -i 's/^storeinfo_rwgt 1/compute_rwgt 1/
        s/[# ]*renscfact .*/renscfact 0.5d0/
        s/[# ]*facscfact .*/facscfact 0.5d0/
        s/lhrwgt_id .*/lhrwgt_id "l"/
        s/lhrwgt_descr .*/lhrwgt_descr "low"/' powheg.input

../pwhg_main-lhrwgt > run.log 2>&1 &
wait

cat <<'EOF'
We have run the program with:
compute_rwgt 1
lhrwgt_group_name 'scales'
lhrwgt_id 'c'
lhrwgt_descr 'central'
EOF
cd ..
echo "check the lhe files are different"
echo diff -b testLHRWGT-1/pwgevents-rwgt.lhe testLHRWGT-2/pwgevents-rwgt.lhe
diff -b testLHRWGT-1/pwgevents-rwgt.lhe testLHRWGT-2/pwgevents-rwgt.lhe
echo "should generate output"
