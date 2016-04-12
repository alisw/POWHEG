#!/bin/bash


mv powheg.input powheg.input-save


# initial run with default scales
cp powheg.input-save  powheg.input
cat <<'EOF' >> powheg.input

lhrwgt_id 'c'
lhrwgt_descr 'Renorm. scale factor = 1, Factor. scale factor = 1'
lhrwgt_group_name 'scales'
lhrwgt_group_combine 'envelope'

EOF
../pwhg_main

# reweight run with low scales
cp powheg.input-save  powheg.input

sed -i 's/.*renscfact.*/renscfact 0.5/ ; s/.*facscfact.*/facscfact 0.5/' powheg.input

cat <<'EOF' >> powheg.input
compute_rwgt 1

lhrwgt_id 'l'
lhrwgt_descr 'Renorm. scale factor = 0.5, Factor. scale factor = 0.5'
lhrwgt_group_name 'scales'
lhrwgt_group_combine 'envelope'

EOF

../pwhg_main

\mv pwgevents-rwgt.lhe pwgevents.lhe 

# reweight run with high scales
cp powheg.input-save  powheg.input

sed -i 's/.*renscfact.*/renscfact 2/ ; s/.*facscfact.*/facscfact 2/' powheg.input

cat <<'EOF' >> powheg.input
compute_rwgt 1

lhrwgt_id 'h'
lhrwgt_descr 'Renorm. scale factor = 2, Factor. scale factor = 2'
lhrwgt_group_name 'scales'
lhrwgt_group_combine 'envelope'

EOF

../pwhg_main

\mv pwgevents-rwgt.lhe pwgevents.lhe 


mv powheg.input-save powheg.input
