#!/usr/bin/env bash

VP=(1 2 3 4)
VT=(40 40 40 40)
np=4
VE=(8 16 32 64 128)
ne=5

Fold=1

for (( ie = 0; ie < ne; ie+=1 ))
do
	for (( ip=0; ip<np; ip+=1 ))
	do
		let NbTimeStep=${VT[$ip]}*Fold
		pref="p${VP[ip]}ne${VE[$ie]}nt${NbTimeStep}"
        
        cat dg.h | \
        sed "s/NbSideElem/${VE[$ie]}/g" | \
        sed "s/OrderOfApprox/${VP[$ip]}/g" > ${pref}.h
        
        cat dg.c | \
        sed "s/NbTimeStep/${NbTimeStep}/g" | \
        sed "s/OwnHeader/${pref}.h/g" > ${pref}.c

        cat submit.pbs | \
        sed "s/TEMPLATE_CASE/${pref}/g" > ${pref}.pbs

        gcc ${pref}.c -lm -o ${pref}.out

        echo "qsub ${pref}.pbs" >> all.pbs
        echo "rm -f ${pref}.c" >> clean.sh
        echo "rm -f ${pref}.h" >> clean.sh
        echo "rm -f ${pref}.pbs" >> clean.sh
        echo "rm -f ${pref}.out" >> clean.sh
	done
	let Fold=Fold*2
done
