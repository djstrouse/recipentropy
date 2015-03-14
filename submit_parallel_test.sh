#!/bin/bash

CODEDIR=/scratch/network/astrandb/recipentropy
BASEDIR=/scratch/network/astrandb/recipentropy

declare subsamps=(1 2)
declare orders=(0 1)

for o in "${orders[@]}" do
	for i in "${subsamps[@]}" do
		echo "addpath(genpath('/scratch/network/astrandb/recipentropy')); a = $i; c = $o; test(a,c,c);" >> ${CODEDIR}/wrapper_$i_$o.m
		cat pbs.info > ${BASEDIR}/$i_$o.sub
		echo "/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -nojvm -r wrapper_$i_$i_$o" >> ${BASEDIR}/$i_$o.sub
		echo "exit" >> ${BASEDIR}/$i_$o.sub
		qsub ${BASEDIR}/$i_$o.sub
	done
done