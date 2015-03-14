#!/bin/sh

# submit_parallel_test_copy.sh
# 
#
# Created by Ari Strandburg-Peshkin on 8/23/13.
# Copyright 2013 __MyCompanyName__. All rights reserved.


CODEDIR=/scratch/network/astrandb/recipentropy
BASEDIR=/scratch/network/astrandb/recipentropy

declare subsamps=(1 2)
declare orders=(0 1)

for o in "${orders[@]}" 
do
	for i in "${subsamps[@]}" 
	do
		echo "addpath(genpath('/scratch/network/astrandb/recipentropy/')); a = ${i}; c=${o}; test(a,c,c);" >> ${CODEDIR}/wrapper_${i}_${i}_${o}.m
		echo "#!/bin/sh" >> ${BASEDIR}/${i}_${o}.sub
		echo "#PBS -l nodes=1:ppn=1,walltime=00:02:00" >> ${BASEDIR}/${i}_${o}.sub
		echo "CODEDIR=/scratch/network/astrandb/recipentropy" >> ${BASEDIR}/${i}_${o}.sub
		echo "DATADIR=/scratch/network/astrandb/recipentropy" >> ${BASEDIR}/${i}_${o}.sub
		echo "cd $CODEDIR" >> ${BASEDIR}/${i}_${o}.sub
		echo "/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -nojvm -r wrapper_${i}_${i}_${o}" >> ${BASEDIR}/${i}_${o}.sub
		echo "exit" >> $BASEDIR/${i}_${o}.sub
		qsub ${BASEDIR}/${i}_${o}.sub
	done
done
		 