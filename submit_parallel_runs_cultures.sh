#!/bin/sh

# submit_parallel_runs.sh
# 
#
# Created by Ari Strandburg-Peshkin on 8/23/13.
# Copyright 2013 __MyCompanyName__. All rights reserved.


CODEDIR=/scratch/network/astrandb/recipentropy
BASEDIR=/scratch/network/astrandb/recipentropy

declare subsamps=(16 17)
declare orders=(0)
declare ningreds=(12)
declare cultures=(3)

for c in "${cultures[@]}"
do
	for n in "${ningreds[@]}"
	do
		for o in "${orders[@]}" 
		do
			for i in "${subsamps[@]}" 
			do
				echo "addpath(genpath('/scratch/network/astrandb/recipentropy/')); RecipEntropyCulturesSingle(${c},${n},${o},${i},1)" >> ${CODEDIR}/wrapper_culture${c}_${n}ingreds_o${o}_s${i}.m
				echo "#!/bin/sh" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				echo "#PBS -l nodes=1:ppn=1,walltime=1:00:00" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				echo "CODEDIR=/scratch/network/astrandb/recipentropy" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				echo "DATADIR=/scratch/network/astrandb/recipentropy" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				echo "cd $CODEDIR" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				echo "/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -nojvm -r wrapper_culture${c}_${n}ingreds_o${o}_s${i}" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				echo "exit" >> ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
				qsub ${BASEDIR}/culture${c}_${n}ingreds_o${o}_s${i}.sub
			done
		done
	done
done
		 