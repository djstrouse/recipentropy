#!/bin/sh

# submit_parallel_runs.sh
# 
#
# Created by Ari Strandburg-Peshkin on 8/23/13.
# Copyright 2013 __MyCompanyName__. All rights reserved.


CODEDIR=/scratch/network/astrandb/recipentropy
BASEDIR=/scratch/network/astrandb/recipentropy

declare subsamps=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)
declare orders=(0 1)
declare ningreds=(4 6 8 10 12 14 16 18)

for n in "${ningreds[@]}"
do
	for o in "${orders[@]}" 
	do
		for i in "${subsamps[@]}" 
		do
			echo "addpath(genpath('/scratch/network/astrandb/recipentropy/')); RecipEntropyMeatsSingle(${n},${o},${i},1)" >> ${CODEDIR}/wrapper_meats_${n}ingreds_o${o}_s${i}.m
			echo "#!/bin/sh" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			echo "#PBS -l nodes=1:ppn=1,walltime=4:00:00:00" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			echo "CODEDIR=/scratch/network/astrandb/recipentropy" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			echo "DATADIR=/scratch/network/astrandb/recipentropy" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			echo "cd $CODEDIR" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			echo "/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -nojvm -r wrapper_meats_${n}ingreds_o${o}_s${i}" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			echo "exit" >> ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
			qsub ${BASEDIR}/meats_${n}ingreds_o${o}_s${i}.sub
		done
	done
done
		 