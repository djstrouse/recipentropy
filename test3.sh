#!/bin/sh

CODEDIR=/scratch/network/astrandb/recipentropy
BASEDIR=/scratch/network/astrandb/recipentropy

declare -a weights=(10 20 30 40 50 75 95)
declare -a rounds=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)

for i in ${weights[@]}; do
    for j in ${rounds[@]}; do
        echo "weight = ${i}/100; round = ${j}; clust_eval(round,weight);" >> ${CODEDIR}/wrapper_${i}_${j}.m
        cat pbs.info > ${BASEDIR}/${i}_${j}.sub
        echo "/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -nojvm -r wrapper_${i}_${j}" >> ${BASEDIR}/${i}_${j}.sub
        echo "exit" >> ${BASEDIR}/${i}_${j}.sub
    done
done