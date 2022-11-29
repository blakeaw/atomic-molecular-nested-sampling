#!/bin/bash

FIRST=0
LAST=0
SKIP=1
NUMREP=4
#obs="fill"
#TYPE="Cv"
#KAPPA="Inf"
TESTLIB="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/TestNumber_Kappa.lib"
TESTdelim="t"

if [ -e "./$cvFILE" ]
then
	echo "$cvFILE found. Cleaning it."
	rm $cvFILE
fi
if [ -e "./$obsFILE" ]
then
	echo "$obsFILE found. Cleaning it."
	rm $obsFILE
fi

for (( i = FIRST; i <= LAST; i=i+SKIP )); do
	echo "--Doing trial: $i"
	KAPPA=`grep -w "${i}${TESTdelim}" $TESTLIB | cut -c7-12`
	cvFILE=`echo "zNS_combine_Cv_kappa${KAPPA}.dat"`
    obsFILE=`echo "zNS_combine_Obs_kappa${KAPPA}.dat"`
	for (( k = 0; k < NUMREP; k++ )); do
		echo "----Rep: $k"
		grep -v "#" zNestOut_int_cube_d2_test_${i}_rep${k}.dat >> $obsFILE
		grep -v "#" zCv_int_cube_d2_test_${i}_rep${k}.dat >> $cvFILE
	done
	echo "Output files are: " 
	echo "$cvFILE"
	echo "$obsFILE"
	

	
done
#rm temp.file
#rm prog.temp
#rm prog_temp.cpp
echo " " 
echo "Script Running Complete!"
#echo "Output file is: $FILE"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "




