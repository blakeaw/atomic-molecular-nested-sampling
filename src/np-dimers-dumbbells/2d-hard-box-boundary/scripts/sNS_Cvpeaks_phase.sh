#!/bin/bash

FIRST=0
LAST=20
SKIP=1
NUMREP=3
NUMDIM=8
TESTLIB="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/TestNumber_Rscale.lib"
TESTdelim="t"
prc="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/2dimensional/d8/bCVmax"
### Combination



for (( i = FIRST; i <= LAST; i=i+SKIP )); do
	echo "--Doing trial: $i"
	RSCALE=`grep -w "${i}${TESTdelim}" $TESTLIB | cut -c7-12`
	cvFILE=`echo "aNest_Cv_avg_rscale${RSCALE}.dat"`
    
	if [ -e "./$cvFILE" ]
	then
		echo "$cvFILE found."
		$prc $cvFILE > temp.file
		cat temp.file
		nmax=`grep "complete" temp.file | cut -c17`
		echo "has $nmax maxima"
		kappa=`grep "kappa" oMMC_cube_d${NUMDIM}_test_${i}_rep0.log | cut -c8-15`
		echo "with kappa $kappa"
		if [[ $nmax -eq 1 ]]; then
			value=`head -n1 temp.file | cut -c24-28`
			echo "$value $kappa" >> Phase_LS.dat
			 
		fi
		if [[ $nmax -eq 2 ]]; then
			value1=`head -n1 temp.file | cut -c24-28`
			value2=`head -n2 temp.file | tail -n1 | cut -c24-28`
			echo "$value1 $kappa" >> Phase_S1S2.dat
			echo "$value2 $kappa" >> Phase_LS1.dat
		fi
		echo "Processing for trial $i with rscale $RSCALE complete."
	fi
	
	
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




