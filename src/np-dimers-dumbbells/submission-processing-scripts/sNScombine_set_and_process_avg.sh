#!/bin/bash

FIRST=14
LAST=17
SKIP=1
NUMREP=4
#obs="fill"
#TYPE="Cv"
#KAPPA="Inf"
TESTLIB="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/TestNumber_Kappa.lib"
TESTdelim="t"

### Combination



for (( i = FIRST; i <= LAST; i=i+SKIP )); do
	echo "--Doing trial: $i"
	KAPPA=`grep -w "${i}${TESTdelim}" $TESTLIB | cut -c7-12`
	cvFILE=`echo "zNS_combine_Cv_kappa${KAPPA}.dat"`
    obsFILE=`echo "zNS_combine_Obs_kappa${KAPPA}.dat"`
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
	for (( k = 0; k < NUMREP; k++ )); do
		echo "----Rep: $k"
		grep -v "#" zNestOut_int_cube_d8_test_${i}_rep${k}.dat >> $obsFILE
		grep -v "#" zCv_int_cube_d8_test_${i}_rep${k}.dat >> $cvFILE
	done
	echo "Combination complete."
	echo "Output files are: " 
	echo "$cvFILE"
	echo "$obsFILE"
	
### Processing
	echo "Beginning Processing..."
####Cv
	echo "Processing average Cv values.."
	INFILE=`echo "$cvFILE"`
	sed -e "s/i_file_i/$INFILE/g" pCvAvg_pre.cpp | sed -e "s/k_val_k/$KAPPA/g" > prog_temp.cpp 
	c++ -O2 prog_temp.cpp -o CVAVG
	./CVAVG
	rm prog_temp.cpp	
 
####Obs
	echo "Processing average observables..."
	CFIRST=2
	CLAST=9
	INFILE=`echo "$obsFILE"`
	
	sed -e "s/i_file_i/$INFILE/g" pObsAvg_pre.cpp | sed -e "s/k_val_k/$KAPPA/g" > prog.temp 

	for (( j = CFIRST; j <= CLAST; j++ )); do
		echo "Doing column number: $j"
			if [ $j = 2 ]; then
				obs="en"
			fi
			if [ $j = 3 ]; then
				obs="es"
			fi
			if [ $j = 4 ]; then
				obs="cv"
			fi
			if [ $j = 5 ]; then
				obs="rg"
			fi
			if [ $j = 6 ]; then
				obs="gd"
			fi
			if [ $j = 7 ]; then
				obs="zo"
			fi
			if [ $j = 8 ]; then
				obs="do"
			fi
			if [ $j = 9 ]; then
				obs="zd"
			fi
		echo "which has obs tag: $obs"
			sed -e "s/otago/$obs/g" prog.temp | sed -e "s/o_col_o/$j/g" > prog_temp.cpp
			c++ -O2 prog_temp.cpp -o OAVG
			./OAVG	 
		#	grep -v "0.000000000000000" zNSval_d2_${obs}_test_${i}.dat > temp.file
		#	cp temp.file zNSval_d2_${obs}_test_${i}.dat
	done 

	


	rm prog.temp
	rm prog_temp.cpp	
	echo "Processing for trial $i with kappa $KAPPA complete."
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




