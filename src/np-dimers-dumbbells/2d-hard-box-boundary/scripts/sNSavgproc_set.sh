#!/bin/bash

CFIRST=3
CLAST=9
INFILE="zNS_combine_Obs_kappaInf.dat"
KAPPA="Inf"


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

	

#rm temp.file
rm prog.temp
rm prog_temp.cpp
echo " " 
echo "Script Running Complete!"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "




