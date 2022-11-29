#!/bin/bash

FIRST=0
LAST=9

#obs="fill"

for (( i = FIRST; i <= LAST; i++ )); do
	echo "Doing trial: $i"
	#grep -v "nan" zNSobs_int_cube_d2_test_$i.dat > temp.file
	cat zNSobs_int_cube_d2_test_$i.dat | sed -e '/nan/,$d' > temp.file
	npt=`wc temp.file | cut -c2-5 `
	echo "npt: $npt"
	sed -e "s/n_pt_n/$npt/g" ./pNSvalues_pre.cpp | sed -e "s/itagi/$i/g" > prog.temp
	for (( j = 3; j < 10; j++ )); do
		echo "Doing column number: $j"
			if [ $j = 3 ]; then
				obs="en"
			fi
			if [ $j = 4 ]; then
				obs="es"
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
			sed -e "s/o_curr_o/$j/g" prog.temp | sed -e "s/otago/$obs/g" > prog_temp.cpp
			c++ -O2 prog_temp.cpp -o NSV
			./NSV	 
		#	grep -v "0.000000000000000" zNSval_d2_${obs}_test_${i}.dat > temp.file
			cp temp.file zNSval_d2_${obs}_test_${i}.dat
	done 

	
done
rm temp.file
rm prog.temp
rm prog_temp.cpp
echo " " 
echo "Script Running Complete!"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "




