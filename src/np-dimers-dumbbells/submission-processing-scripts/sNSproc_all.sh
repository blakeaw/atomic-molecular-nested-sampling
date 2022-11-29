#!/bin/bash

FIRST=0
LAST=17
NUMREP=3
#obs="fill"

for (( i = FIRST; i <= LAST; i++ )); do
	echo "--Doing trial: $i"
	for (( k = 0; k < NUMREP; k++ )); do
		echo "----Rep: $k"
		sed -e "s/n_pt_n/$npt/g" ./pNestProc_pre.cpp | sed -e "s/itagi/$i/g" | sed -e "s/rrtagrr/$k/g" > prog_temp.cpp
	
			c++ -O2 prog_temp.cpp -o NOP
			./NOP	
	done
	 
		
	

	
done
#rm temp.file
#rm prog.temp
rm prog_temp.cpp
# test directory
#if [ ! -d "ProcNestedOut" ]; then
  # Control will enter here if DIRECTORY doesn't exist.
#	mkdir ProcNestedOut
#fi

echo " " 
echo "Script Running Complete!"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "




