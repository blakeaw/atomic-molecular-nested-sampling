#!/bin/bash

FIRST=0
LAST=0
NUMREP=3
ObsTag="do"
NUMELEMENTS=200

for (( i = FIRST; i <= LAST; i++ )); do
	echo "--Doing trial: $i"
	for (( k = 0; k < NUMREP; k++ )); do
		echo "----Rep: $k"
		sed -e "s/otago/$ObsTag/g" ./pProcFreeOrder_pre.cpp | sed -e "s/itagi/$i/g" | sed -e "s/rrtagrr/$k/g" | sed -e "s/n_el_n/$NUMELEMENTS/g" > prog_temp.cpp
	
			c++ -O2 prog_temp.cpp -o NOP
			./NOP	
	done
	 
		
	

	
done
#rm temp.file
#rm prog.temp
#rm prog_temp.cpp
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




