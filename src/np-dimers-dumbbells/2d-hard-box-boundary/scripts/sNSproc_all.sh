#!/bin/bash

FIRST=1
LAST=20
NUMREP=3
NUMDIM=8
#obs="fill"
prc="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/2dimensional/ProcProg"
for (( i = FIRST; i <= LAST; i++ )); do
	echo "--Doing trial: $i"
	for (( k = 0; k < NUMREP; k++ )); do
		echo "----Rep: $k"
		sed -e "s/n_pt_n/$npt/g" $prc/pNestProc_pre.cpp | sed -e "s/itagi/$i/g" | sed -e "s/rrtagrr/$k/g" > prog_temp1.cpp
		sed -e "s/d_dim_d/$NUMDIM/g" prog_temp1.cpp > prog_temp.cpp
			c++ -O2 prog_temp.cpp -o NOP
			./NOP	
	done
	 
		
	

	
done
#rm temp.file
#rm prog.temp
rm prog_temp1.cpp
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




