#!/bin/bash

FIRST=0
LAST=0

#obs="fill"
mkdir SortedConf

for (( i = FIRST; i <= LAST; i++ )); do
	echo "Doing trial: $i"
	
	sed -e "s/itagi/$i/g" /net/uu/nm/cm/bxw109120/Dimers/NestedSampling/pConfigSort_pre.cpp > prog_temp.cpp
	
			c++ -O2 prog_temp.cpp -o SConf
			./SConf	 
		
	

	
done

mv Config*sorted* ./SortedConf
#rm temp.file
#rm prog.temp
rm prog_temp.cpp
echo " " 
echo "Script Running Complete!"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "




