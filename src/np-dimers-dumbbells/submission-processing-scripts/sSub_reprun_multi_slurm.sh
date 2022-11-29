#!/bin/bash

FIRST=0
LAST=9
NUMREP=4
Prefile="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/pNest_dimers_parallel_replicates_pre.cpp"
Wdir=`pwd`
TESTLIB="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/TestNumber_Kappa.lib"
TESTdelim="t"
Slurmsubmit="./scripts/sSlum_submit_ompjob_b.sh"
cpptemp="bin_job"
##SBATCH --nodes=5
#SBATCH --ntasks-per-node=4

for (( i = FIRST; i <= LAST; i++ )); do
	echo "Doing test: $i"
	kappa=`grep -w "${i}${TESTdelim}" $TESTLIB | cut -c7-12`
#	kappa=0.1
	
	echo "with kappa $kappa"
	echo " and numreps $NUMREP"
	binname=`echo "${cpptemp}_${i}"`
	sed -e "s/xxntestnxx/$i/g" $Prefile | sed -e "s/n_rep_n/$NUMREP/g" | sed -e "s/k_val_k/$kappa/g" > prog_temp.cpp
	c++ -O2 -fopenmp -ffast-math -fcx-fortran-rules prog_temp.cpp -o $binname
	#---Slurm--submit to queue
	##SBATCH --ntasks-per-node=4
	sbatch $Slurmsubmit $NUMREP $binname
##	$Slurmsubmit $NUMREP $binname
	rm prog_temp.cpp

	
	
done

#rm $cpptemp
#rm prog_temp.cpp
echo " " 
echo "Script Running Complete!"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "

