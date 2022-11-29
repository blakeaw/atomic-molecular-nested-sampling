#!/bin/bash
#Starting Test Number
FIRST=0
#Ending Test Number
LAST=13
#Number of replicates to run
NUMREP=4
##Number of cpus to split sampling look over
#NUMCPU=4
#Pre NS program file
Prefile="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/pNest_dimers_parallel_replicates_nestedomp_pre.cpp"
Wdir=`pwd`
#TestNumber-kappa Lib file
TESTLIB="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/TestNumber_Kappa.lib"
#delimiting character to pull values from Lib file
TESTdelim="t"
#slurm submit script
Slurmsubmit="./scripts/sSlum_submit_ompjob_b.sh"
#first portion of output binary name
cpptemp="bin_job"
##SBATCH --nodes=5
#SBATCH --ntasks-per-node=4
# Loop over test numbers
for (( i = FIRST; i <= LAST; i++ )); do
	echo "Doing test: $i"
	#get the kappa value from the lib file
	kappa=`grep -w "${i}${TESTdelim}" $TESTLIB | cut -c7-12`
#	kappa=0.1
	
	echo "with kappa $kappa"
	echo " and numreps $NUMREP"
	binname=`echo "${cpptemp}_${i}"`
	#get pre-program and replace test number {i}, replicate number {NUMREP}, and kappa {kappa}
	sed -e "s/xxntestnxx/$i/g" $Prefile | sed -e "s/n_rep_n/$NUMREP/g" | sed -e "s/k_val_k/$kappa/g" > prog_temp.cpp
	#compile
	c++ -O2 -fopenmp -ffast-math -fcx-fortran-rules prog_temp.cpp -o $binname
	#---Slurm--submit to queue
	##SBATCH --ntasks-per-node=4
	#submit to as batch run calling the slurm submit script {Slurmsubmit}
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

