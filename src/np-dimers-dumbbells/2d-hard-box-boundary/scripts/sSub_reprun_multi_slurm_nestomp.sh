#!/bin/bash
#Starting Test Number
FIRST=0
#Ending Test Number
LAST=0
#Number of replicates to run
NUMREP=0
##Number of cpus to split sampling look over
NUMCPU=4
#Pre NS program file
Prefile="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/pNest_dimers_parallel_replicates_nestedomp_2d_mwm_v1.0_pre.cpp"
Wdir=`pwd`
#TestNumber-kappa Lib file
TESTLIB="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/TestNumber_Kappa.lib"
#delimiting character to pull values from Lib file
TESTdelim="t"
#slurm submit script
#Slurmsubmit="/net/uu/nm/cm/bxw109120/Dimers/NestedSampling/2dimensional/scripts/sSlum_submit_ompjob_b.sh"
#first portion of output binary name
cpptemp="bin_job"

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
	sed -e "s/xxntestnxx/$i/g" $Prefile | sed -e "s/n_rep_n/$NUMREP/g" | sed -e "s/k_val_k/$kappa/g" | sed -e "s/n_cpu_n/$NUMCPU/g" > prog_temp.cpp
	#compile
	c++ -fopenmp -O3 prog_temp.cpp -o $binname
	time ./$binname > out.${binname}_rep${NUMREP} &
	#---Slurm--submit to queue
	#submit to as batch run calling the slurm submit script {Slurmsubmit}
#	sbatch -s -n$NUMREP -c$NUMCPU $Slurmsubmit $NUMREP $NUMCPU $binname
##	srun -s -n ${NUMREP} -c ${NUMCPU} -l ./$binname > out.${binname} &
##	$Slurmsubmit $NUMREP $binname
#	rm prog_temp.cpp

	
	
done

#rm $cpptemp
#rm prog_temp.cpp
echo " " 
echo "Script Running Complete!"
echo "Have a nice day!"
echo " " 
echo "...End of Line..."
echo " "

