# Tcl Script to analyse all the chains/clusters
# for a system of dimers where each dimer
# is a new residue

puts "Hello. Performing Cluster Analysis."

# Sources the ProcChainCount_v4.tcl, box_draw.tcl, and GlobalAverages.tcl
set Sdir /net/uu/nm/cm/bxw109120/Dimers/NestedSampling/CAS

# ProcChainCount programs the function for actual cluster analysis
source $Sdir/ProcChainCount_v4.tcl
source $Sdir/box_draw.tcl
exec mkdir ClusterFiles
set topsel [atomselect top all]
# set a top master list
set master [$topsel get residue]
set nmas [llength $master]

# set up the post master list (pmaster) for value replacing
# this is the one that gets used below
set pmaster $master
# set the init val of outer loop counter
set counter 0
set Number 1
set NumberList [list]

puts " "
puts "Working..."

# start the outer loop over the postmaster list
for {set lk 0} {$lk < 1} {incr counter} {
	# pulls the first value from postmaster list
	# that is not the tag value AA
	set uu [lsearch -not -inline $pmaster AA]

	# pulls the index value of of  first non AA 
	# from postmaster list; returns -1 if
	# there are none
	set uuu [lsearch -not $pmaster AA]
	# loop control, prevent infite looping
	if {$counter > 100} {
		break
	}
	# breaks the loop if there are no 
	# values in pmaster that are not AA
	if {$uuu == -1} {
		break
	}
	if {$counter == 15} {
		puts " "
		puts "Working...."
	}
	# calls ChainCount process which returns the 
	# list of residues in the current cluster
	set clusterlist [vmd_chaincount $uu]
#	set clusterlist [list 23 34 12 56 67 45 34 33 31 30]
	set nclusterlist [llength $clusterlist]
	set BoxCoord [vmd_box_molecule $clusterlist]
	draw text $BoxCoord $uu
	lappend NumberList $uu
#	set nclusterlist $outtotal

	# loops over values in clusterlist and edits pmaster
	# to replace with AA 
	for {set y 0} {$y < $nclusterlist} {incr y} {
		# get residue value from clusterlist
		set clval [lindex $clusterlist $y]
		# pulls the index values for all matching in pmaster
		set rindex [list [lsearch -all $pmaster $clval]]
		#forms rindex data into a new list
		set rrindex [lindex $rindex 0]
		set nrrindex [llength $rrindex]
#		puts "rrindex $rrindex"
#		puts "nrrindex $nrrindex"
		# loops over the rrindex list and replaces 
		# corresponding values of index from pmaster
		for {set vb 0} {$vb < $nrrindex} {incr vb} {
			set rval [lindex $rrindex $vb]
#			puts "rval $rval"
			set pmaster [lreplace $pmaster $rval $rval AA] 
		} 
	}
#	puts "pmaster $pmaster"
}
# this portion collects the data files
# from the ChainCount proc and combines
# them into a single file zzTotalClust.dat
set SNumberlist [lsort -increasing $NumberList]
set datlist [glob ./ClusterFiles/zChain*.dat]
set ndatlist [llength $datlist]
set toutfile [open zzTotalClust.dat w]
for {set po 0} {$po < $ndatlist} {incr po} {
	set current [lindex $datlist $po]
	set xcur [exec cat $current]
	set Numb [lindex $SNumberlist $po]
#	puts $toutfile "Chain/Cluster Number: $Numb"
	puts $toutfile "$xcur"
	puts $toutfile " " 
	
} 
puts $toutfile ""
puts $toutfile "Total chains/clusters: $counter"
puts " "
puts "Total chains/clusters: $counter"
puts " "
source $Sdir/GlobalAverages.tcl
close $toutfile
puts " "
puts "Check output file 'zzTotalClust.dat' for detailed data."
puts "Well, my work is done."
puts "Thank You! Have a nice day!"
