# Script to measure the number of dimers in a chain or cluster
# also outputs the time averaged local particle density around
# each of the residues(dimers), derived from CountClust.tcl
## *** User Defined Variables: sr
proc vmd_chaincount {sr} {
#puts "Here we go!"
# Define the residue to start analysis from
# can be any residue in the chain/cluster
#set sr  22

#Define the output file
#set outfile [open zChainSR$sr.dat w]
set outfile [open ./ClusterFiles/zChain$sr.dat w]
# Start the neighbor list with sr initially
set neilist [list $sr]
set tag 0

#puts "Building neighbor list..."

set n [llength $neilist]
#puts " "
# loop over elements in the neighbor list
# to extract all neighbors and neighbor of neighbors.
# neilist is self updated in the loop to include
# all new residues that are within the cutoff
# distance of any other residue in the chain/cluster.
for {set m 0} {$m < $n} {incr m} {

	set aval [lindex $neilist $m]
	# extract the neighbors of the residue lindex m from neilist
	set csel [atomselect top "within 4 of residue $aval"]
	#add the neigbor of neighbors to a new list
	set creslist [$csel get residue]
	set cnreslist [llength $creslist]
	set test $neilist
	set ntest [llength $test]
	# Loop over the neighbor of neighbor list(creslist) elements
	for {set k 0} {$k < $cnreslist} {incr k} {
		set test $neilist
		set ntest [llength $test]
		#extract element k from creslist 
		set bval [lindex $creslist $k]
		# tag is used later for testing
		set tag 0
		# loop over values of list test(updated neilist)
		# and compare the current creslist value
		for {set p 0} {$p < $ntest} {incr p} {
			
			set nval [lindex $test $p]
			# testing condition, compares current creslist
			# value to all the values in the neilist
			# and sets adds to tag if the creslist value
			# is already in the neilist
			if {$bval == $sr || $bval == $nval} {
				set tag [expr $tag + 1]
			
			} 
		}
		# after comparing the current element of creslist
		# to the neilist, this if updates neilist to include
		# the new element if tag is still zero and increments
		# n to increase the neilist length
		if {$tag == 0} {

			lappend neilist $bval
			incr n
		}
	
	# a loop control, to prevent infinite looping
	# should be larger than the number of dimers in
	# the chain or cluster
	}
	if {$n > 50} {
		break
	}
}	

# Start the initial output
set nneilist [llength $neilist]
puts $outfile "Starting Residue: $sr"
puts $outfile "Residues in Chain/Cluster:"
# loops over elements of neilist and extracts the 
# residues, adding them to a new line of the outfile
for {set j 0} {$j < [expr $nneilist -0]} {incr j} {
	set pval [lindex $neilist $j]
	puts $outfile "$pval"
}
# get the total number of dimers in the chain/cluster
set total [expr $nneilist]
puts $outfile "Chain/Cluster Size is $total dimers"
puts $outfile "CST $total"
#puts "Chain Length/Cluster Size is $total dimers"
#puts " "
# Start of the CountClust style analysis
#puts "Starting Local Residue Cluster Analysis..."
#puts " "
puts $outfile "Local Residue Cluster Analysis:"
set nf [molinfo top get numframes]
#set starting frame for analysis
#currently set to halfway of total frames
# to allow some initial formation time
set fstart [expr $nf / 2]
# Initialize the Avg Array
for {set k 0} {$k< $total} {incr k} {
	set fbin($k,Avg) 0.00
	
}
# Initialize Value array for summation over frames
for {set k $fstart} {$k <= $nf} {incr k} {
		set fbin($k,Res) 0.00
		set fbin($k,Gyr) 0.00
}
#loop over the list, neilist, generated in previous part
# and measure the number of particles within the cutoff
# distance of each residue
for {set i 0} {$i < $total} {incr i} {
	set tval [lindex $neilist $i]
	
	set sel [atomselect top "within 4 of residue $tval"]
	set selall [atomselect top "residue $neilist"]
	# loop over frames
	for {set j $fstart} {$j<= $nf} {incr j} {
		$sel frame $j
		$sel update
		$selall frame $j
		$selall update
		set rnum [$sel num]
		set brnum [expr $rnum - 2]
		set fbin($j,Res) $brnum
		if {$i == 0} {
			set RG [measure rgyr $selall]
			set fbin($j,Gyr) $RG
		}
 
	}
	set sumClust 0.00
	

	# loop over frames and sum values in array fbin("",Res)
	for {set p $fstart} {$p <= $nf} {incr p} {
		set sumClust [expr $sumClust + $fbin($p,Res)]
		
	}
	# avgClust is the time averaged local particle number
	set avgClust [expr $sumClust/$fstart]
	
	puts $outfile "$tval $avgClust"
	set fbin($i,Avg) $avgClust
		
	
}
# start calculating the total average value over all residues
set sumAvg 0.0
for {set u 0} {$u < $total} { incr u} {
	set sumAvg [expr $sumAvg + $fbin($u,Avg)]
}
set SumAvg [expr $sumAvg/$total]
# output to screen and outfile
puts $outfile "Average Local Residue Particle Density: $SumAvg"
puts $outfile "ALRPD $SumAvg"
#puts "Average Local Residue Cluster Size: $SumAvg"
# closes the output file
set sumGyr 0.00
for {set p $fstart} {$p <= $nf} {incr p} {
	set sumGyr [expr $sumGyr + $fbin($p,Gyr)]
}
set avgGyr [expr $sumGyr/$fstart]
puts $outfile "Time Averaged Chain/Cluster Radius of Gyration: $avgGyr"
puts $outfile "TACROG $avgGyr"

close $outfile

return $neilist

#puts " "
#puts "Script Complete"
#puts "Have a nice day!"
}
