#Tcl script for addition to ClusterAnalysis
#Calculates and ouputs the global average values for 
# the total cluster size(# of dimers), the local residue particle density,
# and the chain/cluster radius of gyration over all the chains/clusters

set datlist [glob ./ClusterFiles/zChain*.dat]
set ndatlist [llength $datlist]
#puts "Initializing..."
#Initialize the arrays
for {set po 0} {$po < $ndatlist} {incr po} {
	set Cbin(C,$po) 0.00
	set Abin(A,$po) 0.00
	set Tbin(T,$po) 0.00	
} 

#Gets the values from the Cluster files and stores them in an array
#puts "Getting Values...."
for {set po 0} {$po < $ndatlist} {incr po} {
#	puts "Value Loop po $po"
	set current [lindex $datlist $po]
	set Cbin(C,$po) [exec egrep "CST" $current | cut -c5-10]
	set Abin(A,$po) [exec egrep "ALRPD" $current | cut -c7-15]
	set Tbin(T,$po) [exec egrep "TACROG" $current | cut -c7-20]
	
} 
#initialize the sums
set SumC 0.00
set SumA 0.00
set SumT 0.00
# Run the sum over each array
#puts "Summing the Values..."
for {set po 0} {$po < $ndatlist} {incr po} {
	set SumC [expr $SumC + $Cbin(C,$po)]
	set SumA [expr $SumA + $Abin(A,$po)]
	set SumT [expr $SumT + $Tbin(T,$po)]
}
# get the average value	
#puts "Averaging...."
set avgC [expr $SumC / $ndatlist]
set avgA [expr $SumA / $ndatlist]
set avgT [expr $SumT / $ndatlist]
#output
puts $toutfile "Global Averages:"
puts $toutfile "ACST: $avgC"
puts $toutfile "AASRPD: $avgA"
puts $toutfile "ATACROG: $avgT"
puts "Global Averages:"
puts "ACST: $avgC"
puts "AASRPD: $avgA"
puts "ATACROG: $avgT"

