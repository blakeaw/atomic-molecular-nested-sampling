#set R [list 87 86]
proc vmd_box_molecule {R} {
      # get the min and max values for each of the directions
      set sel [atomselect top "residue $R"]

      set coords [lsort -real [$sel get x]]
      set minx [expr [lindex $coords 0] - 10]
      set maxx [expr [lindex [lsort -real -decreasing $coords] 0] + 10]

      set coords [lsort -real [$sel get y]]
      set miny [expr [lindex $coords 0] - 10]
      set maxy [expr [lindex [lsort -real -decreasing $coords] 0] + 10]

      set coords [lsort -real [$sel get z]]
      set minz [expr [lindex $coords 0] - 10]
      set maxz [expr [lindex [lsort -real -decreasing $coords] 0] + 10]

      # and draw the lines
      draw materials off
      draw color yellow
      draw line "$minx $miny $minz" "$maxx $miny $minz"
      draw line "$minx $miny $minz" "$minx $maxy $minz"
      draw line "$minx $miny $minz" "$minx $miny $maxz"

      draw line "$maxx $miny $minz" "$maxx $maxy $minz"
      draw line "$maxx $miny $minz" "$maxx $miny $maxz"

      draw line "$minx $maxy $minz" "$maxx $maxy $minz"
      draw line "$minx $maxy $minz" "$minx $maxy $maxz"

      draw line "$minx $miny $maxz" "$maxx $miny $maxz"
      draw line "$minx $miny $maxz" "$minx $maxy $maxz"

      draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
      draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
      draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
      return "$minx $miny $maxz"
#}
