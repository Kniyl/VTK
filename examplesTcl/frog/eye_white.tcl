source frog.tcl

set NAME eye_white
set TISSUE 5
set START_SLICE 1
set END_SLICE 37
set ZMAX [expr $END_SLICE - $START_SLICE]
set VOI "389 433 183 282 0 $ZMAX"

source segmented8.tcl
