package require vtktcl_interactor

# Source the interactor that we will use for the TkRenderWidget
source C:/vtk32/graphics/examplesTcl/TkInteractor.tcl


# prevent the tk window from showing up then start the event loop
wm withdraw .

# Create the toplevel window
toplevel .top
wm title .top {Sphere Puzzle}

# Create some frames
frame .top.f1 
frame .top.f2
pack .top.f1 .top.f2 -side top -expand 1 -fill both

vtkRenderWindow renWin
vtkTkRenderWidget .top.f1.rw -width 400 -height 400 -rw renWin
BindTkRenderWidget .top.f1.rw
pack .top.f1.rw -expand 1 -fill both

# create a rendering window and renderer
vtkRenderer ren1
    renWin AddRenderer ren1

scale .top.f2.s1 -label " Number Of Pieces: " -orient horizontal \
	-length 200 -from 1 -to 32 -variable NUMBER_OF_PIECES 
scale .top.f2.s2 -label " Number Of Subdivision: " -orient horizontal \
	-length 200 -from 0 -to 5 -variable NUMBER_OF_SUBDIVISIONS 

button .top.f2.reset -text "Reset" \
	-command "puzzle Reset; renWin Render"

pack .top.f2.s1 .top.f2.s2 .top.f2.reset -side top -expand 1 -fill both

button .top.f2.b1 -text "Quit" -command {vtkCommand DeleteAllObjects; exit}
pack .top.f2.b1  -expand 1 -fill x

bind .top.f2.s1 <ButtonRelease> { 
   renWin Render
}

bind .top.f2.s2 <ButtonRelease> { 
   renWin Render
}

proc ToggleGhostCells {} {
   renWin Render
}



proc test {slab} {
    for { set i 0} {$i <= 100} {set i [expr $i + 5]} {
	puzzle MoveVertical $slab $i 1
	renWin Render
    }
}


vtkSpherePuzzle puzzle
vtkPolyDataMapper mapper
   mapper SetInput [puzzle GetOutput]
vtkActor actor
    actor SetMapper mapper

vtkSpherePuzzleArrows arrows
vtkPolyDataMapper mapper2
   mapper2 SetInput [arrows GetOutput]
vtkActor actor2
    actor2 SetMapper mapper2




# Add the actors to the renderer, set the background and size
#
ren1 AddActor actor
ren1 AddActor actor2

ren1 SetBackground 0.1 0.2 0.4


bind .top.f1.rw <Motion> {MotionCallback %x %y}
bind .top.f1.rw <KeyPress-m> {ButtonCallback %x %y}

set LastVal -1
proc MotionCallback {x y} {
    global LastVal

    # Compute display point from Tk display point.
    set WindowY [lindex [.top.f1.rw configure -height] 4]
    set y [expr $WindowY - $y]
    set z [ren1 GetZ $x $y]

    ren1 SetDisplayPoint $x $y $z
    ren1 DisplayToWorld
    set pt [ren1 GetWorldPoint]

    #tk_messageBox -message "$pt"
    set x [lindex $pt 0]
    set y [lindex $pt 1]
    set z [lindex $pt 2]

    set val [puzzle SetPoint $x $y $z]
    if {$val != $LastVal} {
	renWin Render
	set LastVal $val
    }
}


proc ButtonCallback {x y} {

    # Compute display point from Tk display point.
    set WindowY [lindex [.top.f1.rw configure -height] 4]
    set y [expr $WindowY - $y]
    set z [ren1 GetZ $x $y]

    ren1 SetDisplayPoint $x $y $z
    ren1 DisplayToWorld
    set pt [ren1 GetWorldPoint]

    #tk_messageBox -message "$pt"
    set x [lindex $pt 0]
    set y [lindex $pt 1]
    set z [lindex $pt 2]

    for { set i 0} {$i <= 100} {set i [expr $i + 5]} {
	puzzle SetPoint $x $y $z
	puzzle MovePoint $i
	renWin Render
    }
}


