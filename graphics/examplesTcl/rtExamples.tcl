catch {load vtktcl}
#
# This is a regression test script for VTK.
#

# first find all the examples. they can be defined on command line or in
# current directory

if { $argv != ""} {
  set files $argv
  } else {
  set files [lsort [glob {[A-z]*.tcl}]]
}

# remove files that are not suitable for regression tests or simply don't work right now
set noTest {
   rt.tcl rtAll.tcl rib.tcl TkInteractor.tcl TkRenderWidget.tcl RenderWidget.tcl 
   rtExamples.tcl polyViewer.tcl KeyFrame.tcl cameraKey.tcl timing.tcl 
   Decimate.tcl assembly2.tcl connPineRoot.tcl deciHawa.tcl 
   deciPineRoot.tcl deleted.tcl mcTest.tcl viewMCubesFile.tcl vol.tcl
   volTkInteractor.tcl spikeColor.tcl tkwin.tcl 3dsToRIB.tcl backdrop.tcl

   Close.tcl ContinuousClose.tcl HighPassComparison.tcl 
   HybridMedianComparison.tcl IdealHighPass.tcl ShotNoiseInclude.tcl 
   Spectrum.tcl TestContinuousDilate3D.tcl TestContinuousErode3D.tcl 
   TestDilateErode3D.tcl TestFeatureAnd.tcl TestHistogram.tcl 
   TestHistogramEqualization.tcl TestHybridMedian2D.tcl TestIdealHighPass.tcl 
   TestIdealLowPass.tcl TestLogarithmicScale.tcl TestMIPFilter.tcl
   TestOpenClose3D.tcl TestRange3D.tcl TestSkeleton2D.tcl
   TestSobel2D.tcl TestSobel3D.tcl TestSubsample3D.tcl TestVariance3D.tcl
   TestWriter.tcl Timing.tcl VTKSpectrum.tcl WindowLevelInterface.tcl
   vtkImageInclude.tcl ContinuousClose2D.tcl TkViewer2.tcl Pad.tcl
}

for {set i 0} {$i < [llength $noTest]} {incr i} {
   if {[set pos [lsearch $files [lindex $noTest $i]]] != -1} {
      set files [lreplace $files $pos $pos ]
   }
}

# now do the tests
foreach afile $files {
   #
   # only tcl scripts with valid/ images are tested
   if {[catch {set channel [open "valid/$afile.ppm"]}] != 0 } {
      puts "WARNING: There is no valid image for $afile"
      continue
   } else {
      close $channel
   }
   
   vtkMath rtExMath
   rtExMath RandomSeed 6
   
   puts -nonewline "$afile took "
   puts -nonewline "[expr [lindex [time {source $afile} 1] 0] / 1000000.0] seconds "
   
   vtkWindowToImageFilter w2if
   
   # look for a renderWindow ImageWindow or ImageViewer
   # first check for some common names
   if {[info commands renWin] == "renWin"} {
      w2if SetInput renWin
      set threshold 10
   } else {
      set threshold 0
      if {[info commands viewer] == "viewer"} {
	 w2if SetInput [viewer GetImageWindow]
	 viewer Render
      } else {
	 if {[info commands imgWin] == "imgWin"} {
	    w2if SetInput imgWin
	    imgWin Render
	 }
      }
   }

   # run the event loop quickly to map any tkwidget windows
   wm withdraw .
   update

   if {[catch {set channel [open "valid/$afile.ppm"]}] != 0 } {
      puts "\nWARNING: Creating a valid image for $afile"
      vtkPNMWriter rtpnmw
      rtpnmw SetInput [w2if GetOutput]
      rtpnmw SetFileName "valid/$afile.ppm"
      rtpnmw Write
      vtkCommand DeleteAllObjects
      catch {destroy .top}
      continue
   }
   close $channel
   
   vtkPNMReader rtpnm
   rtpnm SetFileName "valid/$afile.ppm"
   
   vtkImageDifference imgDiff
   
   if {$threshold == 0} {imgDiff AllowShiftOff; imgDiff SetThreshold 0}
   imgDiff SetInput [w2if GetOutput]
   imgDiff SetImage [rtpnm GetOutput]
   imgDiff Update

   # a test has to be off by at least threshold pixels for us to care   
   if {[imgDiff GetThresholdedError] <= $threshold} {
       puts "and Passed"
   } else {
       puts "but failed with an error of [imgDiff GetThresholdedError]"
         vtkPNMWriter rtpnmw
         rtpnmw SetInput [imgDiff GetOutput]
         rtpnmw SetFileName "$afile.error.ppm"
         rtpnmw Write
         vtkPNMWriter rtpnmw2
         rtpnmw2 SetInput [w2if GetOutput]
         rtpnmw2 SetFileName "$afile.test.ppm"
         rtpnmw2 Write
   }
   
   vtkCommand DeleteAllObjects
   catch {destroy .top}
}

exit