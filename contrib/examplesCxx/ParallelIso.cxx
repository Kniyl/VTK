// This program test the ports by setting up a simple pipeline.

#include "vtkImageReader.h"
#include "vtkSynchronizedTemplates3D.h"
#include "vtkAppendPolyData.h"
#include "vtkOutputPort.h"
#include "vtkInputPort.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkElevationFilter.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTimerLog.h"
#include "vtkMath.h"




#define ISO_START 3500.0
#define ISO_STEP  -500.0
#define ISO_NUM   6




// callback to test streaming / ports by seeing what extents are being read in.
void reader_start_callback(void *arg)
{
  vtkImageReader *reader = (vtkImageReader*)(arg);
  int *e;
  
  e = reader->GetOutput()->GetUpdateExtent();
  
  cerr << "Reading: " << e[0] << ", " << e[1] << ", " << e[2] << ", " 
       << e[3] << ", " << e[4] << ", " << e[5] << endl; 
}



// call back to set the iso surface value.
void set_iso_val_rmi(void *arg, int id)
{ 
  float val;

  vtkSynchronizedTemplates3D *iso;
  iso = (vtkSynchronizedTemplates3D *)arg;
  val = iso->GetValue(0);
  iso->SetValue(0, val + ISO_STEP);
}

// call back to exit program
// This should really be embedded in the controller.
void exit_callback(void *arg, int id)
{ 
  // clean up controller ?
  exit(0);
}



VTK_THREAD_RETURN_TYPE process( void *vtkNotUsed(arg) )
{
  vtkMultiProcessController *controller;
  vtkImageReader *reader;
  vtkSynchronizedTemplates3D *iso;
  vtkElevationFilter *elev;
  int myid, numProcs;
  float val;
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);
  myid = controller->GetLocalProcessId();
  numProcs = controller->GetNumberOfProcesses();
    
  reader = vtkImageReader::New();
  reader->SetDataByteOrderToLittleEndian();
  reader->SetDataExtent(0, 127, 0, 127, 1, 93);
  reader->SetFilePrefix("../../../vtkdata/headsq/half");
  reader->SetDataSpacing(1.6, 1.6, 1.5);
  
  iso = vtkSynchronizedTemplates3D::New();
  iso->SetInput(reader->GetOutput());
  iso->SetValue(0, ISO_START);
  iso->ComputeScalarsOff();
  iso->ComputeGradientsOff();
  // This should be automatically determined by controller.
  iso->SetNumberOfThreads(1);
  
  // Compute a different color for each process.
  elev = vtkElevationFilter::New();
  elev->SetInput(iso->GetOutput());
  vtkMath::RandomSeed(myid * 100);
  val = vtkMath::Random();
  elev->SetScalarRange(val, val+0.001);

  if (myid != 0)
    {
    // Remote process! Send data throug port.
    vtkOutputPort *upPort = vtkOutputPort::New();
    
    // last, set up a RMI call back to change the iso surface value.
    controller->AddRMI(set_iso_val_rmi, (void *)iso, 300);
    controller->AddRMI(exit_callback, (void *)iso, 666);
  
    upPort->SetInput(elev->GetPolyDataOutput());
    // the different process ids differentiate between sources.
    upPort->SetTag(999);
    // wait for the call back to execute.
    upPort->WaitForUpdate();
    
    // last call never returns, but ...
    upPort->Delete();
    }
  else
    {
    int i, j;
    vtkAppendPolyData *app = vtkAppendPolyData::New();
    vtkInputPort *downPort;
    vtkRenderer *ren = vtkRenderer::New();
    vtkRenderWindow *renWindow = vtkRenderWindow::New();
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    vtkActor *actor = vtkActor::New();
    vtkTimerLog *timer = vtkTimerLog::New();
    vtkCamera *cam = vtkCamera::New();
    
    // This is the main thread: Collect the data and render it.
    app->AddInput(elev->GetPolyDataOutput());
    // ###################### important ####################
    // # this tells the append filter to request pieces from
    // # each of its inputs.
    app->ParallelStreamingOn();
    
    for (i = 1; i < numProcs; ++i)
      {
      downPort = vtkInputPort::New();
      downPort->SetRemoteProcessId(i);
      downPort->SetTag(999);
      app->AddInput(downPort->GetPolyDataOutput());
      // referenced by app ...
      downPort->Delete();
      downPort = NULL;
      }
    
    //putenv("DISPLAY=:0.0");
    
    renWindow->AddRenderer(ren);
    iren->SetRenderWindow(renWindow);
    ren->SetBackground(0.9, 0.9, 0.9);
    renWindow->SetSize( 400, 400);
  
    mapper->SetInput(app->GetOutput());
    actor->SetMapper(mapper);
  
    // assign our actor to the renderer
    ren->AddActor(actor);
  
    cam->SetFocalPoint(100, 100, 65);
    cam->SetPosition(100, 450, 65);
    cam->SetViewUp(0, 0, -1);
    cam->SetViewAngle(30);
    cam->ComputeViewPlaneNormal();
    // this was causing an update.
    //ren->ResetCameraClippingRange();
    //{
    //double *range = ren->GetActiveCamera()->GetClippingRange();
    //cerr << range[0] << ", " << range[1] << endl;
    //}
    cam->SetClippingRange(177.0, 536.0);
    ren->SetActiveCamera(cam);
    
    // loop through some iso surface values.
    for (j = 0; j < ISO_NUM; ++j)
      {
      // set the local value
      set_iso_val_rmi((void*)iso, 0);
      for (i = 1; i < numProcs; ++i)
	{
	// trigger the RMI to change the iso surface value.
	controller->TriggerRMI(i, 300);      
	}
      
      timer->StartTimer();
      app->Update();
      timer->StopTimer();
      cerr << "Update " << val << " took " << timer->GetElapsedTime() 
	   << " seconds\n";
      
      // now render the results
      renWindow->Render();
      //sprintf(filename, "iso%d.ppm", (int)(val));
      //renWindow->SetFileName(filename);
      //renWindow->SaveImageAsPPM();
      val += 400.0;
      }
    
    // just exit
    //for (i = 1; i < numProcs; ++i)
    //  {
    //  // trigger the RMI to exit
    //  controller->TriggerRMI(i, 666);      
    //  }
    //exit(0);
    
    //  Begin mouse interaction
    iren->Start();
    
    // Clean up
    app->Delete();
    ren->Delete();
    renWindow->Delete();
    iren->Delete();
    mapper->Delete();
    actor->Delete();
    }
  
  // clean up objects in all processes.
  controller->UnRegister(NULL);
  reader->Delete();
  iso->Delete();
  elev->Delete();

  return VTK_THREAD_RETURN_VALUE;
}


void main( int argc, char *argv[] )
{
  vtkMultiProcessController *controller;
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);

  controller->Initialize(argc, argv);
  // Needed for threaded controller.
  // controller->SetNumberOfProcesses(2);
  controller->SetSingleMethod(process, NULL);
  controller->SingleMethodExecute();

  controller->UnRegister(NULL);
}





