// This program test the ports by setting up a simple pipeline.

#include "vtkOutputPort.h"
#include "vtkInputPort.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkElevationFilter.h"
#include "vtkRenderWindowInteractor.h"

VTK_THREAD_RETURN_TYPE process_a( void *vtkNotUsed(arg) )
{
  vtkConeSource *cone = vtkConeSource::New();
  vtkElevationFilter *elev = vtkElevationFilter::New();
  vtkOutputPort *upStreamPort = vtkOutputPort::New();
  
  // Set up the pipeline source.
  cone->SetResolution(8);
  elev->SetInput(cone->GetOutput());
  upStreamPort->SetInput(elev->GetPolyDataOutput());
  upStreamPort->SetTag(999);
  
  // wait for the call back to execute.
  upStreamPort->WaitForUpdate();
  
  cone->Delete();
  elev->Delete();
  upStreamPort->Delete();

  return VTK_THREAD_RETURN_VALUE;
}


VTK_THREAD_RETURN_TYPE process_b( void *arg )
{
  vtkMultiProcessController *controller;
  int myid, otherid;
  char *save_filename = (char*)arg;
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);
  myid = controller->GetLocalProcessId();
  if (myid == 0)
    {
    otherid = 1;
    }
  else
    {
    otherid = 0;
    }

  vtkInputPort *downStreamPort = vtkInputPort::New();
  downStreamPort->SetRemoteProcessId(otherid);
  downStreamPort->SetTag(999);
  downStreamPort->GetPolyDataOutput()->SetUpdateExtent(0, 4);
  downStreamPort->Update();  

  vtkPolyDataMapper *coneMapper = vtkPolyDataMapper::New();
  coneMapper->SetInput(downStreamPort->GetPolyDataOutput());

  vtkActor *coneActor = vtkActor::New();
  coneActor->SetMapper(coneMapper);
  
  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor(coneActor);
  ren->SetBackground(0.1, 0.3, 0.5);
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  renWin->SetSize( 300, 300 );

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);  
  
  // draw the resulting scene
  renWin->Render();

  // save for the regression test
  if (save_filename != NULL && save_filename[0] != '\0')
    {
    renWin->SetFileName( save_filename );
    renWin->SaveImageAsPPM();
    // Tell the other process to stop waiting.
    controller->TriggerRMI(otherid, VTK_BREAK_RMI_TAG);
    }
  else
    {
    //  Begin mouse interaction
    iren->Start();
    }
  
  // Clean up
  ren->Delete();
  renWin->Delete();
  iren->Delete();
  downStreamPort->Delete();
  coneMapper->Delete();
  coneActor->Delete();

  return VTK_THREAD_RETURN_VALUE;
}


void main( int argc, char *argv[] )
{
  vtkMultiProcessController *controller;
  char save_filename[100];

  save_filename[0] = '\0';
  if( (argc >= 2) && (strcmp("-S", argv[argc-1]) == 0) )
    {
    sprintf( save_filename, "%s.cxx.ppm", argv[0] );
    }
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);

  controller->Initialize(argc, argv);
  controller->SetNumberOfProcesses(2);
  controller->SetMultipleMethod(1, process_a, NULL);
  controller->SetMultipleMethod(0, process_b, save_filename);
  controller->MultipleMethodExecute();

  controller->UnRegister(NULL);
  
  exit( 1 );
}





