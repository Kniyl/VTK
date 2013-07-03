#include "vtkActor.h"
#include "vtkAssignAttribute.h"
#include "vtkDataSetMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRTAnalyticSource.h"
#include "vtkSMPContourFilter.h"
#include "vtkSMPMergePoints.h"
#include "vtkSMPMinMaxTree.h"
#include "vtkSMPTransform.h"
#include "vtkTestUtilities.h"
#include "vtkThreshold.h"
#include "vtkTimerLog.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkUnstructuredGrid.h"

#include <cstdlib>

const int WaveSize = 100;
const int NumRepetitions = 5;

void testUG(vtkAlgorithm* input, vtkContourFilter* isosurface, bool sequential = true)
  {
  vtkTimerLog *timer = vtkTimerLog::New();
  double t,t0,t1;

  cerr << "/************* ";
  if (sequential)
   cerr << "Sequential";
  else
   cerr << "SMP";
  cerr << " *************/" << endl;

  vtkTransformFilter* transform = vtkTransformFilter::New();
  transform->SetInputConnection( input->GetOutputPort() );
  if (sequential)
    {
    vtkTransform* tr = vtkTransform::New();
    tr->Identity();
    transform->SetTransform(tr);
    tr->Delete();
    }
  else
    {
    vtkSMPTransform* tr = vtkSMPTransform::New();
    tr->Identity();
    transform->SetTransform(tr);
    tr->Delete();
    }

  isosurface->SetInputConnection( transform->GetOutputPort() );
  transform->Delete();
  isosurface->GenerateValues( 11, 30.0, 300.0 );
  isosurface->UseScalarTreeOn();
  cerr << "First " << transform->GetClassName() << " execution" << endl;
  t0 = timer->GetUniversalTime();
  transform->Update();
  t1 = timer->GetUniversalTime();
  cerr << t1-t0 << endl;
  cerr << "Average time for " << NumRepetitions << " other executions" << endl;
  t = 0.0;
  for (int i = 0; i < NumRepetitions; ++i)
    {
    transform->Modified();
    t0 = timer->GetUniversalTime();
    transform->Update();
    t1 = timer->GetUniversalTime();
    t += t1-t0;
    }
  cerr << "Transform: " << (t)/NumRepetitions << endl;

  cerr << "First " << isosurface->GetClassName() << " execution" << endl;
  t0 = timer->GetUniversalTime();
  isosurface->Update();
  t1 = timer->GetUniversalTime();
  cerr << t1-t0 << endl;
  cerr << "Average time for " << NumRepetitions << " other executions" << endl;
  t = 0.0;
  for (int i = 0; i < NumRepetitions; ++i)
    {
    isosurface->Modified();
    t0 = timer->GetUniversalTime();
    isosurface->Update();
    t1 = timer->GetUniversalTime();
    t += t1-t0;
    }
  cerr << "Isosurface: " << (t)/NumRepetitions << endl;
  timer->Delete();
  }

int TestSMPUG( int argc, char * argv [] )
  {
  vtkRTAnalyticSource* wavelet = vtkRTAnalyticSource::New();
  wavelet->SetWholeExtent(-WaveSize,WaveSize,-WaveSize,WaveSize,-WaveSize,WaveSize);
  wavelet->SetCenter(0,0,0);
  wavelet->SetMaximum(255);
  wavelet->SetXFreq(60);
  wavelet->SetYFreq(30);
  wavelet->SetZFreq(40);
  wavelet->SetXMag(10);
  wavelet->SetYMag(18);
  wavelet->SetZMag(5);
  wavelet->SetStandardDeviation(0.5);
  wavelet->SetSubsampleRate(1);

  //convert to unstructured grid
  vtkThreshold* threshold = vtkThreshold::New();
  threshold->SetInputConnection( wavelet->GetOutputPort() );
  threshold->ThresholdBetween(0,3000);
  wavelet->Delete();

  vtkAssignAttribute* aa = vtkAssignAttribute::New();
  aa->SetInputConnection(threshold->GetOutputPort());
  aa->Assign("RTData", vtkDataSetAttributes::SCALARS, vtkAssignAttribute::POINT_DATA);

  /* === Testing contour filter === */

  vtkContourFilter* isosurface1 = vtkContourFilter::New();
  
  vtkSMPContourFilter* isosurface2 = vtkSMPContourFilter::New();
  vtkSMPMergePoints* locator = vtkSMPMergePoints::New();
  isosurface2->SetLocator( locator );
  locator->Delete();
  vtkSMPMinMaxTree* tree = vtkSMPMinMaxTree::New();
  isosurface2->SetScalarTree(tree);
  tree->Delete();

  aa->Update(); // Pull pipeline up-to-date to acurately time further filters
  testUG(aa, isosurface1);
  testUG(aa, isosurface2, false);

  /* === Watching results === */
  
  double b[6];
  isosurface1->GetOutput()->GetBounds(b);

  vtkDataSetMapper* map1 = vtkDataSetMapper::New();
  map1->SetInputConnection( isosurface1->GetOutputPort() );
  vtkActor* actor1 = vtkActor::New();
  actor1->SetMapper( map1 );
  map1->Delete();

  vtkDataSetMapper* map2 = vtkDataSetMapper::New();
  map2->SetInputConnection( isosurface2->GetOutputPort() );
  vtkActor* actor2 = vtkActor::New();
  actor2->SetMapper( map2 );
  map2->Delete();
  actor2->AddPosition( (b[1]-b[0])*1.1, 0., 0. );

  vtkRenderer* viewport = vtkRenderer::New();
  viewport->SetBackground( .5, .5, .5 );
  viewport->AddActor( actor1 );
  viewport->AddActor( actor2 );
  actor1->Delete();
  actor2->Delete();

  vtkRenderWindow* window = vtkRenderWindow::New();
  window->AddRenderer( viewport );
  viewport->Delete();

  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow( window );
  window->Delete();

  iren->Initialize();

  iren->Start();

  iren->Delete();

  isosurface1->Delete();
  isosurface2->Delete();

  return 0;
  }
