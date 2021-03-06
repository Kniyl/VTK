#include "vtkAssignAttribute.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkTransformFilter.h"
#include "vtkTransform.h"
#include "vtkSMPTransform.h"
#include "vtkSMPContourFilter.h"
#include "vtkSMPMergePoints.h"
#include "vtkSMPMinMaxTree.h"

#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include <cstdlib>

#include "vtkUnstructuredGridWriter.h"
#include "vtkSubdivideTetra.h"
#include "vtkMergePoints.h"
#include "vtkCellArray.h"

int main( int argc, char** argv )
{
  if ( argc < 2 )
  {
    cout << "You must provide a file name" << endl;
    return 1;
  }
  int parallel = argc == 2 ? 48 : atoi(argv[2]);
  cerr << "parallel " << parallel << endl;
  /* === Reading 3D model === */

  vtkUnstructuredGridReader* usgReader = 0;
  vtkPolyDataReader* polyReader = vtkPolyDataReader::New();
  polyReader->SetFileName(argv[1]);
  if ( !polyReader->IsFilePolyData() )
    {
    usgReader = vtkUnstructuredGridReader::New();
    usgReader->SetFileName(argv[1]);
    polyReader->CloseVTKFile();
    polyReader->Delete();
    polyReader = 0;
    if ( !usgReader->IsFileUnstructuredGrid() )
      {
      usgReader->Delete();
      cout << argv[1] << " is not a suitable file" << endl;
      return 1;
      }
    else cerr << "UG" << endl;
    }
  else cerr << "PD" << endl;
  cout << "Using file " << argv[1] << endl;

  vtkTransformFilter* pre_transform = vtkTransformFilter::New();
  if ( !usgReader )
    {
    pre_transform->SetInputConnection( polyReader->GetOutputPort() );
    polyReader->Delete();
    }
  else
    {
    pre_transform->SetInputConnection( usgReader->GetOutputPort() );
    usgReader->Delete();
    }

  /*
  vtkAssignAttribute* aa = vtkAssignAttribute::New();
  aa->SetInputConnection(pre_transform->GetOutputPort());
  aa->Assign("RTData", vtkDataSetAttributes::SCALARS, vtkAssignAttribute::POINT_DATA);
  pre_transform->Delete();
*/

  vtkTransformFilter* transform = vtkTransformFilter::New();
  transform->SetInputConnection( usgReader->GetOutputPort() );
  pre_transform->Delete();

  if ( parallel < 3 )
    {
    cerr << "SMPT" << endl;
    vtkSMPTransform* t = vtkSMPTransform::New();
    t->Scale( 1., 2., 3. );
    pre_transform->SetTransform( t );
    transform->SetTransform( t );
    t->Delete();
    }
  else
    {
    cerr << "T" << endl;
    vtkTransform* t = vtkTransform::New();
    t->Scale( 1., 2., 3. );
    pre_transform->SetTransform( t );
    transform->SetTransform( t );
    t->Delete();
    }

  /* === Testing contour filter === */
  transform->Update();
  if ( !usgReader )
    {
    vtkPointSet* data = transform->GetOutput();
    vtkDataArray* s = vtkDoubleArray::New();
    s->SetNumberOfComponents(1);
    s->SetNumberOfTuples(data->GetNumberOfPoints());
    s->SetName("scalars");
    vtkIdType num = data->GetNumberOfCells() / parallel + 1, n;
    vtkGenericCell* cell = vtkGenericCell::New();

    for ( vtkIdType i = 0; i < data->GetNumberOfCells(); ++i )
      {
      data->GetCell( i, cell );
      n = cell->GetNumberOfPoints();
      int v = i < num;
      while ( n-- )
        {
        s->SetTuple1( cell->GetPointId( n ), v );
        }
      }
    cell->Delete();

    data->GetPointData()->SetScalars( s );
    s->Delete();
    }

  vtkContourFilter* isosurface = parallel != 1 ? vtkSMPContourFilter::New() : vtkContourFilter::New();
  cerr << isosurface->GetClassName() << endl;
  if ( parallel < 3 )
    {
    cerr << "SMPMERGE" << endl;
    vtkSMPMergePoints* locator = vtkSMPMergePoints::New();
    isosurface->SetLocator( locator );
    locator->Delete();
    vtkSMPMinMaxTree* tree = vtkSMPMinMaxTree::New();
    isosurface->SetScalarTree(tree);
    tree->Delete();
    }
  else cerr << "NO MERGE" << endl;
  isosurface->SetInputConnection( transform->GetOutputPort() );
  isosurface->SetInputArrayToProcess(0,0,0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RTData");
  isosurface->GenerateValues( 11, 30.0, 300.0 );
  isosurface->UseScalarTreeOn();
  isosurface->Update();

  vtkDataSetMapper* dmap = vtkDataSetMapper::New();
  dmap->SetInputConnection( usgReader->GetOutputPort() );

  vtkActor* dobject = vtkActor::New();
  dobject->SetMapper( dmap );
  dmap->Delete();

  vtkPolyDataMapper* map = vtkPolyDataMapper::New();
  map->SetInputConnection( isosurface->GetOutputPort() );
  isosurface->Delete();

  vtkActor* object = vtkActor::New();
  object->SetMapper( map );
  map->Delete();

  vtkActor* reference = vtkActor::New();
  if ( !usgReader )
    {
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( transform->GetOutputPort() );
    reference->SetMapper( mapper );
    mapper->Delete();
    }
  else
    {
    vtkDataSetMapper* mapper = vtkDataSetMapper::New();
    mapper->SetInputConnection( transform->GetOutputPort() );
    reference->SetMapper( mapper );
    mapper->Delete();
    }
  double b[6];
  transform->GetOutput()->GetBounds(b);
  reference->AddPosition( (b[1]-b[0])*1.02, 0., 0. );

  vtkRenderer* viewport = vtkRenderer::New();
  viewport->SetBackground( .5, .5, .5 );
  //viewport->AddActor( dobject );
  viewport->AddActor( object );
  viewport->AddActor( reference );
  dobject->Delete();
  object->Delete();
  reference->Delete();

  vtkRenderWindow* window = vtkRenderWindow::New();
  window->AddRenderer( viewport );
  viewport->Delete();

  vtkRenderWindowInteractor* eventsCatcher = vtkRenderWindowInteractor::New();
  eventsCatcher->SetRenderWindow( window );
  window->Delete();

  transform->Delete();

  eventsCatcher->Initialize();
  eventsCatcher->Start();

  eventsCatcher->Delete();

  return 0;
}
