#include "vtkLockPointMerger.h"

#include "vtkDummyMergeFunctor.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSMPMergePoints.h"
#include "vtkThreadLocal.h"

//------------------------------------------------------------------------------
vtkStandardNewMacro(vtkLockPointMerger);

//------------------------------------------------------------------------------
void vtkLockPointMerger::PrintSelf(ostream &os, vtkIndent indent)
  {
  this->Superclass::PrintSelf(os,indent);
  }

//------------------------------------------------------------------------------
void vtkLockPointMerger::operator()( vtkRange* r ) const
  {
  vtkRange1D* range = vtkRange1D::SafeDownCast(r);
  vtkThreadLocal<vtkPoints>::iterator itPoints = this->Functor->InPoints->Begin();
  vtkThreadLocal<vtkPointData>::iterator itPd = this->Functor->InPd->Begin();
  vtkThreadLocal<vtkIdList>::iterator itMaps = this->Functor->Maps->Begin();

  vtkIdType id = range->Begin();
  vtkIdType NumberOfPoints = NumberOfPointsFirstThread, NewId;
  while ( id >= NumberOfPoints )
    {
    id -= NumberOfPoints;
    ++itPoints; ++itPd; ++itMaps;
    NumberOfPoints = (*itPoints)->GetNumberOfPoints();
    }

  double* pt = new double[3];
  for (vtkIdType b = range->Begin(); b < range->End(); ++b, ++id)
    {
    (*itPoints)->GetPoint( id, pt );
    if ( this->Functor->outputLocator->SetUniquePoint( pt, NewId ) )
      this->Functor->outputPd->SetTuple( NewId, id, (*itPd) );
    (*itMaps)->SetId( id, NewId );
    }
  delete [] pt;
  }
