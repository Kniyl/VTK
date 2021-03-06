#ifndef _vtkOffsetManager_h_
#define _vtkOffsetManager_h_

//internal functor used in mergeoperator used to deal with 4 cell arrays in vtkPolyData

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkObject.h"
#include <vector>

class vtkCellArray;

class VTKPARALLELSMP_EXPORT vtkOffsetManager : public vtkObject
{
  vtkOffsetManager( const vtkOffsetManager& );
  void operator =( const vtkOffsetManager& );
  std::vector<vtkIdType> cells;
  std::vector<vtkIdType> tuples;
  vtkIdType CellsOffset;
  vtkIdType TuplesOffset;
  std::vector<vtkIdType>::iterator itCells;
  std::vector<vtkIdType>::iterator itTuples;

protected:
  vtkOffsetManager();
  ~vtkOffsetManager() { }
public:
  vtkTypeMacro(vtkOffsetManager,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkOffsetManager* New();

  void InitManageValues ();

  void ManageNextValue ( vtkCellArray* ca );

  vtkIdType GetNumberOfCells() { return CellsOffset; }
  vtkIdType GetNumberOfTuples() { return TuplesOffset; }

  std::vector<vtkIdType>::iterator GetCellsOffset ( ) { return cells.begin(); }
  std::vector<vtkIdType>::iterator GetTuplesOffset ( ) { return tuples.begin(); }
};

#endif //_vtkOffsetManager_h_
