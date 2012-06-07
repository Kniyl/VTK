#include "vtkSMPMinMaxTree.h"
#include "vtkObjectFactory.h"
#include "vtkDataSet.h"
#include "vtkSMP.h"
#include "vtkIdList.h"
#include "vtkGenericCell.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"

class vtkScalarNode {};

template <class TScalar>
class vtkScalarRange : public vtkScalarNode
{
public:
  TScalar min;
  TScalar max;
};

class BuildFunctor : public vtkFunctor
{
  BuildFunctor( const BuildFunctor& );
  void operator =( const BuildFunctor& );

  vtkScalarRange<double>* Tree;
  vtkIdType BF, TreeSize;

protected:
  BuildFunctor() { }
  ~BuildFunctor() { }

public:

  vtkTypeMacro(BuildFunctor,vtkFunctor);
  static BuildFunctor* New();
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  void InitializeData( vtkSMPMinMaxTree* t )
    {
    Tree = static_cast<vtkScalarRange<double>*>(t->Tree);
    BF = t->BranchingFactor;
    TreeSize = t->TreeSize;
    }

  void operator ()( vtkIdType index, vtkSMPThreadID tid ) const
    {
    vtkIdType branchIndex = (index * BF) + 1;
    vtkScalarRange<double> *t = Tree + index, *s = Tree + branchIndex;
    for ( vtkIdType i = 0; i < BF && branchIndex < TreeSize; ++i, ++branchIndex )
      {
      if ( t->min > s[i].min )
        t->min = s[i].min;
      if ( t->max < s[i].max )
        t->max = s[i].max;
      }
    }
};

vtkStandardNewMacro(BuildFunctor);

class BuildLeafFunctor : public vtkFunctor
{
  BuildLeafFunctor( const BuildLeafFunctor& );
  void operator =( const BuildLeafFunctor& );

  vtkScalarRange<double>* Tree;
  vtkIdType BF, TreeSize, NbCells;
  vtkDataSet* DS;
  vtkDataArray* Scalars;

protected:
  BuildLeafFunctor() { }
  ~BuildLeafFunctor() { }

public:

  vtkTypeMacro(BuildLeafFunctor,vtkFunctor);
  static BuildLeafFunctor* New();
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  void InitializeData( vtkSMPMinMaxTree* t )
    {
    Tree = static_cast<vtkScalarRange<double>*>(t->Tree);
    BF = t->BranchingFactor;
    TreeSize = t->LeafOffset * t->BranchingFactor + 1;
    NbCells = t->DataSet->GetNumberOfCells();
    DS = t->DataSet;
    Scalars = t->Scalars;
    }

  void operator ()( vtkIdType index, vtkSMPThreadID tid ) const
    {
    vtkScalarRange<double> *t = Tree + index;
    t->min = VTK_DOUBLE_MAX;
    t->max = -VTK_DOUBLE_MAX;
    vtkIdType cellId = (index * BF) + 1 - TreeSize;
    vtkIdType n;
    vtkGenericCell* cell = vtkGenericCell::New();
    vtkDoubleArray* cellScalars = vtkDoubleArray::New();
    double* pt = new double[3];
    double* s;
    for ( vtkIdType i = 0; i < BF && cellId < NbCells; ++i, ++cellId )
      {
      DS->GetCell(cellId, cell);
      vtkIdList* cellPts = cell->GetPointIds();
      n = cellPts->GetNumberOfIds();
      cellScalars->SetNumberOfTuples( n );
      Scalars->GetTuples(cellPts, cellScalars);
      s = cellScalars->GetPointer(0);

      while ( n != 0 )
        {
        --n;
        if ( s[n] < t->min )
          {
          t->min = s[n];
          }
        if ( s[n] > t->max )
          {
          t->max = s[n];
          }
        }
      }
    delete [] pt;
    cell->Delete();
    cellScalars->Delete();
    }
};

vtkStandardNewMacro(BuildLeafFunctor);

vtkStandardNewMacro(vtkSMPMinMaxTree);

vtkSMPMinMaxTree::vtkSMPMinMaxTree()
  {
  }

vtkSMPMinMaxTree::~vtkSMPMinMaxTree()
  {
  }

void vtkSMPMinMaxTree::PrintSelf(ostream &os, vtkIndent indent)
  {
  this->Superclass::PrintSelf(os, indent);
  }

void vtkSMPMinMaxTree::BuildTree()
  {
  vtkIdType numCells;
  int offset, prod;
  vtkIdType numNodes, numLeafs;

  // Check input...see whether we have to rebuild
  //
  if ( !this->DataSet || (numCells = this->DataSet->GetNumberOfCells()) < 1 )
    {
    vtkErrorMacro( << "No data to build tree with");
    return;
    }

  if ( this->Tree != NULL && this->BuildTime > this->MTime
       && this->BuildTime > this->DataSet->GetMTime() )
    {
    return;
    }

  vtkDebugMacro( << "Building scalar tree..." );

  this->Scalars = this->DataSet->GetPointData()->GetScalars();
  if ( ! this->Scalars )
    {
    vtkErrorMacro( << "No scalar data to build trees with");
    return;
    }

  this->Initialize();

  // Compute the number of levels in the tree
  //
  numLeafs = static_cast<int>(
        ceil(static_cast<double>(numCells)/this->BranchingFactor));
  for (prod=1, numNodes=1, this->Level=0;
       prod < numLeafs && this->Level <= this->MaxLevel; this->Level++ )
    {
    prod *= this->BranchingFactor;
    numNodes += prod;
    }

  this->LeafOffset = offset = numNodes - prod;
  vtkScalarRange<double> *TTree;
  this->TreeSize = numNodes - (prod - numLeafs);
  this->Tree = TTree = new vtkScalarRange<double>[this->TreeSize];

  BuildLeafFunctor* BuildLeaf = BuildLeafFunctor::New();
  BuildLeaf->InitializeData( this );
  vtkIdType StartIndex = this->LeafOffset;
  vtkIdType EndIndex = this->TreeSize;
  vtkSMP::ForEach( StartIndex, EndIndex, BuildLeaf );
  BuildLeaf->Delete();

  BuildFunctor* BuildCells = BuildFunctor::New();
  BuildCells->InitializeData( this );
  while ( StartIndex != 0 )
    {
    EndIndex = StartIndex;
    StartIndex = ( StartIndex - 1 ) / this->BranchingFactor;
    vtkSMP::ForEach( StartIndex, EndIndex, BuildCells );
    }
  BuildCells->Delete();

  this->BuildTime.Modified();
  }

void vtkSMPMinMaxTree::InitTraversal(double scalarValue)
  {
  this->BuildTree();
  vtkScalarRange<double> *TTree =
      static_cast< vtkScalarRange<double> * > (this->Tree);

  this->ScalarValue = scalarValue;
  this->TreeIndex = this->TreeSize;
  }


vtkIdType vtkSMPMinMaxTree::GetAncestor( vtkIdType id, int lvl, int desiredLvl ) const
  {
  vtkIdType result = id;

  while ( lvl < desiredLvl )
    {
    result = ( result - 1 ) / this->BranchingFactor;
    --lvl;
    }

  return result;
  }

vtkIdType vtkSMPMinMaxTree::GetLastDescendant( vtkIdType id, int lvl ) const
  {
  return ( id + 1 ) * this->BranchingFactor;
  }

// In place transformation of indices and levels
void vtkSMPMinMaxTree::GetNextStealableNode( vtkIdType* stealedId, int* stealedLvl ) const
  {
redo:
  if ( *stealedId == 0 )
    {
    *stealedId = -1;
    *stealedLvl = -1;
    return;
    }
  vtkIdType father = *stealedId / this->BranchingFactor;
  vtkIdType previous_father = ( *stealedId - 1 ) / this->BranchingFactor;
  if ( previous_father != father )
    {
    *stealedId = previous_father;
    --(*stealedLvl);
    goto redo;
    }

  ++(*stealedId);
  }

void vtkSMPMinMaxTree::TraverseNode( vtkIdType* index, int* level, vtkFunctor* function, vtkSMPThreadID tid ) const
  {
  if ( *index >= this->TreeSize )
    {
    *index = -1;
    *level = -1;
    }
  else
    {
    vtkScalarRange<double> *t = static_cast<vtkScalarRange<double>*>(this->Tree) + *index;
    if ( t->min > this->ScalarValue || t->max < this->ScalarValue )
      {
      this->GetNextStealableNode( index, level );
      }
    else
      {
      if ( *level == this->Level )
        {
        vtkIdType cell_id = ( *index - this->LeafOffset ) * this->BranchingFactor;
        vtkIdType max_id = this->DataSet->GetNumberOfCells();
        for ( vtkIdType i = 0; i < this->BranchingFactor && cell_id < max_id; ++i, ++cell_id )
          {
          (*function)( cell_id, tid );
          }
        this->GetNextStealableNode( index, level );
        }
      else
        {
        *index = ( *index * this->BranchingFactor ) + 1;
        ++(*level);
        }
      }
    }
  }
