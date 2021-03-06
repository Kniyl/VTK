#include "vtkSMPMinMaxTree.h"
#include "vtkObjectFactory.h"
#include "vtkDataSet.h"
#include "vtkParallelOperators.h"
#include "vtkThreadLocal.h"
#include "vtkTreeFunctor.h"
#include "vtkRangeFunctor.h"
#include "vtkRange1D.h"
#include "vtkTreeFunctor.h"
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

class InitializeFunctor : public vtkRangeFunctor
{
  InitializeFunctor( const InitializeFunctor& );
  void operator =( const InitializeFunctor& );

protected:
  InitializeFunctor()
    {
    this->TLS_Cell = vtkThreadLocal<vtkGenericCell>::New();
    this->TLS_CellScalars = vtkThreadLocal<vtkDoubleArray>::New();
    this->Locks = 0;
    }
  ~InitializeFunctor()
    {
    this->TLS_Cell->Delete();
    this->TLS_CellScalars->Delete();
    if ( this->Locks )
      delete [] this->Locks;
    }

  vtkScalarRange<double> *Tree;
  vtkIdType Size, BF, Offset, Max;
  vtkDataSet* DS;
  vtkDataArray* Scalars;
  vtkThreadLocal<vtkGenericCell>* TLS_Cell;
  vtkThreadLocal<vtkDoubleArray>* TLS_CellScalars;
  vtkIdType* Locks;

public:
  vtkTypeMacro(InitializeFunctor,vtkRangeFunctor);
  static InitializeFunctor* New();
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  void InitializeData( vtkSMPMinMaxTree* t )
    {
    Tree = static_cast<vtkScalarRange<double>*>(t->Tree);
    BF = t->BranchingFactor;
    DS = t->DataSet;
    Scalars = t->Scalars;
    Size = DS->GetNumberOfCells();
    Offset = t->LeafOffset;
    Locks = new vtkIdType[this->Offset];
    memset( this->Locks, 0, sizeof(vtkIdType) * this->Offset );
    Max = t->TreeSize;
    }

  void operator()( vtkRange* r ) const
    {
    vtkRange1D* range = vtkRange1D::SafeDownCast(r);
    int tid = range->GetTid();
    vtkGenericCell* cell = this->TLS_Cell->GetLocal(tid);
    if (!cell) cell = this->TLS_Cell->NewLocal(tid);
    vtkDoubleArray* cellScalars = this->TLS_CellScalars->GetLocal(tid);
    if (!cellScalars) cellScalars = this->TLS_CellScalars->GetLocal(tid);
    double my_min = VTK_DOUBLE_MAX;
    double my_max = -VTK_DOUBLE_MAX;
    double* s;

    vtkIdType start = range->Begin();
    vtkIdType cellId = ( start - this->Offset ) * this->BF;
    for (; start < range->End(); ++start)
      {
      vtkIdType index = start;
      if ( cellId < this->Size )
        {
        for ( vtkIdType i = 0; i < this->BF && cellId < this->Size; ++i, ++cellId )
          {
          this->DS->GetCell( cellId, cell );
          vtkIdList* cellPts = cell->GetPointIds();
          vtkIdType n = cellPts->GetNumberOfIds();
          cellScalars->SetNumberOfTuples( n );
          this->Scalars->GetTuples( cellPts, cellScalars );
          s = cellScalars->GetPointer( 0 );

          while ( n-- )
            {
            if ( s[n] < my_min )
              {
              my_min = s[n];
              }
            if ( s[n] > my_max )
              {
              my_max = s[n];
              }
            }
          }
        this->Tree[index].max = my_max;
        this->Tree[index].min = my_min;
        }

      while ( index )
        {
        index = ( index - 1 ) / this->BF;
        if ( __sync_add_and_fetch(&(this->Locks[index]), 1) != this->BF )
          break;
        for ( vtkIdType i = index * this->BF + 1; i < ( index + 1 ) * this->BF && i < this->Max; ++i )
          {
          if ( this->Tree[i].min < my_min )
            {
            my_min = this->Tree[i].min;
            }
          if ( this->Tree[i].max > my_max )
            {
            my_max = this->Tree[i].max;
            }
          }
        this->Tree[index].max = my_max;
        this->Tree[index].min = my_min;
        }
      }
    }
};

vtkStandardNewMacro(InitializeFunctor);

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
  vtkIdType numCells, cellId, i, j, numScalars;
  int level, offset, parentOffset, prod;
  vtkIdType numNodes, node, numLeafs, leaf, numParentLeafs;
  vtkCell *cell;
  vtkIdList *cellPts;
  vtkScalarRange<double> *tree, *parent;
  double *s;
  vtkDoubleArray *cellScalars;

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
  cellScalars = vtkDoubleArray::New();
  cellScalars->Allocate(100);

  // Compute the number of levels in the tree
  //
  numLeafs = ((numCells - 1) / this->BranchingFactor) + 1;
  for (prod=1, numNodes=1, this->Level=0;
       prod < numLeafs && this->Level <= this->MaxLevel; this->Level++ )
    {
    prod *= this->BranchingFactor;
    numNodes += prod;
    }

  this->LeafOffset = offset = numNodes - prod;
  vtkScalarRange<double> *TTree;
  this->TreeSize = offset + numLeafs;
  this->Tree = TTree = new vtkScalarRange<double>[this->TreeSize];
/*
  InitializeFunctor* InitTree = InitializeFunctor::New();
  InitTree->InitializeData( this );
  vtkParallelOperators::ForEach( offset, numNodes, InitTree );
  InitTree->Delete();
*/
  for ( i=0; i < this->TreeSize; i++ )
    {
    TTree[i].min = VTK_DOUBLE_MAX;
    TTree[i].max = -VTK_DOUBLE_MAX;
    }

  // Loop over all cells getting range of scalar data and place into leafs
  //
  for ( cellId=0, node=0; node < numLeafs; node++ )
    {
    tree = TTree + offset + node;
    for ( i=0; i < this->BranchingFactor && cellId < numCells; i++, cellId++ )
      {
      cell = this->DataSet->GetCell(cellId);
      cellPts = cell->GetPointIds();
      numScalars = cellPts->GetNumberOfIds();
      cellScalars->SetNumberOfTuples(numScalars);
      this->Scalars->GetTuples(cellPts, cellScalars);
      s = cellScalars->GetPointer(0);

      for ( j=0; j < numScalars; j++ )
        {
        if ( s[j] < tree->min )
          {
          tree->min = s[j];
          }
        if ( s[j] > tree->max )
          {
          tree->max = s[j];
          }
        }
      }
    }

  // Now build top levels of tree in bottom-up fashion
  //
  for ( level=this->Level; level > 0; level-- )
    {
    parentOffset = offset - prod/this->BranchingFactor;
    prod /= this->BranchingFactor;
    numParentLeafs = static_cast<int>(
      ceil(static_cast<double>(numLeafs)/this->BranchingFactor));

    for ( leaf=0, node=0; node < numParentLeafs; node++ )
      {
      parent = TTree + parentOffset + node;
      for ( i=0; i < this->BranchingFactor && leaf < numLeafs; i++, leaf++ )
        {
        tree = TTree + offset + leaf;
        if ( tree->min < parent->min )
          {
          parent->min = tree->min;
          }
        if ( tree->max > parent->max )
          {
          parent->max = tree->max;
          }
        }
      }

    numLeafs = numParentLeafs;
    offset = parentOffset;
    }

  this->BuildTime.Modified();
  cellScalars->Delete();
  }

void vtkSMPMinMaxTree::InitTraversal(double scalarValue)
  {
  this->BuildTree();
  vtkScalarRange<double> *TTree =
      static_cast< vtkScalarRange<double> * > (this->Tree);

  this->ScalarValue = scalarValue;
  this->TreeIndex = this->TreeSize;
  }

int vtkSMPMinMaxTree::TraverseNode( vtkIdType id, int lvl, vtkTreeFunctor* function, vtkLocalData* data ) const
  {
  if ( id >= this->TreeSize )
    {
    return 0;
    }

  vtkScalarRange<double> *t = static_cast<vtkScalarRange<double>*>(this->Tree) + id;
  if ( t->min > this->ScalarValue || t->max < this->ScalarValue )
    {
    return 0;
    }

  if ( lvl == this->Level ) //leaf
    {
    vtkIdType cell_id = ( id - this->LeafOffset ) * this->BranchingFactor;
    vtkIdType max_id = this->DataSet->GetNumberOfCells();
    for ( vtkIdType i = 0; i < this->BranchingFactor && cell_id < max_id; ++i, ++cell_id )
      {
      (*function)( cell_id, data );
      }
    return 0;
    }
  else //node
    {
    return 1;
    }
  }

void vtkSMPMinMaxTree::GetTreeSize( int& max_level, vtkIdType& branching_factor ) const
  {
  max_level = this->Level;
  branching_factor = this->BranchingFactor;
  }
