#include "vtkSMPContourFilter.h"
#include "vtkObjectFactory.h"
#include "vtkMergeDataSets.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCommand.h"
#include "vtkContourGrid.h"
#include "vtkCutter.h"
#include "vtkGenericCell.h"
#include "vtkGridSynchronizedTemplates3D.h"
#include "vtkImageData.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMergePoints.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkRectilinearSynchronizedTemplates.h"
#include "vtkScalarTree.h"
#include "vtkSimpleScalarTree.h"
#include "vtkParallelOperators.h"
#include "vtkRangeFunctorInitializable.h"
#include "vtkRange1D.h"
#include "vtkTreeFunctorInitializable.h"
#include "vtkSMPMergePoints.h"
#include "vtkSMPMinMaxTree.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkSynchronizedTemplates2D.h"
#include "vtkSynchronizedTemplates3D.h"
#include "vtkUniformGrid.h"

vtkStandardNewMacro(vtkSMPContourFilter);

vtkSMPContourFilter::vtkSMPContourFilter() : vtkContourFilter() { }
vtkSMPContourFilter::~vtkSMPContourFilter() { }

void vtkSMPContourFilter::PrintSelf(ostream& os, vtkIndent indent)
  {
  this->Superclass::PrintSelf(os,indent);
  }

/* ================================================================================
  Local data for tree traversal
 ================================================================================ */
class vtkSMPContourFilterLocalData : public vtkLocalData
  {
    vtkSMPContourFilterLocalData(const vtkSMPContourFilterLocalData&);
    void operator=(const vtkSMPContourFilterLocalData&);
  protected:
    vtkSMPContourFilterLocalData() :
        Locator(0), newVerts(0), newLines(0),
        newPolys(0), outPd(0), outCd(0),
        Cells(0), CellsScalars(0)
      {
      }
    ~vtkSMPContourFilterLocalData(){}
  public:
    vtkTypeMacro(vtkSMPContourFilterLocalData,vtkLocalData);
    static vtkSMPContourFilterLocalData* New();
    virtual void PrintSelf(ostream& os, vtkIndent indent)
      {
      this->Superclass::PrintSelf(os,indent);
      os << indent << "Locator: (" << Locator << ")" << endl;
      if (Locator) Locator->PrintSelf(os,indent.GetNextIndent());
      os << indent << "newVerts: (" << newVerts << ")" << endl;
      if (newVerts) newVerts->PrintSelf(os,indent.GetNextIndent());
      os << indent << "newLines: (" << newLines << ")" << endl;
      if (newLines) newLines->PrintSelf(os,indent.GetNextIndent());
      os << indent << "newPolys: (" << newPolys << ")" << endl;
      if (newPolys) newPolys->PrintSelf(os,indent.GetNextIndent());
      os << indent << "outPd: (" << outPd << ")" << endl;
      if (outPd) outPd->PrintSelf(os,indent.GetNextIndent());
      os << indent << "outCd: (" << outCd << ")" << endl;
      if (outCd) outCd->PrintSelf(os,indent.GetNextIndent());
      os << indent << "Cells: (" << Cells << ")" << endl;
      if (Cells) Cells->PrintSelf(os,indent.GetNextIndent());
      os << indent << "CellsScalars: (" << CellsScalars << ")" << endl;
      if (CellsScalars) CellsScalars->PrintSelf(os,indent.GetNextIndent());
      }

    vtkIncrementalPointLocator* Locator;
    vtkCellArray* newVerts;
    vtkCellArray* newLines;
    vtkCellArray* newPolys;
    vtkPointData* outPd;
    vtkCellData* outCd;
    vtkGenericCell* Cells;
    vtkDataArray* CellsScalars;
  };

vtkStandardNewMacro(vtkSMPContourFilterLocalData);

/* ================================================================================
  Generic contouring: Functors for parallel execution without ScalarTree
 ================================================================================ */
class ContourRangeFunctor : public vtkRangeFunctorInitializable
{
  ContourRangeFunctor( const ContourRangeFunctor& );
  void operator =( const ContourRangeFunctor& );

protected:
  ContourRangeFunctor()
    {
    Locator = vtkThreadLocal<vtkIncrementalPointLocator>::New();
    newPts = vtkThreadLocal<vtkPoints>::New();
    newVerts = vtkThreadLocal<vtkCellArray>::New();
    newLines = vtkThreadLocal<vtkCellArray>::New();
    newPolys = vtkThreadLocal<vtkCellArray>::New();
    outPd = vtkThreadLocal<vtkPointData>::New();
    outCd = vtkThreadLocal<vtkCellData>::New();
    Cells = vtkThreadLocal<vtkGenericCell>::New();
    CellsScalars = vtkThreadLocal<vtkDataArray>::New();
    }

  ~ContourRangeFunctor()
    {
    Locator->Delete();
    newPts->Delete();
    newVerts->Delete();
    newLines->Delete();
    newPolys->Delete();
    outPd->Delete();
    outCd->Delete();
    Cells->Delete();
    CellsScalars->Delete();
    }

public:
  vtkTypeMacro(ContourRangeFunctor,vtkRangeFunctorInitializable);
  static ContourRangeFunctor* New();
  virtual void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  vtkThreadLocal<vtkIncrementalPointLocator>* Locator;
  vtkThreadLocal<vtkPoints>* newPts;
  vtkThreadLocal<vtkCellArray>* newVerts;
  vtkThreadLocal<vtkCellArray>* newLines;
  vtkThreadLocal<vtkCellArray>* newPolys;
  vtkThreadLocal<vtkPointData>* outPd;
  vtkThreadLocal<vtkCellData>* outCd;
  vtkThreadLocal<vtkGenericCell>* Cells;
  vtkThreadLocal<vtkDataArray>* CellsScalars;

  vtkDataArray* inScalars;
  vtkDataSet* input;
  vtkPointData* inPd;
  vtkCellData* inCd;
  vtkIdType estimatedSize;
  double* values;
  int numContours;
  int computeScalars;

  vtkIncrementalPointLocator* refLocator;
  vtkCellArray* outputVerts;
  vtkCellArray* outputLines;
  vtkCellArray* outputPolys;
  vtkCellData* outputCd;
  vtkPointData* outputPd;

  unsigned char cellTypeDimensions[VTK_NUMBER_OF_CELL_TYPES];

  int dimensionality;

  void SetData( vtkDataSet* _input, vtkPoints* _inPts, vtkCellData* _incd,
                  vtkPointData* _inpd, vtkIncrementalPointLocator* _locator,
                  vtkIdType& _size, double* _values, int _number,
                  vtkDataArray* _scalars, int _compute, vtkCellArray* _outputVerts,
                  vtkCellArray* _outputLines, vtkCellArray* _outputPolys,
                  vtkCellData* _outputCd, vtkPointData* _outputPd )
    {
    vtkCutter::GetCellTypeDimensions(cellTypeDimensions);

    this->inCd = _incd;
    this->inPd = _inpd;
    this->input = _input;
    this->refLocator = _locator;
    this->estimatedSize = _size;
    this->values = _values;
    this->numContours = _number;
    this->inScalars = _scalars;
    this->computeScalars = _compute;

    outputVerts = _outputVerts;
    outputLines = _outputLines;
    outputPolys = _outputPolys;
    outputCd = _outputCd;
    outputPd = _outputPd;

    int tid = this->MasterThreadId;
    if ( !computeScalars )
      {
      _outputPd->CopyScalarsOff();
      }
    _outputPd->InterpolateAllocate( inPd, estimatedSize, estimatedSize );
    _outputCd->CopyAllocate( inCd, estimatedSize, estimatedSize );
    Locator->SetLocal( _locator, tid );
    newPts->SetLocal( _inPts, tid );
    newVerts->SetLocal( _outputVerts, tid );
    newLines->SetLocal( _outputLines, tid );
    newPolys->SetLocal( _outputPolys, tid );
    outPd->SetLocal( _outputPd, tid );
    outCd->SetLocal( _outputCd, tid );
    Cells->NewLocal(tid);
    vtkDataArray* cScalars = CellsScalars->NewLocal( inScalars, tid );
    cScalars->SetNumberOfComponents( inScalars->GetNumberOfComponents() );
    cScalars->Allocate( cScalars->GetNumberOfComponents() * VTK_CELL_SIZE );

    Initialized(tid);
    }

  void operator()( vtkRange* r ) const
    {
    vtkRange1D* range = vtkRange1D::SafeDownCast(r);
    int tid = range->GetTid();
    vtkIncrementalPointLocator* Locator = this->Locator->GetLocal(tid);
    vtkCellArray* newVerts = this->newVerts->GetLocal(tid);
    vtkCellArray* newLines = this->newLines->GetLocal(tid);
    vtkCellArray* newPolys = this->newPolys->GetLocal(tid);
    vtkPointData* outPd = this->outPd->GetLocal(tid);
    vtkCellData* outCd = this->outCd->GetLocal(tid);
    vtkGenericCell* cell = this->Cells->GetLocal(tid);
    vtkDataArray* cellScalars = this->CellsScalars->GetLocal(tid);

    for (vtkIdType cellId = range->Begin(); cellId < range->End(); ++cellId)
      {
      int cellType = input->GetCellType(cellId);
      if (cellType >= VTK_NUMBER_OF_CELL_TYPES)
        { // Protect against new cell types added.
  //      vtkErrorMacro("Unknown cell type " << cellType);
        return;
        }
      if (cellTypeDimensions[cellType] != dimensionality)
        {
        return;
        }
      input->GetCell(cellId,cell);
      vtkIdList* cellPts = cell->GetPointIds();
      if ( cellScalars->GetSize() / cellScalars->GetNumberOfComponents() < cellPts->GetNumberOfIds() )
        {
        cellScalars->Allocate(cellScalars->GetNumberOfComponents()*cellPts->GetNumberOfIds());
        }
      inScalars->GetTuples(cellPts,cellScalars);

      for (int i = 0; i < numContours; i++)
        {
        double v = values[i];
        cell->Contour( v, cellScalars, Locator,
                       newVerts, newLines, newPolys,
                       inPd, outPd, inCd, cellId, outCd );
        }
      }
    }

  virtual void Init(int tid) const
    {
    vtkPoints* pts = this->newPts->NewLocal(tid);
    pts->Allocate( this->estimatedSize, this->estimatedSize );
    vtkIncrementalPointLocator* l = this->Locator->NewLocal( this->refLocator, tid );
    l->InitPointInsertion( pts, this->input->GetBounds(), this->estimatedSize );

    vtkCellArray* c = this->newVerts->NewLocal(tid);
    c->Allocate( this->estimatedSize, this->estimatedSize );
    c = this->newLines->NewLocal(tid);
    c->Allocate( this->estimatedSize, this->estimatedSize );
    c = this->newPolys->NewLocal(tid);
    c->Allocate( this->estimatedSize, this->estimatedSize );

    vtkPointData* pd = this->outPd->NewLocal(tid);
    if ( !this->computeScalars )
      {
      pd->CopyScalarsOff();
      }
    pd->InterpolateAllocate( this->inPd, this->estimatedSize, this->estimatedSize );

    vtkCellData* cd = this->outCd->NewLocal(tid);
    cd->CopyAllocate( this->inCd, this->estimatedSize, this->estimatedSize );

    this->Cells->NewLocal(tid);

    vtkDataArray* cScalars = this->CellsScalars->NewLocal( this->inScalars, tid );
    cScalars->SetNumberOfComponents( this->inScalars->GetNumberOfComponents() );
    cScalars->Allocate( cScalars->GetNumberOfComponents() * VTK_CELL_SIZE );

    Initialized(tid);
    }
};

vtkStandardNewMacro(ContourRangeFunctor);

class ContourRangeAndSplit : public ContourRangeFunctor
{
    ContourRangeAndSplit(const ContourRangeAndSplit&);
    void operator=(const ContourRangeAndSplit&);
  protected:
    ContourRangeAndSplit()
      {
      PolyOutputs = vtkThreadLocal<vtkPolyData>::New();
      }
    ~ContourRangeAndSplit()
      {
      PolyOutputs->Delete();
      }

    vtkThreadLocal<vtkDataObject>* outputs;
  public:
    vtkTypeMacro(ContourRangeAndSplit,ContourRangeFunctor);
    static ContourRangeAndSplit* New();
    virtual void PrintSelf(ostream& os, vtkIndent indent)
      {
      this->Superclass::PrintSelf(os,indent);
      }

    vtkThreadLocal<vtkPolyData>* PolyOutputs;

    void SetData( vtkDataSet* _input, vtkIncrementalPointLocator* _locator,
                  vtkIdType& _size, double* _values, int _number,
                  vtkDataArray* _scalars, int _compute,
                  vtkThreadLocal<vtkDataObject>* _outputs)
      {
      vtkCutter::GetCellTypeDimensions(cellTypeDimensions);

      this->input = _input;
      this->inPd = _input->GetPointData();
      this->inCd = _input->GetCellData();
      this->refLocator = _locator;
      this->estimatedSize = _size;
      this->values = _values;
      this->numContours = _number;
      this->inScalars = _scalars;
      this->computeScalars = _compute;

      this->outputs = _outputs;

      this->Init(this->MasterThreadId);
      }
  
    virtual void Init(int tid) const
      {
      this->Superclass::Init(tid);

      vtkPolyData* output = vtkPolyData::SafeDownCast(
          this->outputs->NewLocal(tid));
      this->PolyOutputs->SetLocal(output,tid);

      vtkPointData* pd = output->GetPointData();
      if (!this->computeScalars)
      {
        pd->CopyScalarsOff();
      }
      pd->InterpolateAllocate(
          this->inPd,this->estimatedSize,this->estimatedSize);
      this->outPd->SetLocal(pd,tid);

      vtkCellData* cd = output->GetCellData();
      cd->CopyAllocate(this->inCd,this->estimatedSize,this->estimatedSize);
      this->outCd->SetLocal(cd,tid);
      }
};

vtkStandardNewMacro(ContourRangeAndSplit);

/* ================================================================================
  Generic contouring: Functors for parallel execution with ScalarTree
 ================================================================================ */
class ContourTraversalFunctor : public vtkTreeFunctorInitializable
{
  ContourTraversalFunctor( const ContourTraversalFunctor& );
  void operator =( const ContourTraversalFunctor& );

protected:
  ContourTraversalFunctor()
    {
    Locator = vtkThreadLocal<vtkIncrementalPointLocator>::New();
    newPts = vtkThreadLocal<vtkPoints>::New();
    newVerts = vtkThreadLocal<vtkCellArray>::New();
    newLines = vtkThreadLocal<vtkCellArray>::New();
    newPolys = vtkThreadLocal<vtkCellArray>::New();
    outPd = vtkThreadLocal<vtkPointData>::New();
    outCd = vtkThreadLocal<vtkCellData>::New();
    Cells = vtkThreadLocal<vtkGenericCell>::New();
    CellsScalars = vtkThreadLocal<vtkDataArray>::New();
    }
  ~ContourTraversalFunctor()
    {
    Locator->Delete();
    newPts->Delete();
    newVerts->Delete();
    newLines->Delete();
    newPolys->Delete();
    outPd->Delete();
    outCd->Delete();
    Cells->Delete();
    CellsScalars->Delete();
    }

public:
  vtkTypeMacro(ContourTraversalFunctor,vtkTreeFunctorInitializable);
  static ContourTraversalFunctor* New();
  virtual void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  vtkThreadLocal<vtkIncrementalPointLocator>* Locator;
  vtkThreadLocal<vtkPoints>* newPts;
  vtkThreadLocal<vtkCellArray>* newVerts;
  vtkThreadLocal<vtkCellArray>* newLines;
  vtkThreadLocal<vtkCellArray>* newPolys;
  vtkThreadLocal<vtkPointData>* outPd;
  vtkThreadLocal<vtkCellData>* outCd;
  vtkThreadLocal<vtkGenericCell>* Cells;
  vtkThreadLocal<vtkDataArray>* CellsScalars;

  vtkDataArray* inScalars;
  vtkDataSet* input;
  vtkPointData* inPd;
  vtkCellData* inCd;
  vtkIdType estimatedSize;
  double* values;
  int numContours;
  int computeScalars;

  vtkIncrementalPointLocator* refLocator;
  vtkCellArray* outputVerts;
  vtkCellArray* outputLines;
  vtkCellArray* outputPolys;
  vtkCellData* outputCd;
  vtkPointData* outputPd;

  double ScalarValue;

  void SetData( vtkDataSet* _input, vtkPoints* _inPts, vtkCellData* _incd,
                  vtkPointData* _inpd, vtkIncrementalPointLocator* _locator,
                  vtkIdType& _size, double* _values, int _number,
                  vtkDataArray* _scalars, int _compute, vtkCellArray* _outputVerts,
                  vtkCellArray* _outputLines, vtkCellArray* _outputPolys,
                  vtkCellData* _outputCd, vtkPointData* _outputPd )
    {
    this->inCd = _incd;
    this->inPd = _inpd;
    this->input = _input;
    this->refLocator = _locator;
    this->estimatedSize = _size;
    this->values = _values;
    this->numContours = _number;
    this->inScalars = _scalars;
    this->computeScalars = _compute;

    outputVerts = _outputVerts;
    outputLines = _outputLines;
    outputPolys = _outputPolys;
    outputCd = _outputCd;
    outputPd = _outputPd;

    int tid = this->MasterThreadId;
    if ( !computeScalars )
      {
      _outputPd->CopyScalarsOff();
      }
    _outputPd->InterpolateAllocate( inPd, estimatedSize, estimatedSize );
    _outputCd->CopyAllocate( inCd, estimatedSize, estimatedSize );
    Locator->SetLocal( _locator, tid );
    newPts->SetLocal( _inPts, tid );
    newVerts->SetLocal( _outputVerts, tid );
    newLines->SetLocal( _outputLines, tid );
    newPolys->SetLocal( _outputPolys, tid );
    outPd->SetLocal( _outputPd, tid );
    outCd->SetLocal( _outputCd, tid );
    Cells->NewLocal(tid);
    vtkDataArray* cScalars = CellsScalars->NewLocal( inScalars, tid );
    cScalars->SetNumberOfComponents( inScalars->GetNumberOfComponents() );
    cScalars->Allocate( cScalars->GetNumberOfComponents() * VTK_CELL_SIZE );

    Initialized(tid);
    }

  vtkLocalData* getLocal(int tid) const
    {
    vtkSMPContourFilterLocalData* data =
      vtkSMPContourFilterLocalData::New();
    data->Locator = this->Locator->GetLocal(tid);
    data->newVerts = this->newVerts->GetLocal(tid);
    data->newLines = this->newLines->GetLocal(tid);
    data->newPolys = this->newPolys->GetLocal(tid);
    data->outPd = this->outPd->GetLocal(tid);
    data->outCd = this->outCd->GetLocal(tid);
    data->Cells = this->Cells->GetLocal(tid);
    data->CellsScalars = this->CellsScalars->GetLocal(tid);
    return data;
    }

  void operator()( vtkIdType id, vtkLocalData* d ) const
    {
    vtkSMPContourFilterLocalData* data =
      static_cast<vtkSMPContourFilterLocalData*>(d);
    vtkGenericCell* cell = data->Cells;
    vtkDataArray* scalars = data->CellsScalars;

    this->input->GetCell( id, cell );
    vtkIdList* cellPts = cell->GetPointIds();
    scalars->SetNumberOfTuples( cellPts->GetNumberOfIds() );
    this->inScalars->GetTuples( cellPts, scalars );

    cell->Contour( ScalarValue, scalars, data->Locator,
                   data->newVerts, data->newLines, data->newPolys,
                   this->inPd, data->outPd,
                   this->inCd, id, data->outCd);
    }

  void Init (int tid) const
    {
    vtkPoints* pts = this->newPts->NewLocal(tid);
    pts->Allocate( this->estimatedSize, this->estimatedSize );
    vtkIncrementalPointLocator* l = this->Locator->NewLocal( this->refLocator, tid );
    l->InitPointInsertion( pts, this->input->GetBounds(), this->estimatedSize );

    vtkCellArray* c = this->newVerts->NewLocal(tid);
    c->Allocate( this->estimatedSize, this->estimatedSize );
    c = this->newLines->NewLocal(tid);
    c->Allocate( this->estimatedSize, this->estimatedSize );
    c = this->newPolys->NewLocal(tid);
    c->Allocate( this->estimatedSize, this->estimatedSize );

    vtkPointData* pd = this->outPd->NewLocal(tid);
    if ( !this->computeScalars )
      {
      pd->CopyScalarsOff();
      }
    pd->InterpolateAllocate( this->inPd, this->estimatedSize, this->estimatedSize );

    vtkCellData* cd = this->outCd->NewLocal(tid);
    cd->CopyAllocate( this->inCd, this->estimatedSize, this->estimatedSize );

    this->Cells->NewLocal(tid);

    vtkDataArray* cScalars = this->CellsScalars->NewLocal( this->inScalars, tid );
    cScalars->SetNumberOfComponents( this->inScalars->GetNumberOfComponents() );
    cScalars->Allocate( cScalars->GetNumberOfComponents() * VTK_CELL_SIZE );

    Initialized(tid);
    }
};

vtkStandardNewMacro(ContourTraversalFunctor);

class ContourTraverseAndSplit : public ContourTraversalFunctor
{
    ContourTraverseAndSplit(const ContourTraverseAndSplit&);
    void operator=(const ContourTraverseAndSplit&);
  protected:
    ContourTraverseAndSplit()
      {
      PolyOutputs = vtkThreadLocal<vtkPolyData>::New();
      }
    ~ContourTraverseAndSplit()
      {
      PolyOutputs->Delete();
      }

    vtkThreadLocal<vtkDataObject>* outputs;
  public:
    vtkTypeMacro(ContourTraverseAndSplit,ContourTraversalFunctor);
    static ContourTraverseAndSplit* New();
    virtual void PrintSelf(ostream& os, vtkIndent indent)
      {
      this->Superclass::PrintSelf(os,indent);
      }

    vtkThreadLocal<vtkPolyData>* PolyOutputs;

    void SetData( vtkDataSet* _input, vtkIncrementalPointLocator* _locator,
                  vtkIdType& _size, double* _values, int _number,
                  vtkDataArray* _scalars, int _compute,
                  vtkThreadLocal<vtkDataObject>* _outputs)
      {
      this->input = _input;
      this->inPd = _input->GetPointData();
      this->inCd = _input->GetCellData();
      this->refLocator = _locator;
      this->estimatedSize = _size;
      this->values = _values;
      this->numContours = _number;
      this->inScalars = _scalars;
      this->computeScalars = _compute;

      this->outputs = _outputs;

      this->Init(this->MasterThreadId);
      }
  
    virtual void Init(int tid) const
      {
      this->Superclass::Init(tid);

      vtkPolyData* output = vtkPolyData::SafeDownCast(
          this->outputs->NewLocal(tid));
      this->PolyOutputs->SetLocal(output,tid);

      vtkPointData* pd = output->GetPointData();
      if (!this->computeScalars)
        {
        pd->CopyScalarsOff();
        }
      pd->InterpolateAllocate(
          this->inPd,this->estimatedSize,this->estimatedSize);
      this->outPd->SetLocal(pd,tid);

      vtkCellData* cd = output->GetCellData();
      cd->CopyAllocate(this->inCd,this->estimatedSize,this->estimatedSize);
      this->outCd->SetLocal(cd,tid);
      }
};

vtkStandardNewMacro(ContourTraverseAndSplit);



int vtkSMPContourFilter::ProcessRequest(vtkInformation* request,
                                        vtkInformationVector** inVector,
                                        vtkInformationVector* outVector)
  {
  if(!this->vtkSMPAlgorithm::ProcessRequest(
        request,inVector,outVector,this->GetNumberOfOutputPorts()))
    {
    return this->Superclass::ProcessRequest(request,inVector,outVector);
    }
  return 1;
  }

/* ================================================================================
 General contouring filter.  Handles arbitrary input.
 ================================================================================ */
int vtkSMPContourFilter::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
  {
  // get the input
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!input)
    {
    return 0;
    }

  // get the contours
  int numContours=this->ContourValues->GetNumberOfContours();
  double *values=this->ContourValues->GetValues();
  int i;

  // is there data to process?
  if (!this->GetInputArrayToProcess(0, inputVector))
    {
    return 1;
    }

  int sType = this->GetInputArrayToProcess(0, inputVector)->GetDataType();

  // handle 2D images
  if (vtkImageData::SafeDownCast(input) && sType != VTK_BIT &&
      !vtkUniformGrid::SafeDownCast(input))
    {
    int dim = 3;
    int *uExt = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
    if (uExt[0] == uExt[1])
      {
      --dim;
      }
    if (uExt[2] == uExt[3])
      {
      --dim;
      }
    if (uExt[4] == uExt[5])
      {
      --dim;
      }

    if ( dim == 2 )
      {
      this->SynchronizedTemplates2D->SetNumberOfContours(numContours);
      for (i=0; i < numContours; i++)
        {
        this->SynchronizedTemplates2D->SetValue(i,values[i]);
        }
      this->SynchronizedTemplates2D->
        SetInputArrayToProcess(0,this->GetInputArrayInformation(0));
      return
        this->SynchronizedTemplates2D->ProcessRequest(request,inputVector,outputVector);
      }
    else if ( dim == 3 )
      {
      this->SynchronizedTemplates3D->SetNumberOfContours(numContours);
      for (i=0; i < numContours; i++)
        {
        this->SynchronizedTemplates3D->SetValue(i,values[i]);
        }
      this->SynchronizedTemplates3D->SetComputeNormals(this->ComputeNormals);
      this->SynchronizedTemplates3D->SetComputeGradients(this->ComputeGradients);
      this->SynchronizedTemplates3D->SetComputeScalars(this->ComputeScalars);
      this->SynchronizedTemplates3D->
        SetInputArrayToProcess(0,this->GetInputArrayInformation(0));

      return this->SynchronizedTemplates3D->ProcessRequest(request,inputVector,outputVector);
      }
    } //if image data

  // handle 3D RGrids
  if (vtkRectilinearGrid::SafeDownCast(input) && sType != VTK_BIT)
    {
    int *uExt = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
    // if 3D
    if (uExt[0] < uExt[1] && uExt[2] < uExt[3] && uExt[4] < uExt[5])
      {
      this->RectilinearSynchronizedTemplates->SetNumberOfContours(numContours);
      for (i=0; i < numContours; i++)
        {
        this->RectilinearSynchronizedTemplates->SetValue(i,values[i]);
        }
      this->RectilinearSynchronizedTemplates->SetComputeNormals(this->ComputeNormals);
      this->RectilinearSynchronizedTemplates->SetComputeGradients(this->ComputeGradients);
      this->RectilinearSynchronizedTemplates->SetComputeScalars(this->ComputeScalars);
      this->RectilinearSynchronizedTemplates->
        SetInputArrayToProcess(0,this->GetInputArrayInformation(0));
      return this->RectilinearSynchronizedTemplates->
        ProcessRequest(request,inputVector,outputVector);
      }
    } // if 3D Rgrid

  // handle 3D SGrids
  if (vtkStructuredGrid::SafeDownCast(input) && sType != VTK_BIT)
    {
    int *uExt = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
    // if 3D
    if (uExt[0] < uExt[1] && uExt[2] < uExt[3] && uExt[4] < uExt[5])
      {
      this->GridSynchronizedTemplates->SetNumberOfContours(numContours);
      for (i=0; i < numContours; i++)
        {
        this->GridSynchronizedTemplates->SetValue(i,values[i]);
        }
      this->GridSynchronizedTemplates->SetComputeNormals(this->ComputeNormals);
      this->GridSynchronizedTemplates->SetComputeGradients(this->ComputeGradients);
      this->GridSynchronizedTemplates->SetComputeScalars(this->ComputeScalars);
      this->GridSynchronizedTemplates->
        SetInputArrayToProcess(0,this->GetInputArrayInformation(0));
      return this->GridSynchronizedTemplates->
        ProcessRequest(request,inputVector,outputVector);
      }
    } //if 3D SGrid

  vtkIdType cellId;
//  int abortExecute=0;
  vtkIdList *cellPts;
  vtkDataArray *inScalars;
  vtkCellArray *newVerts, *newLines, *newPolys;
  vtkPoints *newPts;
  vtkIdType numCells, estimatedSize;
  vtkDataArray *cellScalars;

  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    info->Get(vtkDataObject::DATA_OBJECT()));
  if (!output) {return 0;}


  vtkPointData *inPd=input->GetPointData(), *outPd=output->GetPointData();
  vtkCellData *inCd=input->GetCellData(), *outCd=output->GetCellData();

  vtkDebugMacro(<< "Executing contour filter");
  // Do not handle UnStructuredGrid since ForEach can't apply.
  // vtkContourGrid iterates over cells in a non-independant way

  numCells = input->GetNumberOfCells();
  inScalars = this->GetInputArrayToProcess(0,inputVector);
  if ( ! inScalars || numCells < 1 )
    {
    vtkDebugMacro(<<"No data to contour");
    return 1;
    }

  // Create objects to hold output of contour operation. First estimate
  // allocation size.
  //
  estimatedSize=
      static_cast<vtkIdType>(pow(static_cast<double>(numCells),.75));
  estimatedSize *= numContours;
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024)
    {
    estimatedSize = 1024;
    }

  newPts = vtkPoints::New();
  newPts->Allocate(estimatedSize,estimatedSize);
  newVerts = vtkCellArray::New();
  newLines = vtkCellArray::New();
  newPolys = vtkCellArray::New();
  cellScalars = inScalars->NewInstance();

  // locator used to merge potentially duplicate points
  if ( this->Locator == NULL )
    {
    this->CreateDefaultLocator();
    }
  this->Locator->InitPointInsertion (newPts, input->GetBounds(),estimatedSize);

  // interpolate data along edge
  // if we did not ask for scalars to be computed, don't copy them
  if (!this->ComputeScalars)
    {
    outPd->CopyScalarsOff();
    }

  // If enabled, build a scalar tree to accelerate search
  //
  vtkSMPMinMaxTree* parallelTree = vtkSMPMinMaxTree::SafeDownCast(this->ScalarTree);
  if ( !this->UseScalarTree || parallelTree )
    {
    vtkThreadLocal<vtkIncrementalPointLocator>* computedLocator;
    vtkThreadLocal<vtkPoints>* computedNewPts;
    vtkThreadLocal<vtkCellArray>* computedNewVerts;
    vtkThreadLocal<vtkCellArray>* computedNewLines;
    vtkThreadLocal<vtkCellArray>* computedNewPolys;
    vtkThreadLocal<vtkPointData>* computedOutPd;
    vtkThreadLocal<vtkCellData>* computedOutCd;

    // Init (thread local init is drown into first ForEach)
    input->GetCellType( 0 ); // Build cell representation so that Threads can access them safely
    // Exec
    if ( this->UseScalarTree )
      {
      ContourTraversalFunctor* TreeContour = ContourTraversalFunctor::New();
      TreeContour->SetData( input, newPts, inCd, inPd, this->Locator,
                            estimatedSize, values, numContours, inScalars, this->ComputeScalars,
                            newVerts, newLines, newPolys, outCd, outPd );
      parallelTree->SetDataSet(input);
      for ( i = 0; i < numContours; ++i )
        {
        TreeContour->ScalarValue = values[i];
        parallelTree->InitTraversal( values[i] );
        vtkParallelOperators::Traverse( parallelTree, TreeContour );
        }
      computedLocator = TreeContour->Locator;
      computedNewPts = TreeContour->newPts;
      computedNewVerts = TreeContour->newVerts;
      computedNewLines = TreeContour->newLines;
      computedNewPolys = TreeContour->newPolys;
      computedOutPd = TreeContour->outPd;
      computedOutCd = TreeContour->outCd;
      computedLocator->Register(this);
      computedNewPts->Register(this);
      computedNewVerts->Register(this);
      computedNewLines->Register(this);
      computedNewPolys->Register(this);
      computedOutPd->Register(this);
      computedOutCd->Register(this);
      TreeContour->Delete();
      }
    else
      {
      ContourRangeFunctor* contourFunctor = ContourRangeFunctor::New();
      contourFunctor->SetData( input, newPts, inCd, inPd, this->Locator,
                           estimatedSize, values, numContours, inScalars, this->ComputeScalars,
                           newVerts, newLines, newPolys, outCd, outPd );
      for ( contourFunctor->dimensionality = 1; contourFunctor->dimensionality <= 3; ++(contourFunctor->dimensionality) )
        {
        vtkParallelOperators::ForEach( 0, numCells, contourFunctor );
        }
      computedLocator = contourFunctor->Locator;
      computedNewPts = contourFunctor->newPts;
      computedNewVerts = contourFunctor->newVerts;
      computedNewLines = contourFunctor->newLines;
      computedNewPolys = contourFunctor->newPolys;
      computedOutPd = contourFunctor->outPd;
      computedOutCd = contourFunctor->outCd;
      computedLocator->Register(this);
      computedNewPts->Register(this);
      computedNewVerts->Register(this);
      computedNewLines->Register(this);
      computedNewPolys->Register(this);
      computedOutPd->Register(this);
      computedOutCd->Register(this);
      contourFunctor->Delete();
      }
    // Merge
    vtkMergeDataSets* mergeOp = vtkMergeDataSets::New();
    mergeOp->MasterThreadPopulatedOutputOn();
    vtkSMPMergePoints* parallelLocator = vtkSMPMergePoints::SafeDownCast( this->Locator );
    if ( parallelLocator )
      {
      vtkThreadLocal<vtkSMPMergePoints>* SMPLocator = vtkThreadLocal<vtkSMPMergePoints>::New();
      computedLocator->FillDerivedThreadLocal(SMPLocator);
      mergeOp->MergePolyData(
          parallelLocator, SMPLocator,
          outPd, computedOutPd,
          newVerts, computedNewVerts,
          newLines, computedNewLines,
          newPolys, computedNewPolys,
          0, 0, outCd, computedOutCd);
      SMPLocator->Delete();
      }
    else
      {
      mergeOp->MergePolyData(
          newPts, computedNewPts, input->GetBounds(),
          outPd, computedOutPd,
          newVerts, computedNewVerts,
          newLines, computedNewLines,
          newPolys, computedNewPolys,
          0, 0, outCd, computedOutCd);
      }
    mergeOp->Delete();
    computedLocator->UnRegister(this);
    computedNewPts->UnRegister(this);
    computedNewVerts->UnRegister(this);
    computedNewLines->UnRegister(this);
    computedNewPolys->UnRegister(this);
    computedOutPd->UnRegister(this);
    computedOutCd->UnRegister(this);
    } //if using scalar tree
  else
    {
    // Move of previously deleted Allocations
    newVerts->Allocate( estimatedSize, estimatedSize );
    newLines->Allocate( estimatedSize, estimatedSize );
    newPolys->Allocate( estimatedSize, estimatedSize );
    outPd->InterpolateAllocate( inPd, estimatedSize, estimatedSize );
    outCd->CopyAllocate( inCd, estimatedSize, estimatedSize );
    cellScalars->SetNumberOfComponents( inScalars->GetNumberOfComponents() );
    cellScalars->Allocate( cellScalars->GetNumberOfComponents() * VTK_CELL_SIZE );

    vtkCell *cell;
    if ( this->ScalarTree == NULL )
      {
      this->ScalarTree = vtkSimpleScalarTree::New();
      }
    this->ScalarTree->SetDataSet(input);
    // Note: This will have problems when input contains 2D and 3D cells.
    // CellData will get scrabled because of the implicit ordering of
    // verts, lines and polys in vtkPolyData.  The solution
    // is to convert this filter to create unstructured grid.
    //
    // Loop over all contour values.  Then for each contour value,
    // loop over all cells.
    //
    for ( i=0; i < numContours; i++ )
      {
      for ( this->ScalarTree->InitTraversal(values[i]);
            (cell=this->ScalarTree->GetNextCell(cellId,cellPts,cellScalars)) != NULL; )
        {
        cell->Contour(values[i], cellScalars, this->Locator,
                      newVerts, newLines, newPolys, inPd, outPd,
                      inCd, cellId, outCd);

        } //for all cells
      } //for all contour values
    } //using scalar tree

  vtkDebugMacro(<<"Created: "
                << newPts->GetNumberOfPoints() << " points, "
                << newVerts->GetNumberOfCells() << " verts, "
                << newLines->GetNumberOfCells() << " lines, "
                << newPolys->GetNumberOfCells() << " triangles");

  // Update ourselves.  Because we don't know up front how many verts, lines,
  // polys we've created, take care to reclaim memory.
  //
  output->SetPoints(newPts);
  newPts->Delete();
  cellScalars->Delete();

  if (newVerts->GetNumberOfCells())
    {
    output->SetVerts(newVerts);
    }
  newVerts->Delete();

  if (newLines->GetNumberOfCells())
    {
    output->SetLines(newLines);
    }
  newLines->Delete();

  if (newPolys->GetNumberOfCells())
    {
    output->SetPolys(newPolys);
    }
  newPolys->Delete();

  this->Locator->Initialize();//releases leftover memory
  output->Squeeze();

  return 1;
  }

int vtkSMPContourFilter::SplitData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector,
    vtkThreadLocal<vtkDataObject>** outputs)
  {
  // get the input
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!input)
    {
    return 0;
    }

  // get the contours
  int numContours=this->ContourValues->GetNumberOfContours();
  double *values=this->ContourValues->GetValues();
  int i;

  // is there data to process?
  if (!this->GetInputArrayToProcess(0, inputVector))
    {
    return 1;
    }

  int sType = this->GetInputArrayToProcess(0, inputVector)->GetDataType();

  // handle 2D images
  if (vtkImageData::SafeDownCast(input) && sType != VTK_BIT &&
      !vtkUniformGrid::SafeDownCast(input))
    {
    vtkWarningMacro(<< "vtkImageData input not supported");
    return 0;
    } //if image data

  // handle 3D RGrids
  if (vtkRectilinearGrid::SafeDownCast(input) && sType != VTK_BIT)
    {
    vtkWarningMacro(<< "vtkRectilinearGrid input not supported");
    return 0;
    } // if 3D Rgrid

  // handle 3D SGrids
  if (vtkStructuredGrid::SafeDownCast(input) && sType != VTK_BIT)
    {
    vtkWarningMacro(<< "vtkStructuredGrid input not supported");
    } //if 3D SGrid

  vtkDataArray *inScalars;
  vtkIdType numCells, estimatedSize;

  if (!outputs[0]) {return 0;}

  vtkDebugMacro(<< "Executing contour filter");
  if (input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID && !this->UseScalarTree)
    {
    vtkDebugMacro(<< "vtkUnstructuredGrid input not supported");
    return 0;
    } //if type VTK_UNSTRUCTURED_GRID
  else
    {
    numCells = input->GetNumberOfCells();
    inScalars = this->GetInputArrayToProcess(0,inputVector);
    if ( ! inScalars || numCells < 1 )
      {
      vtkDebugMacro(<<"No data to contour");
      return 1;
      }

    // Create objects to hold output of contour operation. First estimate
    // allocation size.
    //
    estimatedSize=
        static_cast<vtkIdType>(pow(static_cast<double>(numCells),.75));
    estimatedSize *= numContours;
    estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
    if (estimatedSize < 1024)
      {
      estimatedSize = 1024;
      }

    // locator used to merge potentially duplicate points
    if ( this->Locator == NULL )
      {
      this->CreateDefaultLocator();
      }

    // If enabled, build a scalar tree to accelerate search
    //
    vtkSMPMinMaxTree* parallelTree = vtkSMPMinMaxTree::SafeDownCast(this->ScalarTree);
    if ( !this->UseScalarTree || parallelTree )
      {
      vtkThreadLocal<vtkPolyData>* outDatasets;
      vtkThreadLocal<vtkPoints>* outPoints;
      vtkThreadLocal<vtkCellArray>* outVerts;
      vtkThreadLocal<vtkCellArray>* outLines;
      vtkThreadLocal<vtkCellArray>* outPolys;

      // Init (thread local init is drown into first ForEach)
      input->GetCellType( 0 ); // Build cell representation so that Threads can access them safely
      if ( this->UseScalarTree )
        {
        ContourTraverseAndSplit* treeFunctor = ContourTraverseAndSplit::New();
        treeFunctor->SetData( input, this->Locator, estimatedSize, values,
                           numContours, inScalars, this->ComputeScalars, outputs[0] );
        parallelTree->SetDataSet(input);
        for ( i = 0; i < numContours; ++i )
          {
          treeFunctor->ScalarValue = values[i];
          parallelTree->InitTraversal(values[i]);
          vtkParallelOperators::Traverse(parallelTree,treeFunctor);
          }
        outDatasets = treeFunctor->PolyOutputs;
        outPoints = treeFunctor->newPts;
        outVerts = treeFunctor->newVerts;
        outLines = treeFunctor->newLines;
        outPolys = treeFunctor->newPolys;
        outDatasets->Register(this);
        outPoints->Register(this);
        outVerts->Register(this);
        outLines->Register(this);
        outPolys->Register(this);
        treeFunctor->Delete();
        }
      else
        {
        ContourRangeAndSplit* contourFunctor = ContourRangeAndSplit::New();
        contourFunctor->SetData( input, this->Locator, estimatedSize, values,
                           numContours, inScalars, this->ComputeScalars, outputs[0] );
        for ( contourFunctor->dimensionality = 1; contourFunctor->dimensionality <= 3; ++(contourFunctor->dimensionality) )
          {
          vtkParallelOperators::ForEach( 0, numCells, contourFunctor );
          }
        outDatasets = contourFunctor->PolyOutputs;
        outPoints = contourFunctor->newPts;
        outVerts = contourFunctor->newVerts;
        outLines = contourFunctor->newLines;
        outPolys = contourFunctor->newPolys;
        outDatasets->Register(this);
        outPoints->Register(this);
        outVerts->Register(this);
        outLines->Register(this);
        outPolys->Register(this);
        contourFunctor->Delete();
        }

      // No more merge in V2 but we finalize each output
      // Easier in sequential
      vtkThreadLocal<vtkPolyData>::iterator itOutput;
      vtkThreadLocal<vtkPoints>::iterator newPts;
      vtkThreadLocal<vtkCellArray>::iterator newVerts;
      vtkThreadLocal<vtkCellArray>::iterator newLines;
      vtkThreadLocal<vtkCellArray>::iterator newPolys;
      for ( itOutput = outDatasets->Begin(),
            newPts = outPoints->Begin(),
            newVerts = outVerts->Begin(),
            newLines = outLines->Begin(),
            newPolys = outPolys->Begin();
            itOutput != outDatasets->End();
            ++itOutput, ++newPts, ++newVerts, ++newLines, ++newPolys)
        {
        vtkDebugMacro(<<"Created: "
                      << (*newPts)->GetNumberOfPoints() << " points, "
                      << (*newVerts)->GetNumberOfCells() << " verts, "
                      << (*newLines)->GetNumberOfCells() << " lines, "
                      << (*newPolys)->GetNumberOfCells() << " triangles");
        (*itOutput)->SetPoints(*newPts);

        if ((*newVerts)->GetNumberOfCells())
          {
          (*itOutput)->SetVerts(*newVerts);
          }

        if ((*newLines)->GetNumberOfCells())
          {
          (*itOutput)->SetLines(*newLines);
          }

        if ((*newPolys)->GetNumberOfCells())
          {
          (*itOutput)->SetPolys(*newPolys);
          }
        (*itOutput)->Squeeze();
        }

      outDatasets->UnRegister(this);
      outPoints->UnRegister(this);
      outVerts->UnRegister(this);
      outLines->UnRegister(this);
      outPolys->UnRegister(this);
      } //if using scalar tree
    else
      {
      vtkWarningMacro(<< "Scalar trees that are not a vtkSMPMinMaxTree not supported");
      return 0;
      } //using scalar tree

    this->Locator->Initialize();//releases leftover memory
    } //else if not vtkUnstructuredGrid

  return 1;
  }
