#include "vtkSMPClipDataSet.h"
#include "vtkParallelOperators.h"
#include "vtkRangeFunctor.h"
#include "vtkRangeFunctorInitializable.h"
#include "vtkRange1D.h"
#include "vtkMergeDataSets.h"
#include "vtkThreadLocal.h"

#include "vtkCallbackCommand.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkClipVolume.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkImageData.h"
#include "vtkImplicitFunction.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkSMPMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkPolyhedron.h"

#include <math.h>

class GenerateClipValueExecutor : public vtkRangeFunctor
{
  GenerateClipValueExecutor(const GenerateClipValueExecutor&);
  void operator =(const GenerateClipValueExecutor&);

protected:
  GenerateClipValueExecutor(){}
  ~GenerateClipValueExecutor(){}

  vtkFloatArray* scalars;
  vtkDataSet* input;
  vtkImplicitFunction* clipFunction;

public:
  vtkTypeMacro(GenerateClipValueExecutor,vtkRangeFunctor);
  static GenerateClipValueExecutor* New();
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  void SetData(vtkFloatArray* s, vtkDataSet* i, vtkImplicitFunction* f)
    {
    scalars = s;
    input = i;
    clipFunction = f;
    }

  void operator()(vtkRange* r) const
    {
    vtkRange1D* range = vtkRange1D::SafeDownCast(r);
    double p[3];
    for (vtkIdType i = range->Begin(); i < range->End(); ++i)
      {
      this->input->GetPoint(i,p);
      double s = this->clipFunction->FunctionValue(p);
      this->scalars->SetTuple1(i,s);
      }
    }
};

class ClipDataSetRange : public vtkRangeFunctorInitializable
{
  ClipDataSetRange(const ClipDataSetRange&);
  void operator =(const ClipDataSetRange&);

protected:
  ClipDataSetRange() : input(0), clipScalars(0), inPD(0), inCD(0), estimatedSize(0),
             value(0.0), insideOut(0), numOutputs(0), refLocator(0)
    {
    cell = vtkThreadLocal<vtkGenericCell>::New();
    cellScalars = vtkThreadLocal<vtkFloatArray>::New();
    newPoints = vtkThreadLocal<vtkPoints>::New();
    outPD = vtkThreadLocal<vtkPointData>::New();
    locator = vtkThreadLocal<vtkIncrementalPointLocator>::New();

    conn[0] = conn[1] = 0;
    outCD[0] = outCD[1] = 0;
    types[0] = types[1] = 0;
    locs[0] = locs[1] = 0;
    }
  ~ClipDataSetRange()
    {
    cell->Delete();
    cellScalars->Delete();
    newPoints->Delete();
    outPD->Delete();
    locator->Delete();
    for (int i = 0; i < numOutputs; ++i)
      {
      conn[i]->Delete();
      outCD[i]->Delete();
      types[i]->Delete();
      locs[i]->Delete();
      }
    }

  vtkDataSet* input;
  vtkDataArray* clipScalars;
  vtkPointData* inPD;
  vtkCellData* inCD;
  double value;
  int insideOut;
  int numOutputs;
  vtkIdType estimatedSize;
  vtkIncrementalPointLocator* refLocator;

public:
  vtkThreadLocal<vtkGenericCell>* cell;
  vtkThreadLocal<vtkFloatArray>* cellScalars;
  vtkThreadLocal<vtkPoints>* newPoints;
  vtkThreadLocal<vtkPointData>* outPD;
  vtkThreadLocal<vtkIncrementalPointLocator>* locator;

  vtkThreadLocal<vtkCellArray>* conn[2];
  vtkThreadLocal<vtkCellData>* outCD[2];
  vtkThreadLocal<vtkUnsignedCharArray>* types[2];
  vtkThreadLocal<vtkIdTypeArray>* locs[2];

  vtkTypeMacro(ClipDataSetRange,vtkRangeFunctorInitializable);
  static ClipDataSetRange* New();
  virtual void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  void SetData(vtkDataSet* i, vtkDataArray* c, vtkPointData* pd, vtkCellData* cd,
               vtkIdType e, double v, int inside, int num, vtkPoints* pts,
               vtkPointData* out_pd, vtkIncrementalPointLocator* l, vtkCellArray** _conn,
               vtkCellData** _out, vtkUnsignedCharArray** _types, vtkIdTypeArray** _locs)
    {
    input = i; clipScalars = c; inPD = pd; inCD = cd; refLocator = l;
    estimatedSize = e; value = v; insideOut = inside; numOutputs = num;

    int tid = this->MasterThreadId;
    input->GetCell(0,cell->NewLocal(tid));
    cellScalars->NewLocal(tid)->Allocate(VTK_CELL_SIZE);
    newPoints->SetLocal(pts,tid);
    locator->SetLocal(l,tid);
    outPD->SetLocal(out_pd,tid);
    for (int i = 0; i < num; ++i)
      {
      conn[i] = vtkThreadLocal<vtkCellArray>::New();
      conn[i]->SetLocal(_conn[i],tid);
      outCD[i] = vtkThreadLocal<vtkCellData>::New();
      outCD[i]->SetLocal(_out[i],tid);
      types[i] = vtkThreadLocal<vtkUnsignedCharArray>::New();
      types[i]->SetLocal(_types[i],tid);
      locs[i] = vtkThreadLocal<vtkIdTypeArray>::New();
      locs[i]->SetLocal(_locs[i],tid);
      }

    Initialized(tid);
    }

  virtual void Init(int tid) const
    {
    cell->NewLocal(tid);
    cellScalars->NewLocal(tid)->Allocate(VTK_CELL_SIZE);
    vtkPoints* pts = newPoints->NewLocal(tid);
    pts->Allocate(estimatedSize,estimatedSize/2);
    outPD->NewLocal(tid)->InterpolateAllocate(inPD,estimatedSize,estimatedSize/2);
    locator->NewLocal(refLocator,tid)->InitPointInsertion (pts, input->GetBounds());
    for (int i = 0; i < numOutputs; ++i)
      {
      vtkCellArray* _conn = conn[i]->NewLocal(tid);
      _conn->Allocate(estimatedSize,estimatedSize/2);
      _conn->InitTraversal();
      outCD[i]->NewLocal(tid)->CopyAllocate(inCD,estimatedSize,estimatedSize/2);
      types[i]->NewLocal(tid)->Allocate(estimatedSize,estimatedSize/2);
      locs[i]->NewLocal(tid)->Allocate(estimatedSize,estimatedSize/2);
      }
    Initialized(tid);
    }

  void operator()(vtkRange* r) const
    {
    vtkRange1D* range = vtkRange1D::SafeDownCast(r);
    int tid = range->GetTid();
    vtkGenericCell* cell = this->cell->GetLocal(tid);
    vtkFloatArray* cellScalars = this->cellScalars->GetLocal(tid);
    vtkPointData* outPD = this->outPD->GetLocal(tid);
    vtkIncrementalPointLocator* locator = this->locator->GetLocal(tid);
    vtkCellArray* conn[2] = {this->conn[0]->GetLocal(tid),0};
    vtkCellData* outCD[2] = {this->outCD[0]->GetLocal(tid),0};
    vtkUnsignedCharArray* types[2] = {this->types[0]->GetLocal(tid),0};
    vtkIdTypeArray* locs[2] = {this->locs[0]->GetLocal(tid),0};
    if ( numOutputs > 1 )
      {
      conn[1] = this->conn[1]->GetLocal(tid);
      outCD[1] = this->outCD[1]->GetLocal(tid);
      types[1] = this->types[1]->GetLocal(tid);
      locs[1] = this->locs[1]->GetLocal(tid);
      }
    
    for (vtkIdType cellId = range->Begin(); cellId < range->End(); ++cellId)
      {
      input->GetCell(cellId,cell);
      vtkPoints* cellPts = cell->GetPoints();
      vtkIdList* cellIds = cell->GetPointIds();
      vtkIdType npts = cellPts->GetNumberOfPoints();

      // evaluate implicit cutting function
      for (int i=0; i < npts; i++ )
        {
        double s = clipScalars->GetComponent(cellIds->GetId(i), 0);
        cellScalars->InsertTuple(i, &s);
        }

      // perform the clipping
      int num[2];
      num[0] = conn[0]->GetNumberOfCells();
      cell->Clip(value, cellScalars, locator,
          conn[0], inPD, outPD, inCD, cellId,
          outCD[0], this->insideOut);
      num[0] = conn[0]->GetNumberOfCells() - num[0];

      if ( numOutputs > 1 )
        {
        num[1] = conn[1]->GetNumberOfCells();
        cell->Clip(value, cellScalars, locator,
            conn[1], inPD, outPD, inCD, cellId,
            outCD[1], this->insideOut);
        num[1] = conn[1]->GetNumberOfCells() - num[0];
        }

      for (int i=0; i<numOutputs; i++) //for both outputs
        {
        for (int j=0; j < num[i]; j++)
          {
          if (cell->GetCellType() == VTK_POLYHEDRON)
            {
            //Polyhedron cells have a special cell connectivity format
            //(nCell0Faces, nFace0Pts, i, j, k, nFace1Pts, i, j, k, ...).
            //But we don't need to deal with it here. The special case is handled
            //by vtkUnstructuredGrid::SetCells(), which will be called next.
            types[i]->InsertNextValue(VTK_POLYHEDRON);
            }
          else
            {
            locs[i]->InsertNextValue(conn[i]->GetTraversalLocation());
            vtkIdType* pts;
            conn[i]->GetNextCell(npts,pts);

            int cellType = 0;
            //For each new cell added, got to set the type of the cell
            switch ( cell->GetCellDimension() )
              {
              case 0: //points are generated--------------------------------
                cellType = (npts > 1 ? VTK_POLY_VERTEX : VTK_VERTEX);
                break;

              case 1: //lines are generated---------------------------------
                cellType = (npts > 2 ? VTK_POLY_LINE : VTK_LINE);
                break;

              case 2: //polygons are generated------------------------------
                cellType = (npts == 3 ? VTK_TRIANGLE :
                                        (npts == 4 ? VTK_QUAD : VTK_POLYGON));
                break;

              case 3: //tetrahedra or wedges are generated------------------
                cellType = (npts == 4 ? VTK_TETRA : VTK_WEDGE);
                break;
              } //switch

            types[i]->InsertNextValue(cellType);
            }
          } //for each new cell
        } //for both outputs
      }
    }
};

class ClipDataSetSplit : public ClipDataSetRange
{
    ClipDataSetSplit(const ClipDataSetSplit&);
    void operator=(const ClipDataSetSplit&);
  protected:
    ClipDataSetSplit() : ClipDataSetRange(), outputs(0), generatedOutputs(0)
      {
      USOutputs = vtkThreadLocal<vtkUnstructuredGrid>::New();
      generatedUSOutputs = vtkThreadLocal<vtkUnstructuredGrid>::New();
      }
    ~ClipDataSetSplit()
      {
      USOutputs->Delete();
      generatedUSOutputs->Delete();
      }

    vtkThreadLocal<vtkDataObject>* outputs;
    vtkThreadLocal<vtkDataObject>* generatedOutputs;

    vtkThreadLocal<vtkUnstructuredGrid>* USOutputs;
    vtkThreadLocal<vtkUnstructuredGrid>* generatedUSOutputs;

  public:
    vtkTypeMacro(ClipDataSetSplit,ClipDataSetRange);
    static ClipDataSetSplit* New();
    virtual void PrintSelf(ostream& os, vtkIndent indent)
      {
      this->Superclass::PrintSelf(os,indent);
      }

    void SetData(vtkDataSet* i, vtkDataArray* c, vtkPointData* pd, vtkCellData* cd,
                 vtkIdType e, double v, int inside, int num,
                 vtkIncrementalPointLocator* l, vtkThreadLocal<vtkDataObject>* o,
                 vtkThreadLocal<vtkDataObject>* go)
      {
      input = i; clipScalars = c; inPD = pd; inCD = cd; refLocator = l;
      estimatedSize = e; value = v; insideOut = inside; numOutputs = num;

      this->outputs = o;
      this->generatedOutputs = go;

      this->Init(this->MasterThreadId);
      }

    virtual void Init(int tid) const
      {
      vtkUnstructuredGrid* output[2] = { vtkUnstructuredGrid::SafeDownCast(
          this->outputs->GetLocal(tid)), NULL };
      USOutputs->SetLocal(output[0],tid);
      if (numOutputs > 1)
        {
        output[1] = vtkUnstructuredGrid::SafeDownCast(
          this->generatedOutputs->GetLocal(tid));
        generatedUSOutputs->SetLocal(output[1],tid);
        }
      cell->NewLocal(tid);
      cellScalars->NewLocal(tid)->Allocate(VTK_CELL_SIZE);
      vtkPoints* pts = newPoints->NewLocal(tid);
      pts->Allocate(estimatedSize,estimatedSize/2);
      vtkPointData* pd = output[0]->GetPointData();
      vtkDataSetAttributes* tempDSA = vtkDataSetAttributes::New();
      tempDSA->InterpolateAllocate(inPD, 1, 2);
      pd->InterpolateAllocate(inPD,estimatedSize,estimatedSize/2);
      tempDSA->Delete();
      outPD->SetLocal(pd,tid);
      locator->NewLocal(refLocator,tid)->InitPointInsertion (pts, input->GetBounds());
      vtkCellData* cd;
      for (int i = 0; i < numOutputs; ++i)
        {
        vtkCellArray* _conn = conn[i]->NewLocal(tid);
        _conn->Allocate(estimatedSize,estimatedSize/2);
        _conn->InitTraversal();
        cd = output[i]->GetCellData();
        cd->CopyAllocate(inCD,estimatedSize,estimatedSize/2);
        outCD[i]->SetLocal(cd,tid);
        types[i]->NewLocal(tid)->Allocate(estimatedSize,estimatedSize/2);
        locs[i]->NewLocal(tid)->Allocate(estimatedSize,estimatedSize/2);
        }
      Initialized(tid);
      }

    void FinalizeOutput()
      {
      vtkThreadLocal<vtkUnstructuredGrid>::iterator itOutput;
      vtkThreadLocal<vtkPoints>::iterator itPts;
      vtkThreadLocal<vtkUnsignedCharArray>::iterator itTypes;
      vtkThreadLocal<vtkIdTypeArray>::iterator itLocs;
      vtkThreadLocal<vtkCellArray>::iterator itConn;
      for ( itOutput = this->USOutputs->Begin(),
            itPts = this->newPoints->Begin(),
            itTypes = this->types[0]->Begin(),
            itLocs = this->locs[0]->Begin(),
            itConn = this->conn[0]->Begin();
            itOutput != this->USOutputs->End();
            ++itOutput, ++itPts, ++itTypes, ++itLocs, ++itConn )
        {
        (*itOutput)->SetPoints(*itPts);
        (*itOutput)->SetCells(*itTypes, *itLocs, *itConn);
        (*itOutput)->Squeeze();
        }
      if (numOutputs > 1)
        {
        for ( itOutput = this->generatedUSOutputs->Begin(),
              itPts = this->newPoints->Begin(),
              itTypes = this->types[1]->Begin(),
              itLocs = this->locs[1]->Begin(),
              itConn = this->conn[1]->Begin();
              itOutput != this->generatedUSOutputs->End();
              ++itOutput, ++itPts, ++itTypes, ++itLocs, ++itConn )
          {
          (*itOutput)->SetPoints(*itPts);
          (*itOutput)->SetCells(*itTypes, *itLocs, *itConn);
          (*itOutput)->Squeeze();
          }
        }
      }
};

vtkStandardNewMacro(ClipDataSetRange);
vtkStandardNewMacro(ClipDataSetSplit);
vtkStandardNewMacro(GenerateClipValueExecutor);
vtkStandardNewMacro(vtkSMPClipDataSet);

vtkSMPClipDataSet::vtkSMPClipDataSet() : vtkClipDataSet() { }
vtkSMPClipDataSet::~vtkSMPClipDataSet() { }

void vtkSMPClipDataSet::PrintSelf(ostream& os, vtkIndent indent)
  {
  this->Superclass::PrintSelf(os,indent);
  }

int vtkSMPClipDataSet::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
  {
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *realInput = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  // We have to create a copy of the input because clip requires being
  // able to InterpolateAllocate point data from the input that is
  // exactly the same as output. If the input arrays and output arrays
  // are different vtkCell3D's Clip will fail. By calling InterpolateAllocate
  // here, we make sure that the output will look exactly like the output
  // (unwanted arrays will be eliminated in InterpolateAllocate). The
  // last argument of InterpolateAllocate makes sure that arrays are shallow
  // copied from realInput to input.
  vtkSmartPointer<vtkDataSet> input;
  input.TakeReference(realInput->NewInstance());
  input->CopyStructure(realInput);
  input->GetCellData()->PassData(realInput->GetCellData());
  input->GetPointData()->InterpolateAllocate(realInput->GetPointData(), 0, 0, 1);

  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkUnstructuredGrid *clippedOutput = this->GetClippedOutput();

  vtkIdType numPts = input->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();
  vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
  vtkCellData *inCD=input->GetCellData();
  vtkCellData *outCD[2];
  vtkPoints *newPoints;
  vtkDataArray *clipScalars;
  vtkPoints *cellPts;
  vtkIdList *cellIds;
  double s;
  vtkIdType npts;
  vtkIdType *pts;
  int cellType = 0;
  vtkIdType i;
  int j;
  vtkIdType estimatedSize;
  vtkUnsignedCharArray *types[2];
  types[0] = types[1] = 0;
  vtkIdTypeArray *locs[2];
  locs[0] = locs[1] = 0;
  int numOutputs = 1;

  outCD[0] = 0;
  outCD[1] = 0;

  vtkDebugMacro(<< "Clipping dataset");

  int inputObjectType = input->GetDataObjectType();

  // if we have volumes
  if (inputObjectType == VTK_STRUCTURED_POINTS ||
      inputObjectType == VTK_IMAGE_DATA)
    {
    int dimension;
    int *dims = vtkImageData::SafeDownCast(input)->GetDimensions();
    for (dimension=3, i=0; i<3; i++)
      {
      if ( dims[i] <= 1 )
        {
        dimension--;
        }
      }
    if ( dimension >= 3 )
      {
      this->ClipVolume(input, output);
      return 1;
      }
     }

  // Initialize self; create output objects
  //
  if ( numPts < 1 )
    {
    vtkDebugMacro(<<"No data to clip");
    return 1;
    }

  if ( !this->ClipFunction && this->GenerateClipScalars )
    {
    vtkErrorMacro(<<"Cannot generate clip scalars if no clip function defined");
    return 1;
    }

  if ( numCells < 1 )
    {
    return this->ClipPoints(input, output, inputVector);
    }

  // allocate the output and associated helper classes
  estimatedSize = numCells;
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024)
    {
    estimatedSize = 1024;
    }
  vtkCellArray *conn[2];
  conn[0] = conn[1] = 0;
  conn[0] = vtkCellArray::New();
  conn[0]->Allocate(estimatedSize,estimatedSize/2);
  conn[0]->InitTraversal();
  types[0] = vtkUnsignedCharArray::New();
  types[0]->Allocate(estimatedSize,estimatedSize/2);
  locs[0] = vtkIdTypeArray::New();
  locs[0]->Allocate(estimatedSize,estimatedSize/2);
  if ( this->GenerateClippedOutput )
    {
    numOutputs = 2;
    conn[1] = vtkCellArray::New();
    conn[1]->Allocate(estimatedSize,estimatedSize/2);
    conn[1]->InitTraversal();
    types[1] = vtkUnsignedCharArray::New();
    types[1]->Allocate(estimatedSize,estimatedSize/2);
    locs[1] = vtkIdTypeArray::New();
    locs[1]->Allocate(estimatedSize,estimatedSize/2);
    }
  newPoints = vtkPoints::New();
  newPoints->Allocate(numPts,numPts/2);

  // locator used to merge potentially duplicate points
  if ( this->Locator == NULL )
    {
    this->CreateDefaultLocator();
    }
  this->Locator->InitPointInsertion (newPoints, input->GetBounds());

  // Determine whether we're clipping with input scalars or a clip function
  // and do necessary setup.
  if ( this->ClipFunction )
    {
    vtkFloatArray *tmpScalars = vtkFloatArray::New();
    tmpScalars->SetNumberOfTuples(numPts);
    tmpScalars->SetName("ClipDataSetScalars");
    inPD = vtkPointData::New();
    inPD->ShallowCopy(input->GetPointData());//copies original
    if ( this->GenerateClipScalars )
      {
      inPD->SetScalars(tmpScalars);
      }
    GenerateClipValueExecutor* generateClipScalarFunctor =
        GenerateClipValueExecutor::New();
    generateClipScalarFunctor->SetData(tmpScalars,input,this->ClipFunction);
    vtkParallelOperators::ForEach(0,numPts,generateClipScalarFunctor);
    generateClipScalarFunctor->Delete();
    clipScalars = tmpScalars;
    }
  else //using input scalars
    {
    clipScalars = this->GetInputArrayToProcess(0,inputVector);
    if ( !clipScalars )
      {
      for ( i=0; i<2; i++ )
        {
        if (conn[i])
          {
          conn[i]->Delete();
          }
        if (types[i])
          {
          types[i]->Delete();
          }
        if (locs[i])
          {
          locs[i]->Delete();
          }
        }
      newPoints->Delete();
      // When processing composite datasets with partial arrays, this warning is
      // not applicable, hence disabling it.
      // vtkErrorMacro(<<"Cannot clip without clip function or input scalars");
      return 1;
      }
    }

  // Refer to BUG #8494 and BUG #11016. I cannot see any reason why one would
  // want to turn CopyScalars Off. My understanding is that this was done to
  // avoid copying of "ClipDataSetScalars" to the output when
  // this->GenerateClipScalars is false. But, if GenerateClipScalars is false,
  // then "ClipDataSetScalars" is not added as scalars to the input at all
  // (refer to code above) so it's a non-issue. Leaving CopyScalars untouched
  // i.e. ON avoids dropping of arrays (#8484) as well as segfaults (#11016).
  //if ( !this->GenerateClipScalars &&
  //  !this->GetInputArrayToProcess(0,inputVector))
  //  {
  //  outPD->CopyScalarsOff();
  //  }
  //else
  //  {
  //  outPD->CopyScalarsOn();
  //  }
  vtkDataSetAttributes* tempDSA = vtkDataSetAttributes::New();
  tempDSA->InterpolateAllocate(inPD, 1, 2);
  outPD->InterpolateAllocate(inPD,estimatedSize,estimatedSize/2);
  tempDSA->Delete();
  outCD[0] = output->GetCellData();
  outCD[0]->CopyAllocate(inCD,estimatedSize,estimatedSize/2);
  if ( this->GenerateClippedOutput )
    {
    outCD[1] = clippedOutput->GetCellData();
    outCD[1]->CopyAllocate(inCD,estimatedSize,estimatedSize/2);
    }

  //Process all cells and clip each in turn
  //
  ClipDataSetRange* functor = ClipDataSetRange::New();
  functor->SetData(input, clipScalars, inPD, inCD, estimatedSize,
                   this->UseValueAsOffset || !this->ClipFunction ? this->Value : 0.0,
                   this->InsideOut, numOutputs, newPoints, outPD, this->Locator,
                   conn, outCD, types, locs);
  vtkParallelOperators::ForEach(0,numCells,functor);

  vtkSMPMergePoints* parallelLocator = vtkSMPMergePoints::SafeDownCast(this->Locator);
  vtkMergeDataSets* mergeOp = vtkMergeDataSets::New();
  mergeOp->MasterThreadPopulatedOutputOn();
  if (parallelLocator)
    {
    vtkThreadLocal<vtkSMPMergePoints>* partialMeshes;
    functor->locator->FillDerivedThreadLocal(partialMeshes);
    mergeOp->MergeUnstructuredGrid(
        parallelLocator, partialMeshes,
        outPD, functor->outPD,
        conn[0], functor->conn[0],
        outCD[0], functor->outCD[0],
        locs[0], functor->locs[0],
        types[0], functor->types[0]);
    partialMeshes->Delete();
    }
  else
    {
    mergeOp->MergeUnstructuredGrid(
        newPoints, functor->newPoints, input->GetBounds(),
        outPD, functor->outPD,
        conn[0], functor->conn[0],
        outCD[0], functor->outCD[0],
        locs[0], functor->locs[0],
        types[0], functor->types[0]);
    }

  if ( this->ClipFunction )
    {
    clipScalars->Delete();
    inPD->Delete();
    }

  output->SetPoints(newPoints);
  output->SetCells(types[0], locs[0], conn[0]);
  conn[0]->Delete();
  types[0]->Delete();
  locs[0]->Delete();

  if ( this->GenerateClippedOutput )
    {
    mergeOp->MergeUnstructuredGrid(
        0,0,0,0,
        conn[1], functor->conn[1],
        outCD[1], functor->outCD[1],
        locs[1], functor->locs[1],
        types[1], functor->types[1]);
    clippedOutput->SetPoints(newPoints);
    clippedOutput->SetCells(types[1], locs[1], conn[1]);
    conn[1]->Delete();
    types[1]->Delete();
    locs[1]->Delete();
    }

  mergeOp->Delete();
  functor->Delete();

  newPoints->Delete();
  this->Locator->Initialize();//release any extra memory
  output->Squeeze();

  return 1;
  }

vtkSplittingFilterStandardProcessRequest(vtkSMPClipDataSet);

int vtkSMPClipDataSet::SplitData(vtkInformation* request,
    vtkInformationVector** inVector,
    vtkInformationVector* outVector,
    vtkThreadLocal<vtkDataObject>** outputs)
  {
  // get the info objects
  vtkInformation *inInfo = inVector[0]->GetInformationObject(0);

  // get the input and output
  vtkDataSet *realInput = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  // We have to create a copy of the input because clip requires being
  // able to InterpolateAllocate point data from the input that is
  // exactly the same as output. If the input arrays and output arrays
  // are different vtkCell3D's Clip will fail. By calling InterpolateAllocate
  // here, we make sure that the output will look exactly like the output
  // (unwanted arrays will be eliminated in InterpolateAllocate). The
  // last argument of InterpolateAllocate makes sure that arrays are shallow
  // copied from realInput to input.
  vtkSmartPointer<vtkDataSet> input;
  input.TakeReference(realInput->NewInstance());
  input->CopyStructure(realInput);
  input->GetCellData()->PassData(realInput->GetCellData());
  input->GetPointData()->InterpolateAllocate(realInput->GetPointData(), 0, 0, 1);

  vtkIdType numPts = input->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();
  vtkPointData *inPD=input->GetPointData();
  vtkCellData *inCD=input->GetCellData();
  vtkDataArray *clipScalars;
  vtkPoints *cellPts;
  vtkIdList *cellIds;
  double s;
  vtkIdType npts;
  vtkIdType *pts;
  int cellType = 0;
  vtkIdType i;
  int j;
  vtkIdType estimatedSize;
  int numOutputs = 1;
  vtkThreadLocal<vtkDataObject>* gOutputs = NULL;

  vtkDebugMacro(<< "Clipping dataset");

  int inputObjectType = input->GetDataObjectType();

  // if we have volumes
  if (inputObjectType == VTK_STRUCTURED_POINTS)
    {
    vtkErrorMacro(<< "vtkStructuredPoints not supported");
    return 1;
    }
  if (inputObjectType == VTK_IMAGE_DATA)
    {
    vtkErrorMacro(<< "vtkImageData not supported");
    return 1;
    }

  // Initialize self; create output objects
  //
  if ( numPts < 1 )
    {
    vtkDebugMacro(<<"No data to clip");
    return 1;
    }

  if ( !this->ClipFunction && this->GenerateClipScalars )
    {
    vtkErrorMacro(<<"Cannot generate clip scalars if no clip function defined");
    return 1;
    }

  if ( numCells < 1 )
    {
    vtkErrorMacro(<< "Points without cells are not supported");
    return 1;
    }

  // allocate the output and associated helper classes
  estimatedSize = numCells;
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024)
    {
    estimatedSize = 1024;
    }
  if ( this->GenerateClippedOutput )
    {
    numOutputs = 2;
    gOutputs = outputs[1];
    }

  // locator used to merge potentially duplicate points
  if ( this->Locator == NULL )
    {
    this->CreateDefaultLocator();
    }

  // Determine whether we're clipping with input scalars or a clip function
  // and do necessary setup.
  if ( this->ClipFunction )
    {
    vtkFloatArray *tmpScalars = vtkFloatArray::New();
    tmpScalars->SetNumberOfTuples(numPts);
    tmpScalars->SetName("ClipDataSetScalars");
    inPD = vtkPointData::New();
    inPD->ShallowCopy(input->GetPointData());//copies original
    if ( this->GenerateClipScalars )
      {
      inPD->SetScalars(tmpScalars);
      }
    GenerateClipValueExecutor* generateClipScalarFunctor =
        GenerateClipValueExecutor::New();
    generateClipScalarFunctor->SetData(tmpScalars,input,this->ClipFunction);
    vtkParallelOperators::ForEach(0,numPts,generateClipScalarFunctor);
    generateClipScalarFunctor->Delete();
    clipScalars = tmpScalars;
    }
  else //using input scalars
    {
    clipScalars = this->GetInputArrayToProcess(0,inVector);
    if ( !clipScalars )
      {
      // When processing composite datasets with partial arrays, this warning is
      // not applicable, hence disabling it.
      // vtkErrorMacro(<<"Cannot clip without clip function or input scalars");
      return 1;
      }
    }

  ClipDataSetSplit* functor = ClipDataSetSplit::New();
  functor->SetData(input, clipScalars, inPD, inCD, estimatedSize,
                   this->UseValueAsOffset || !this->ClipFunction ? this->Value : 0.0,
                   this->InsideOut, numOutputs, this->Locator, outputs[0], gOutputs);
  vtkParallelOperators::ForEach(0,numCells,functor);

  // No merge but we finalize each output
  functor->FinalizeOutput();

  functor->Delete();
  return 1;
  }
