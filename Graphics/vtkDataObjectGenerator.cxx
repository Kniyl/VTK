/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDataObjectGenerator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkDataObjectGenerator.h"

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkDataObjectTypes.h>
#include <vtkImageData.h>
#include <vtkUniformGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkHierarchicalDataSet.h>
#include <vtkHierarchicalBoxDataSet.h>
#include <vtkAMRBox.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <vtkStreamingDemandDrivenPipeline.h>

#include <vtkstd/vector>

vtkCxxRevisionMacro(vtkDataObjectGenerator, "1.1");
vtkStandardNewMacro(vtkDataObjectGenerator);

//============================================================================
enum vtkDataObjectGeneratorTokenCodes
  {
  ID1, UF1, RG1, SG1, PD1, US1, GS, GE, HDS, HDE, HBS, HBE, MBS, MBE, NUMTOKENS
  };

const char vtkDataObjectGeneratorTokenStrings[NUMTOKENS][4] =
{ "ID1",
  "UF1",
  "RG1",
  "SG1",
  "PD1",
  "US1",
  "(",
  ")",
  "HD<",
  ">",
  "HB[",
  "]",
  "MB{",
  "}",
};

const char vtkDataObjectGeneratorTypeStrings[NUMTOKENS][30] =
{ "vtkImageData",
  "vtkUniformGrid",
  "vtkRectilinearGrid",
  "vtkStructuredGrid",
  "vtkPolyData",
  "vtkUnstructuredGrid",
  "NA",
  "NA",
  "vtkHierarchicalDataSet",
  "NA",
  "vtkHierarchicalBoxDataSet",
  "NA",
  "vtkMultiBlockDataSet",
  "NA",
};

//============================================================================
class vtkInternalStructureCache
{
  //a class to keep the overall structure in memory in. It is a simple tree 
  //where each node has a data set type flag and pointers to children
public:
  vtkInternalStructureCache()
  {
    type = -1;
    parent = NULL;
  }
  ~vtkInternalStructureCache()
  {
    vtkstd::vector<vtkInternalStructureCache *>::iterator it;
    for (it = this->children.begin();
         it != this->children.end();
         it++)
      {
      delete *it;
      }
  }
  vtkInternalStructureCache *add_dataset(int t)
  {
    vtkInternalStructureCache *child = new vtkInternalStructureCache();
    child->type = t;
    child->parent = this;
    children.push_back(child);    
    return child;
  }
  void print(int level=0)
  {
    for (int i = 0; i < level; i++)
      {
      cerr << " ";
      }
    if (type >= 0)
      {
      cerr << vtkDataObjectGeneratorTokenStrings[type] << endl;
      }
    else
      {
      cerr << "HOLDER" << endl;
      }
    vtkstd::vector<vtkInternalStructureCache *>::iterator it;
    for (it = this->children.begin();
         it != this->children.end();
         it++)
      {
      (*it)->print(level+1);
      }
    if (type = GS || type == HDS || type == HBS || type == MBS)
      {
      for (int i = 0; i < level; i++)
        {
        cerr << " ";
        }
      switch (type)
        {
        case GS:
          cerr << vtkDataObjectGeneratorTokenStrings[GE] << endl;
          break;
        case HDS:
          cerr << vtkDataObjectGeneratorTokenStrings[HDE] << endl;
          break;
        case HBS:
          cerr << vtkDataObjectGeneratorTokenStrings[HBE] << endl;
          break;
        case MBS:
          cerr << vtkDataObjectGeneratorTokenStrings[MBE] << endl;
          break;
        }        
      }
  }

  int type;
  vtkInternalStructureCache *parent;
  vtkstd::vector<vtkInternalStructureCache *> children;
};


//----------------------------------------------------------------------------
//search the head of the input string for one of the tokens we know how to
//do something with. If we see something, bump char ptr passed it, and return
//a code that says what we found. Skip over chars we don't recognize. When 
//nothing is left in the string return -1.
int vtkDataObjectGeneratorGetNextToken(char **str)
{
  if (!str || !*str)
    {
    return 0;
    }

  int len = strlen(*str);
  int l;
  while (len && *str)
    {
    for (int i = 0; i < NUMTOKENS; i++)
      {
      l = strlen(vtkDataObjectGeneratorTokenStrings[i]);
      if (len >= l
          && 
          !strncmp(*str, vtkDataObjectGeneratorTokenStrings[i], l))
        {
        *str+=l;
        return i;
        }
      }
    len--;
    *str+=1;
    }
  return -1;
}

//----------------------------------------------------------------------------
vtkInternalStructureCache *vtkDataObjectGeneratorParseStructure(char *Program)
{
  vtkInternalStructureCache *structure = new vtkInternalStructureCache();
  vtkInternalStructureCache *sptr = structure;

  char *ptr = Program;
  bool done = false;
  while (!done)
    {
    switch (vtkDataObjectGeneratorGetNextToken(&ptr))
      {
      case ID1:
        sptr->add_dataset(ID1);
        break;
      case UF1:
        sptr->add_dataset(UF1);
        break;
      case RG1:
        sptr->add_dataset(RG1);
        break;
      case SG1:
        sptr->add_dataset(SG1);
        break;
      case PD1:
        sptr->add_dataset(PD1);
        break;
      case US1:
        sptr->add_dataset(US1);
        break;
      case GS:
        {
        vtkInternalStructureCache *cptr = sptr->add_dataset(GS);
        sptr = cptr;
        break;
        }
      case HDS:
        {
        vtkInternalStructureCache *cptr = sptr->add_dataset(HDS);
        sptr = cptr;
        break;
        }
      case HBS:
        {
        vtkInternalStructureCache *cptr = sptr->add_dataset(HBS);
        sptr = cptr;
        break;
        }
      case MBS:
        {
        vtkInternalStructureCache *cptr = sptr->add_dataset(MBS);
        sptr = cptr;
        break;
        }
      case GE:
      case HDE:
      case HBE:
      case MBE:
        sptr = sptr->parent;
        break;
      default:
        done = true;
      }
    }

  return structure;
}

//---------------------------------------------------------------------------
vtkDataObject * vtkDataObjectGenerator::CreateOutputDataObjects(
  vtkInternalStructureCache *structure)
{
  vtkDataObject *outData;
  switch (structure->type)
    {
    case -1: //top holder it should hold a single data set, use it
    {
    if (!structure->children.size())
      {
      return NULL;
      }
    return this->CreateOutputDataObjects(structure->children.front());
    }
    case ID1:
    case UF1:
    case RG1:
    case SG1:
    case PD1:
    case US1:
    { 
    /*
    cerr 
      << "Creating " 
      << vtkDataObjectGeneratorTypeStrings[structure->type] 
      << endl;
    */
    outData = 
      vtkDataObjectTypes::NewDataObject(
        vtkDataObjectGeneratorTypeStrings[structure->type]);
    return outData;
    }
    case HDS:
    case HBS:
    case MBS:
    {
    //only create top level struct in RequestDataObject, do not recurse
    //the contents of the structure is cleared before RequestData anyway
    /*
    cerr 
      << "Creating " 
      << vtkDataObjectGeneratorTypeStrings[structure->type] 
      << endl;
    */
    outData = 
      vtkDataObjectTypes::NewDataObject(
        vtkDataObjectGeneratorTypeStrings[structure->type]);
    return outData;
    }
    case HDE: //should never be created
    case HBE: //should never be created
    case MBE: //should never be created
    case GS: //should be skipped over by MBS
    case GE: //should never be created
    default:
    //cerr << "UH OH" << endl;
    return NULL;
    }
}

//---------------------------------------------------------------------------
vtkDataObject * vtkDataObjectGenerator::FillOutputDataObjects(
  vtkInternalStructureCache *structure,
  int level,
  int stripe
  )
{
  vtkDataObject *outData;
  int t = structure->type;
  if (t != -1 &&
      t != HDE &&
      t != HBE &&
      t != MBE &&
      t != GS  &&
      t != GE)
    {    
    if (level==1 && ((stripe%this->Processors) != this->Rank)) 
      {
      //for parallel processing, each processor gets a different set of 
      //stripes of the data sets within the groups in the first level of 
      //composite data sets
      /*
      cerr << "Ignoring "
           << vtkDataObjectGeneratorTypeStrings[t] 
           << endl;
      */
      }
    else
      {
      /*
      cerr << "Filling " 
           << vtkDataObjectGeneratorTypeStrings[t] 
           << endl;
      */
      outData = vtkDataObjectTypes::NewDataObject(
        vtkDataObjectGeneratorTypeStrings[t]);
      }
    }

  switch (t)
    {
    case -1: //top level is a holder, use the single data set inside instead
    {
    if (!structure->children.size())
      {
      return NULL;
      }
    return this->FillOutputDataObjects(structure->children.front(), level+1);
    }
    case ID1:
    {
    this->MakeImageData1(vtkDataSet::SafeDownCast(outData));
    return outData;
    }
    case UF1:
    {
    this->MakeUniformGrid1(vtkDataSet::SafeDownCast(outData));
    return outData;
    }
    case RG1:
    {
    this->MakeRectilinearGrid1(vtkDataSet::SafeDownCast(outData));
    return outData;
    }
    case SG1:
    {
    this->MakeStructuredGrid1(vtkDataSet::SafeDownCast(outData));
    return outData;
    }
    case PD1:
    { 
    this->MakePolyData1(vtkDataSet::SafeDownCast(outData));
    return outData;
    }
    case US1:
    {
    this->MakeUnstructuredGrid1(vtkDataSet::SafeDownCast(outData));
    return outData;
    }
    case HDS:
    {
    //TODO: set metadata that describes parent child relationships
    vtkHierarchicalDataSet *hdo = 
      vtkHierarchicalDataSet::SafeDownCast(outData);

    this->YOffset += 1.0;
    //fill in the contents of this multiblockdataset
    //by iterate over the children of my children (which must be groups)
    hdo->SetNumberOfGroups(structure->children.size());
    vtkstd::vector<vtkInternalStructureCache *>::iterator git;
    vtkIdType gcnt = 0;
    for (git = structure->children.begin();
         git != structure->children.end();
         git++)
      {
      this->ZOffset += 1.0;
      vtkInternalStructureCache *gptr = *git; 
      //gpt->type should be a group
      vtkIdType nds = gptr->children.size();
      hdo->SetNumberOfDataSets(gcnt, nds);

      vtkstd::vector<vtkInternalStructureCache *>::iterator dit;
      vtkIdType dcnt = 0;
      for (dit = gptr->children.begin();
           dit != gptr->children.end();
           dit++)
        {
        this->XOffset += 1.0;
        vtkInternalStructureCache *dptr = *dit; 
        vtkDataObject *dobj = NULL;
        if (dptr->type == ID1 ||
            dptr->type == UF1 ||
            dptr->type == RG1 ||
            dptr->type == SG1 ||
            dptr->type == PD1 ||
            dptr->type == US1)
          {
          //restrict HierarchicalDataSets to be one level deep
          dobj = this->FillOutputDataObjects(dptr, level+1, dcnt);
          }
        hdo->SetDataSet(gcnt, dcnt, dobj);
        if (dobj)
          {
          dobj->Delete();
          }
        dcnt++;
        }
      this->XOffset -= dcnt;
      gcnt++;
      }
    this->ZOffset -= gcnt;
    this->YOffset -= 1.0;
    return outData;
    }
    case HBS:
    {
    //Making octrees, structured can grid up space arbitratily though

    vtkHierarchicalBoxDataSet *hbo = 
      vtkHierarchicalBoxDataSet::SafeDownCast(outData);

    hbo->SetNumberOfGroups(structure->children.size());
    vtkstd::vector<vtkInternalStructureCache *>::iterator git;
    vtkIdType gcnt = 0;
    for (git = structure->children.begin();
         git != structure->children.end();
         git++)
      {
      //cerr << "LVL=" << gcnt  << endl;

      vtkInternalStructureCache *gptr = *git; 
      //gptr->type should be a group

      vtkIdType nds = gptr->children.size();
      hbo->SetNumberOfDataSets(gcnt, nds); 

      //each of the dimensions of each parent cell are broken into this 
      //many pieces this must be the inverse of the spacing for the geometry 
      //to line up
      int refinement = 2; 
      hbo->SetRefinementRatio(gcnt,refinement);       

      vtkstd::vector<vtkInternalStructureCache *>::iterator dit;
      vtkIdType dcnt = 0; //TODO: read in a location to create sparse trees

      int maxchildren = (int)pow((double)8.0,gcnt); 
      //making octrees, this is total number of possible children in this level

      int r2 = (int)pow((double)refinement,gcnt); 
      //how many children across each dimension

      for (dit = gptr->children.begin();
           dit != gptr->children.end() 
             && dcnt<maxchildren //ignore extra children
             ;
           dit++)
        {
        //cerr << "DS=" << dcnt  << endl;
        vtkInternalStructureCache *dptr = *dit; 
        //dptr->type should be UF1

        //Figure out where in the parent level the cells of this new data 
        //set resides, this is used to create blanking parent child 
        //relationships
        //*2 is because each child is 2 cells across
        //+1 (in hi) is because we are counting cells inclusively. 
        //If children were 3x3x3 it would be *3+2
        //Note, other orderings are equally valid
        const int lo[3] = {dcnt/(r2*r2)%r2*2, dcnt/r2%r2*2, dcnt%r2*2};
        const int hi[3] = {dcnt/(r2*r2)%r2*2+1, dcnt/r2%r2*2+1, dcnt%r2*2+1};

        /*
        cerr << "LO=" << lo[0] << "," << lo[1] << "," << lo[2] << " "
             << "HI=" << hi[0] << "," << hi[1] << "," << hi[2] << endl;
        */
        vtkDataObject *dobj = NULL;
        //restrict HierarchicalBoxes's to contain only UniformGrids
        //until I make it read location to make sparse creation easy, put 
        //dummy dataobjects in as placeholders 
        if (dptr->type == UF1)
          {
          dobj = this->FillOutputDataObjects(dptr, level+1, dcnt);
          vtkUniformGrid *uf = vtkUniformGrid::SafeDownCast(dobj);
          //scale and translate the children to align with the parent the
          //blanking information
          double spacing = pow(0.5, gcnt+1);  //==1.0/(2*r2)
          uf->SetSpacing(spacing, spacing, spacing);
          double spa[3];
          uf->GetSpacing(spa);          
          //cerr << "SPACE=" <<spa[0] <<"," <<spa[1] <<"," <<spa[2] <<endl;
          double org[3];
          uf->SetOrigin(lo[0]*spacing, lo[1]*spacing, lo[2]*spacing);
          uf->GetOrigin(org);
          //cerr << "ORIGIN=" <<org[0] <<"," <<org[1] <<"," <<org[2] <<endl;
          uf->SetExtent(0,2,0,2,0,2); //octrees, 2 cells == 3 points across
          int ex[6];
          uf->GetExtent(ex);
          }        
        
        vtkAMRBox region(3,lo,hi);
        hbo->SetDataSet(gcnt, dcnt, 
                        region, vtkUniformGrid::SafeDownCast(dobj));
        if (dobj)
          {
          dobj->Delete();
          }
        dcnt++;
        }
      gcnt++;
      }
    hbo->GenerateVisibilityArrays();
    return outData;
    }
    case MBS:
    {
    vtkMultiBlockDataSet *mbo = 
      vtkMultiBlockDataSet::SafeDownCast(outData);

    this->YOffset += 1.0;
    //fill in the contents of this multiblockdataset
    //by iterating over the children of all my children (which must be groups)
    mbo->SetNumberOfGroups(structure->children.size());
    vtkstd::vector<vtkInternalStructureCache *>::iterator git;
    vtkIdType gcnt = 0;

    for (git = structure->children.begin();
         git != structure->children.end();
         git++)
      {
      this->ZOffset += 1.0;
      vtkInternalStructureCache *gptr = *git; 
      //gptr->type should be group
      vtkIdType nds = gptr->children.size();
      mbo->SetNumberOfDataSets(gcnt, nds);

      vtkstd::vector<vtkInternalStructureCache *>::iterator dit;
      vtkIdType dcnt = 0;
      for (dit = gptr->children.begin();
           dit != gptr->children.end();
           dit++)
        {
        this->XOffset += 1.0;

        vtkInternalStructureCache *dptr = *dit; 
        //dptr->type should not be group
        vtkDataObject *dobj = this->FillOutputDataObjects(dptr, level+1, dcnt);
        mbo->SetDataSet(gcnt, dcnt, dobj);
        if (dobj)
          {
          dobj->Delete();
          }
        dcnt++;
        }
      this->XOffset -= dcnt;

      gcnt++;
      }
    this->ZOffset -= gcnt;

    this->YOffset -= 1.0;

    return outData;
    }
    case HDE: //should never be created
    case HBE: //should never be created
    case MBE: //should never be created
    case GS: //should be skipped over by MBS
    case GE: //should never be created
    default:
    //cerr << "UH OH" << endl;
    return NULL;
    }
}

//----------------------------------------------------------------------------
vtkDataObjectGenerator::vtkDataObjectGenerator()
{
  this->SetNumberOfInputPorts(0);
  this->CellIdCounter = 0;
  this->PointIdCounter = 0;
  this->Program=NULL;
  this->Structure = NULL;
  this->XOffset = 0.0;
  this->YOffset = 0.0;
  this->ZOffset = 0.0;
}

//----------------------------------------------------------------------------
vtkDataObjectGenerator::~vtkDataObjectGenerator()
{
  this->SetProgram(NULL);
  if (this->Structure != NULL)
    {
    delete this->Structure;
    }
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Program: " << this->Program << "\n";
}

//----------------------------------------------------------------------------
int vtkDataObjectGenerator::RequestDataObject(vtkInformation *req,
                                          vtkInformationVector **inV,
                                          vtkInformationVector *outV)
{
  if (!this->Program)
     {
     vtkErrorMacro("No string to generate data objects for");
     return 1;
     }
  
  vtkInformation *outInfo = outV->GetInformationObject(0);
  vtkDataObject *outData = outInfo->Get(vtkDataObject::DATA_OBJECT());

  if (this->Structure != NULL)
    {
    delete this->Structure;
    }
  this->Structure = vtkDataObjectGeneratorParseStructure(this->Program);
  outData = this->CreateOutputDataObjects(this->Structure);
  if (outData)
    {
    outData->SetPipelineInformation(outInfo);
    }

  return VTK_OK;
}

//----------------------------------------------------------------------------
int vtkDataObjectGenerator::RequestInformation(vtkInformation *req,
                                               vtkInformationVector **inV,
                                               vtkInformationVector *outV)
{
  if (!this->Structure)
    {
    vtkErrorMacro("Program has not been parsed");
    return 1;
    }

  if (!this->Structure->children.size())
    {
    return 1;
    }

  //Say that this filter can break up its output into any number of pieces
  vtkInformation *outInfo = outV->GetInformationObject(0);
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  //If my output is an atomic structured type, fill in the whole extent info
  vtkInternalStructureCache *top = this->Structure->children.front();
  int t = top->type;
  if (t == ID1 ||
      t == RG1 ||
      t == SG1)
    {
    int ext[6];
    ext[0] = 0;
    ext[1] = 1;
    ext[2] = 0;
    ext[3] = 1;
    ext[4] = 0;
    ext[5] = 1;
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext ,6);
    double spacing[3];
    spacing[0] = 1.0;
    spacing[1] = 1.0;
    spacing[2] = 1.0;
    outInfo->Set(vtkDataObject::SPACING(),spacing,3);
    double origin[3];
    origin[0] = 0.0;
    origin[1] = 0.0;
    origin[2] = 0.0;
    outInfo->Set(vtkDataObject::ORIGIN(),origin,3);    
    }
  if (t == UF1)
    {
    int ext[6];
    ext[0] = 0;
    ext[1] = 2;
    ext[2] = 0;
    ext[3] = 2;
    ext[4] = 0;
    ext[5] = 2;
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext ,6);
    double spacing[3];
    spacing[0] = 0.5;
    spacing[1] = 0.5;
    spacing[2] = 0.5;
    outInfo->Set(vtkDataObject::SPACING(),spacing,3);
    double origin[3];
    origin[0] = 0.0;
    origin[1] = 0.0;
    origin[2] = 0.0;
    outInfo->Set(vtkDataObject::ORIGIN(),origin,3);    
    }

  //Could create vtkCompositeDataInformation here.
  return this->Superclass::RequestInformation(req, inV, outV);
}

//----------------------------------------------------------------------------
int vtkDataObjectGenerator::RequestUpdateExtent(vtkInformation *req,
                                                vtkInformationVector **inV,
                                                vtkInformationVector *outV)
{  
  //This is a source and doesn't have any inputs. 
  //I can defer this to the parent class because it does not have any 
  //inputs to request extent/pieces from dependent on what is requested by
  //my outputs.
  return this->Superclass::RequestUpdateExtent(req, inV, outV);
}

//----------------------------------------------------------------------------
int vtkDataObjectGenerator::RequestData(vtkInformation *req,
                                        vtkInformationVector **inV,
                                        vtkInformationVector *outV)
{
  if (!this->Structure)
    {
    vtkErrorMacro("Program has not been parsed");
    return 1;
    }

  //For parallel processing, this will stripe the datasets contained
  //in the first level of composite data sets.
  vtkInformation *outInfo = outV->GetInformationObject(0);
  this->Rank = 0;
  if (outInfo->Has(
        vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()))
    {
    this->Rank = 
      outInfo->Get(
        vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    }
  this->Processors = 1;
  if (outInfo->Has(
        vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()))
    {
    this->Processors =
      outInfo->Get(
        vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    }

  this->CellIdCounter = 0;
  this->PointIdCounter = 0;
  vtkDataObject *outData = this->FillOutputDataObjects(this->Structure, -1);
  outInfo->Set(vtkDataObject::DATA_OBJECT(), outData);     
  
  return VTK_OK;
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakeValues(vtkDataSet *ds)
{
  vtkIdTypeArray *ids;
  vtkDoubleArray *xcoords;
  vtkDoubleArray *ycoords;
  vtkDoubleArray *zcoords;
  vtkIdType num;

  num = ds->GetNumberOfCells();

  //give each cell a unique id and record its centroid
  ids=vtkIdTypeArray::New();
  ids->SetName("Cell Ids");
  ids->SetNumberOfComponents(1);
  ids->SetNumberOfTuples(num);
  xcoords = vtkDoubleArray::New();
  xcoords->SetName("Cell X");
  xcoords->SetNumberOfComponents(1);
  xcoords->SetNumberOfTuples(num);
  ycoords = vtkDoubleArray::New();
  ycoords->SetName("Cell Y");
  ycoords->SetNumberOfComponents(1);
  ycoords->SetNumberOfTuples(num);
  zcoords = vtkDoubleArray::New();
  zcoords->SetName("Cell Z");
  zcoords->SetNumberOfComponents(1);
  zcoords->SetNumberOfTuples(num);
  
  for (vtkIdType i = 0; i < num; i++)
    {
    ids->SetValue(i, this->CellIdCounter++);
    double *bds = ds->GetCell(i)->GetBounds();
    xcoords->SetValue(i, (bds[0]+bds[1])*0.5);
    ycoords->SetValue(i, (bds[2]+bds[3])*0.5);
    zcoords->SetValue(i, (bds[4]+bds[5])*0.5);
    }
  ds->GetCellData()->SetGlobalIds(ids);
  ds->GetCellData()->AddArray(xcoords);
  ds->GetCellData()->AddArray(ycoords);
  ds->GetCellData()->AddArray(zcoords);

  ids->Delete();  
  xcoords->Delete();
  ycoords->Delete();
  zcoords->Delete();

  //give each point a unique id and record its location
  num = ds->GetNumberOfPoints();
  ids=vtkIdTypeArray::New();
  ids->SetName("Point Ids");
  ids->SetNumberOfComponents(1);
  ids->SetNumberOfTuples(num);
  xcoords = vtkDoubleArray::New();
  xcoords->SetName("Point X");
  xcoords->SetNumberOfComponents(1);
  xcoords->SetNumberOfTuples(num);
  ycoords = vtkDoubleArray::New();
  ycoords->SetName("Point Y");
  ycoords->SetNumberOfComponents(1);
  ycoords->SetNumberOfTuples(num);
  zcoords = vtkDoubleArray::New();
  zcoords->SetName("Point Z");
  zcoords->SetNumberOfComponents(1);
  zcoords->SetNumberOfTuples(num);

  for (vtkIdType i = 0; i < num; i++)
    {
    ids->SetValue(i, this->PointIdCounter++);
    double *coords = ds->GetPoint(i);
    xcoords->SetValue(i, coords[0]);
    ycoords->SetValue(i, coords[1]);
    zcoords->SetValue(i, coords[2]);
    }
  ds->GetPointData()->SetGlobalIds(ids);
  ds->GetPointData()->AddArray(xcoords);
  ds->GetPointData()->AddArray(ycoords);
  ds->GetPointData()->AddArray(zcoords);

  ids->Delete();  
  xcoords->Delete();
  ycoords->Delete();
  zcoords->Delete();
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakeImageData1(vtkDataSet *ids)
{
  //ID1 == an ImageData of 1 voxel
  vtkImageData *ds = vtkImageData::SafeDownCast(ids);
  if (!ds)
    {
    return;
    }

  ds->Initialize();
  ds->SetDimensions(2,2,2); //1 cell
  ds->SetOrigin(this->XOffset,this->YOffset,this->ZOffset);
  ds->SetSpacing(1.0,1.0,1.0);

  this->MakeValues(ds);
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakeUniformGrid1(vtkDataSet *ids)
{
  //UF1 == an UniformGrid of 8 voxels
  vtkUniformGrid *ds = vtkUniformGrid::SafeDownCast(ids);
  if (!ds)
    {
    return;
    }

  ds->Initialize();
  ds->SetDimensions(3,3,3); //8 cells to make octrees
  ds->SetOrigin(this->XOffset,this->YOffset,this->ZOffset);
  ds->SetSpacing(0.5,0.5,0.5); 

  this->MakeValues(ds);
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakeRectilinearGrid1(vtkDataSet *ids)
{
  //RG1 = a RectilnearGrid of 1 voxel
  vtkRectilinearGrid *ds = vtkRectilinearGrid::SafeDownCast(ids);
  if (!ds)
    {
    return;
    }

  ds->Initialize();
  ds->SetDimensions(2,2,2); //1 cell
  vtkDoubleArray *da;

  da = vtkDoubleArray::New();
  da->SetNumberOfComponents(1);
  da->SetNumberOfTuples(2);
  da->SetName("X Coords");
  da->SetValue(0, this->XOffset);
  da->SetValue(1, this->XOffset+1.0);
  ds->SetXCoordinates(da);
  da->Delete();

  da = vtkDoubleArray::New();
  da->SetNumberOfComponents(1);
  da->SetNumberOfTuples(2);
  da->SetName("Y Coords");
  da->SetValue(0, this->YOffset);
  da->SetValue(1, this->YOffset+1.0);
  ds->SetYCoordinates(da);
  da->Delete();

  da = vtkDoubleArray::New();
  da->SetNumberOfComponents(1);
  da->SetNumberOfTuples(2);
  da->SetName("Z Coords");
  da->SetValue(0, this->ZOffset);
  da->SetValue(1, this->ZOffset+1.0);
  ds->SetZCoordinates(da);
  da->Delete();

  this->MakeValues(ds);
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakeStructuredGrid1(vtkDataSet *ids)
{
  //SG1 = a StructuredGrid of 1 voxel
  vtkStructuredGrid *ds = vtkStructuredGrid::SafeDownCast(ids);
  if (!ds)
    {
    return;
    }

  ds->Initialize();
  ds->SetDimensions(2,2,2); //1 cell
  vtkPoints *pts = vtkPoints::New();
  const double &XO = this->XOffset;
  const double &YO = this->YOffset;
  const double &ZO = this->ZOffset;
  pts->InsertNextPoint(XO+0.0, YO+0.0, ZO+0.0);
  pts->InsertNextPoint(XO+0.0, YO+0.0, ZO+1.0);
  pts->InsertNextPoint(XO+0.0, YO+1.0, ZO+0.0);
  pts->InsertNextPoint(XO+0.0, YO+1.0, ZO+1.0);
  pts->InsertNextPoint(XO+1.0, YO+0.0, ZO+0.0);
  pts->InsertNextPoint(XO+1.0, YO+0.0, ZO+1.0);
  pts->InsertNextPoint(XO+1.0, YO+1.0, ZO+0.0);
  pts->InsertNextPoint(XO+1.0, YO+1.0, ZO+1.0);
  ds->SetPoints(pts);
  pts->Delete();

  this->MakeValues(ds);
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakePolyData1(vtkDataSet *ids)
{
  //PD1 = a PolyData of 1 triangle
  vtkPolyData *ds = vtkPolyData::SafeDownCast(ids);
  if (!ds)
    {
    return;
    }

  ds->Initialize();
  vtkPoints *pts = vtkPoints::New();
  const double &XO = this->XOffset;
  const double &YO = this->YOffset;
  const double &ZO = this->ZOffset;
  pts->InsertNextPoint(XO+0.0, YO+0.0, ZO+0.0);
  pts->InsertNextPoint(XO+0.0, YO+0.0, ZO+1.0);
  pts->InsertNextPoint(XO+0.0, YO+1.0, ZO+0.0);
  ds->SetPoints(pts);
  pts->Delete();
  ds->Allocate();
  int verts[3] = {0,1,2};
  ds->InsertNextCell(VTK_TRIANGLE, 3, verts);

  this->MakeValues(ds);
}

//----------------------------------------------------------------------------
void vtkDataObjectGenerator::MakeUnstructuredGrid1(vtkDataSet *ids)
{
  //UG1 = an UnstructuredGrid of 1 triangle
  vtkUnstructuredGrid *ds = vtkUnstructuredGrid::SafeDownCast(ids);
  if (!ds)
    {
    return;
    }
  ds->Initialize();
  vtkPoints *pts = vtkPoints::New();
  const double &XO = this->XOffset;
  const double &YO = this->YOffset;
  const double &ZO = this->ZOffset;
  pts->InsertNextPoint(XO+0.0, YO+0.0, ZO+0.0);
  pts->InsertNextPoint(XO+0.0, YO+0.0, ZO+1.0);
  pts->InsertNextPoint(XO+0.0, YO+1.0, ZO+0.0);
  ds->SetPoints(pts);
  pts->Delete();
  ds->Allocate();
  int verts[3] = {0,1,2};
  ds->InsertNextCell(VTK_TRIANGLE, 3, verts);

  this->MakeValues(ds);
}
