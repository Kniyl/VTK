/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkHyperStreamline.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkHyperStreamline.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkHyperStreamline* vtkHyperStreamline::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkHyperStreamline");
  if(ret)
    {
    return (vtkHyperStreamline*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkHyperStreamline;
}




//
// Special classes for manipulating data
//BTX
class vtkHyperPoint { //;prevent man page generation
public:
    vtkHyperPoint(); // method sets up storage
    vtkHyperPoint &operator=(const vtkHyperPoint& hp); //for resizing
    
    float   X[3];    // position 
    int     CellId;  // cell
    int     SubId; // cell sub id
    float   P[3];    // parametric coords in cell 
    float   W[3];    // eigenvalues (sorted in decreasing value)
    float   *V[3];   // pointers to eigenvectors (also sorted)
    float   V0[3];   // storage for eigenvectors
    float   V1[3];
    float   V2[3];
    float   S;       // scalar value 
    float   D;       // distance travelled so far 
};
//ETX

class vtkHyperArray { //;prevent man page generation
public:
  vtkHyperArray();
  ~vtkHyperArray()
    {
      if (this->Array)
	{
	delete [] this->Array;
	}
    };
  int GetNumberOfPoints() {return this->MaxId + 1;};
  vtkHyperPoint *GetHyperPoint(int i) {return this->Array + i;};
  vtkHyperPoint *InsertNextHyperPoint() 
    {
    if ( ++this->MaxId >= this->Size )
      {
      this->Resize(this->MaxId);
      }
    return this->Array + this->MaxId;
    }
  vtkHyperPoint *Resize(int sz); //reallocates data
  void Reset() {this->MaxId = -1;};

  vtkHyperPoint *Array;  // pointer to data
  int MaxId;             // maximum index inserted thus far
  int Size;              // allocated size of data
  int Extend;            // grow array by this amount
  float Direction;       // integration direction
};

#define VTK_START_FROM_POSITION 0
#define VTK_START_FROM_LOCATION 1

vtkHyperPoint::vtkHyperPoint()
{
  this->V[0] = this->V0;
  this->V[1] = this->V1;
  this->V[2] = this->V2;
}

vtkHyperPoint& vtkHyperPoint::operator=(const vtkHyperPoint& hp)
{
  int i, j;

  for (i=0; i<3; i++) 
    {
    this->X[i] = hp.X[i];
    this->P[i] = hp.P[i];
    this->W[i] = hp.W[i];
    for (j=0; j<3; j++)
      {
      this->V[j][i] = hp.V[j][i];
      }
    }
  this->CellId = hp.CellId;
  this->SubId = hp.SubId;
  this->S = hp.S;
  this->D = hp.D;

  return *this;
}

vtkHyperArray::vtkHyperArray()
{
  this->MaxId = -1; 
  this->Array = new vtkHyperPoint[1000];
  this->Size = 1000;
  this->Extend = 5000;
  this->Direction = VTK_INTEGRATE_FORWARD;
}

vtkHyperPoint *vtkHyperArray::Resize(int sz)
{
  vtkHyperPoint *newArray;
  int newSize, i;

  if (sz >= this->Size)
    {
    newSize = this->Size + 
      this->Extend*(((sz-this->Size)/this->Extend)+1);
    
}
  else
    {
    newSize = sz;
    }

  newArray = new vtkHyperPoint[newSize];

  for (i=0; i<sz; i++)
    {
    newArray[i] = this->Array[i];
    }

  this->Size = newSize;
  delete [] this->Array;
  this->Array = newArray;

  return this->Array;
}

// Construct object with initial starting position (0,0,0); integration step 
// length 0.2; step length 0.01; forward integration; terminal eigenvalue 0.0;
// number of sides 6; radius 0.5; and logarithmic scaling off.
vtkHyperStreamline::vtkHyperStreamline()
{
  this->StartFrom = VTK_START_FROM_POSITION;
  this->StartPosition[0] = this->StartPosition[1] = this->StartPosition[2] = 0.0;

  this->StartCell = 0;
  this->StartSubId = 0;
  this->StartPCoords[0] = this->StartPCoords[1] = this->StartPCoords[2] = 0.5;

  this->Streamers = NULL;

  this->MaximumPropagationDistance = 100.0;
  this->IntegrationStepLength = 0.2;
  this->StepLength = 0.01;
  this->IntegrationDirection = VTK_INTEGRATE_FORWARD;
  this->TerminalEigenvalue = 0.0;
  this->NumberOfSides = 6;
  this->Radius = 0.5;
  this->LogScaling = 0;
  this->IntegrationEigenvector = 0; //Major eigenvector
}

vtkHyperStreamline::~vtkHyperStreamline()
{
  if ( this->Streamers )
    {
    delete [] this->Streamers;
    }
}

// Specify the start of the hyperstreamline in the cell coordinate system. 
// That is, cellId and subId (if composite cell), and parametric coordinates.
void vtkHyperStreamline::SetStartLocation(int cellId, int subId, float pcoords[3])
{
  if ( cellId != this->StartCell || subId != this->StartSubId ||
  pcoords[0] !=  this->StartPCoords[0] || 
  pcoords[1] !=  this->StartPCoords[1] || 
  pcoords[2] !=  this->StartPCoords[2] )
    {
    this->Modified();
    this->StartFrom = VTK_START_FROM_LOCATION;

    this->StartCell = cellId;
    this->StartSubId = subId;
    this->StartPCoords[0] = pcoords[0];
    this->StartPCoords[1] = pcoords[1];
    this->StartPCoords[2] = pcoords[2];
    }
}

// Specify the start of the hyperstreamline in the cell coordinate system. 
// That is, cellId and subId (if composite cell), and parametric coordinates.
void vtkHyperStreamline::SetStartLocation(int cellId, int subId, float r, float s, float t)
{
  float pcoords[3];
  pcoords[0] = r;
  pcoords[1] = s;
  pcoords[2] = t;

  this->SetStartLocation(cellId, subId, pcoords);
}

// Get the starting location of the hyperstreamline in the cell coordinate
// system. Returns the cell that the starting point is in.
int vtkHyperStreamline::GetStartLocation(int& subId, float pcoords[3])
{
  subId = this->StartSubId;
  pcoords[0] = this->StartPCoords[0];
  pcoords[1] = this->StartPCoords[1];
  pcoords[2] = this->StartPCoords[2];
  return this->StartCell;
}

// Specify the start of the hyperstreamline in the global coordinate system. 
// Starting from position implies that a search must be performed to find 
// initial cell to start integration from.
void vtkHyperStreamline::SetStartPosition(float x[3])
{
  if ( x[0] != this->StartPosition[0] || x[1] != this->StartPosition[1] || 
  x[2] != this->StartPosition[2] )
    {
    this->Modified();
    this->StartFrom = VTK_START_FROM_POSITION;

    this->StartPosition[0] = x[0];
    this->StartPosition[1] = x[1];
    this->StartPosition[2] = x[2];
    }
}

// Specify the start of the hyperstreamline in the global coordinate system. 
// Starting from position implies that a search must be performed to find 
// initial cell to start integration from.
void vtkHyperStreamline::SetStartPosition(float x, float y, float z)
{
  float pos[3];
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;

  this->SetStartPosition(pos);
}

// Get the start position of the hyperstreamline in global x-y-z coordinates.
float *vtkHyperStreamline::GetStartPosition()
{
  return this->StartPosition;
}

// Use the major eigenvector field as the vector field through which to 
// integrate. The major eigenvector is the eigenvector whose corresponding
// eigenvalue is closest to positive infinity.
void vtkHyperStreamline::IntegrateMajorEigenvector()
{
  if ( this->IntegrationEigenvector != 0 )
    {
    this->Modified();
    this->IntegrationEigenvector = 0;
  }
}

// Use the major eigenvector field as the vector field through which to 
// integrate. The major eigenvector is the eigenvector whose corresponding
// eigenvalue is between the major and minor eigenvalues.
void vtkHyperStreamline::IntegrateMediumEigenvector()
{
  if ( this->IntegrationEigenvector != 1 )
    {
    this->Modified();
    this->IntegrationEigenvector = 1;
  }
}

// Use the major eigenvector field as the vector field through which to 
// integrate. The major eigenvector is the eigenvector whose corresponding
// eigenvalue is closest to negative infinity.
void vtkHyperStreamline::IntegrateMinorEigenvector()
{
  if ( this->IntegrationEigenvector != 2 )
    {
    this->Modified();
    this->IntegrationEigenvector = 2;
  }
}

// Make sure coordinate systems are consistent
static void FixVectors(float **prev, float **current, int iv, int ix, int iy)
{
  float p0[3], p1[3], p2[3];
  float v0[3], v1[3], v2[3];
  float temp[3];
  int i;

  for (i=0; i<3; i++)
    {
    v0[i] = current[i][iv];
    v1[i] = current[i][ix];
    v2[i] = current[i][iy];
    }

  if ( prev == NULL ) //make sure coord system is right handed
    {
    vtkMath::Cross(v0,v1,temp);
    if ( vtkMath::Dot(v2,temp) < 0.0 )
      {
      for (i=0; i<3; i++)
	{
	current[i][iy] *= -1.0;
	}
      }
    }

  else //make sure vectors consistent from one point to the next
    {
    for (i=0; i<3; i++)
      {
      p0[i] = prev[i][iv];
      p1[i] = prev[i][ix];
      p2[i] = prev[i][iy];
      }
    if ( vtkMath::Dot(p0,v0) < 0.0 )
      {
      for (i=0; i<3; i++)
	{
	current[i][iv] *= -1.0;
	}
      }
    if ( vtkMath::Dot(p1,v1) < 0.0 )
      {
      for (i=0; i<3; i++)
	{
	current[i][ix] *= -1.0;
	}
      }
    if ( vtkMath::Dot(p2,v2) < 0.0 )
      {
      for (i=0; i<3; i++)
	{
	current[i][iy] *= -1.0;
	}
      }
    }
}

void vtkHyperStreamline::Execute()
{
  vtkDataSet *input = this->GetInput();
  vtkPointData *pd=input->GetPointData();
  vtkScalars *inScalars;
  vtkTensors *inTensors;
  vtkTensor *tensor;
  vtkHyperPoint *sNext, *sPtr;
  int i, j, k, ptId, subId, iv, ix, iy;
  vtkCell *cell;
  float ev[3], xNext[3];
  float d, step, dir, tol2, p[3];
  float *w=new float[input->GetMaxCellSize()], dist2;
  float closestPoint[3];
  float *m[3], *v[3];
  float m0[3], m1[3], m2[3];
  float v0[3], v1[3], v2[3];
  vtkTensors *cellTensors;
  vtkScalars *cellScalars;
  // set up working matrices
  v[0] = v0; v[1] = v1; v[2] = v2; 
  m[0] = m0; m[1] = m1; m[2] = m2; 

  vtkDebugMacro(<<"Generating hyperstreamline(s)");
  this->NumberOfStreamers = 0;

  if ( ! (inTensors=pd->GetTensors()) )
    {
    vtkErrorMacro(<<"No tensor data defined!");
    return;
    }

  cellTensors = vtkTensors::New();
  cellScalars = vtkScalars::New();
  cellTensors->Allocate(VTK_CELL_SIZE);
  cellScalars->Allocate(VTK_CELL_SIZE);
  
  
  inScalars = pd->GetScalars();
  tol2 = input->GetLength() / 1000.0;
  tol2 = tol2 * tol2;
  iv = this->IntegrationEigenvector;
  ix = (iv + 1) % 3;
  iy = (iv + 2) % 3;
  //
  // Create starting points
  //
  this->NumberOfStreamers = 1;
 
  if ( this->IntegrationDirection == VTK_INTEGRATE_BOTH_DIRECTIONS )
    {
    this->NumberOfStreamers *= 2;
    }

  this->Streamers = new vtkHyperArray[this->NumberOfStreamers];

  if ( this->StartFrom == VTK_START_FROM_POSITION )
    {
    sPtr = this->Streamers[0].InsertNextHyperPoint();
    for (i=0; i<3; i++)
      {
      sPtr->X[i] = this->StartPosition[i];
      }
    sPtr->CellId = input->FindCell(this->StartPosition, NULL, (-1), 0.0, 
                                   sPtr->SubId, sPtr->P, w);
    }

  else //VTK_START_FROM_LOCATION
    {
    sPtr = this->Streamers[0].InsertNextHyperPoint();
    cell =  input->GetCell(sPtr->CellId);
    cell->EvaluateLocation(sPtr->SubId, sPtr->P, sPtr->X, w);
    }
  //
  // Finish initializing each hyperstreamline
  //
  this->Streamers[0].Direction = 1.0;
  sPtr = this->Streamers[0].GetHyperPoint(0);
  sPtr->D = 0.0;
  if ( sPtr->CellId >= 0 ) //starting point in dataset
    {
    cell = input->GetCell(sPtr->CellId);
    cell->EvaluateLocation(sPtr->SubId, sPtr->P, xNext, w);

    inTensors->GetTensors(cell->PointIds, cellTensors);

    // interpolate tensor, compute eigenfunctions
    for (j=0; j<3; j++)
      {
      for (i=0; i<3; i++)
	{
	m[i][j] = 0.0;
	}
      }
    for (k=0; k < cell->GetNumberOfPoints(); k++)
      {
      tensor = cellTensors->GetTensor(k);
      for (j=0; j<3; j++) 
        {
        for (i=0; i<3; i++) 
          {
          m[i][j] += tensor->GetComponent(i,j) * w[k];
          }
        }
      }

    vtkMath::Jacobi(m, sPtr->W, sPtr->V);
    FixVectors(NULL, sPtr->V, iv, ix, iy);

    if ( inScalars ) 
      {
      inScalars->GetScalars(cell->PointIds, cellScalars);
      for (sPtr->S=0, i=0; i < cell->GetNumberOfPoints(); i++)
	{
        sPtr->S += cellScalars->GetScalar(i) * w[i];
	}
      }

    if ( this->IntegrationDirection == VTK_INTEGRATE_BOTH_DIRECTIONS )
      {
      this->Streamers[1].Direction = -1.0;
      sNext = this->Streamers[1].InsertNextHyperPoint();
      *sNext = *sPtr;
      }
    else if ( this->IntegrationDirection == VTK_INTEGRATE_BACKWARD )
      {
      this->Streamers[0].Direction = -1.0;
      }
    } //for hyperstreamline in dataset
  //
  // For each hyperstreamline, integrate in appropriate direction (using RK2).
  //
  for (ptId=0; ptId < this->NumberOfStreamers; ptId++)
    {
    //get starting step
    sPtr = this->Streamers[ptId].GetHyperPoint(0);
    if ( sPtr->CellId < 0 )
      {
      continue;
      }

    dir = this->Streamers[ptId].Direction;
    cell = input->GetCell(sPtr->CellId);
    cell->EvaluateLocation(sPtr->SubId, sPtr->P, xNext, w);
    step = this->IntegrationStepLength * sqrt((double)cell->GetLength2());
    inTensors->GetTensors(cell->PointIds, cellTensors);
    if ( inScalars ) {inScalars->GetScalars(cell->PointIds, cellScalars);}

    //integrate until distance has been exceeded
    while ( sPtr->CellId >= 0 && fabs(sPtr->W[0]) > this->TerminalEigenvalue &&
	    sPtr->D < this->MaximumPropagationDistance )
      {

      //compute updated position using this step (Euler integration)
      for (i=0; i<3; i++)
	{
        xNext[i] = sPtr->X[i] + dir * step * sPtr->V[i][iv];
	}

      //compute updated position using updated step
      cell->EvaluatePosition(xNext, closestPoint, subId, p, dist2, w);

      //interpolate tensor
      for (j=0; j<3; j++)
	{
	for (i=0; i<3; i++)
	  {
	  m[i][j] = 0.0;
	  }
	}
      for (k=0; k < cell->GetNumberOfPoints(); k++)
        {
        tensor = cellTensors->GetTensor(k);
        for (j=0; j<3; j++) 
          {
          for (i=0; i<3; i++) 
            {
            m[i][j] += tensor->GetComponent(i,j) * w[k];
            }
          }
        }

      vtkMath::Jacobi(m, ev, v);
      FixVectors(sPtr->V, v, iv, ix, iy);

      //now compute final position
      for (i=0; i<3; i++)
	{
        xNext[i] = sPtr->X[i] + 
                   dir * (step/2.0) * (sPtr->V[i][iv] + v[i][iv]);
	}
      sNext = this->Streamers[ptId].InsertNextHyperPoint();

      if ( cell->EvaluatePosition(xNext, closestPoint, sNext->SubId, 
      sNext->P, dist2, w) )
        { //integration still in cell
        for (i=0; i<3; i++)
	  {
	  sNext->X[i] = closestPoint[i];
	  }
        sNext->CellId = sPtr->CellId;
        sNext->SubId = sPtr->SubId;
        }
      else
        { //integration has passed out of cell
        sNext->CellId = input->FindCell(xNext, cell, sPtr->CellId, tol2, 
                                        sNext->SubId, sNext->P, w);
        if ( sNext->CellId >= 0 ) //make sure not out of dataset
          {
          for (i=0; i<3; i++)
	    {
	    sNext->X[i] = xNext[i];
	    }
          cell = input->GetCell(sNext->CellId);
          inTensors->GetTensors(cell->PointIds, cellTensors);
          if (inScalars){inScalars->GetScalars(cell->PointIds, cellScalars);}
          step = this->IntegrationStepLength * sqrt((double)cell->GetLength2());
          }
        }

      if ( sNext->CellId >= 0 )
        {
        cell->EvaluateLocation(sNext->SubId, sNext->P, xNext, w);
        for (j=0; j<3; j++)
	  {
	  for (i=0; i<3; i++)
	    {
	    m[i][j] = 0.0;
	    }
	  }
        for (k=0; k < cell->GetNumberOfPoints(); k++)
          {
          tensor = cellTensors->GetTensor(k);
          for (j=0; j<3; j++) 
            {
            for (i=0; i<3; i++) 
              {
              m[i][j] += tensor->GetComponent(i,j) * w[k];
              }
            }
          }

        vtkMath::Jacobi(m, sNext->W, sNext->V);
        FixVectors(sPtr->V, sNext->V, iv, ix, iy);

        if ( inScalars )
	  {
          for (sNext->S=0.0, i=0; i < cell->GetNumberOfPoints(); i++)
	    {
            sNext->S += cellScalars->GetScalar(i) * w[i];
	    }
	  }
        d = sqrt((double)vtkMath::Distance2BetweenPoints(sPtr->X,sNext->X));
        sNext->D = sPtr->D + d;
        }

      sPtr = sNext;

      }//for elapsed time

    } //for each hyperstreamline

  this->BuildTube();

  delete [] w;
  cellTensors->Delete();
  cellScalars->Delete();  
}

void vtkHyperStreamline::BuildTube()
{
  vtkHyperPoint *sPrev, *sPtr;
  vtkPoints *newPts;
  vtkVectors *newVectors;
  vtkNormals *newNormals;
  vtkScalars *newScalars=NULL;
  vtkCellArray *newStrips;
  int i, ptId, j, id, k, i1, i2;
  int npts, ptOffset=0;
  float dOffset, x[3], v[3], s, r, r1[3], r2[3], stepLength;
  float xT[3], sFactor, normal[3], w[3];
  float theta=2.0*vtkMath::Pi()/this->NumberOfSides;
  vtkPointData *outPD;
  vtkDataSet *input = this->GetInput();
  vtkPolyData *output = this->GetOutput();
  int iv, ix, iy, numIntPts;
  //
  // Initialize
  //
  vtkDebugMacro(<<"Creating hyperstreamline tube");
  if ( this->NumberOfStreamers <= 0 )
    {
    return;
    }

  stepLength = input->GetLength() * this->StepLength;
  outPD = output->GetPointData();

  iv = this->IntegrationEigenvector;
  ix = (iv+1) % 3;
  iy = (iv+2) % 3;
  //
  // Allocate
  //
  newPts  = vtkPoints::New();
  newPts ->Allocate(2500);
  if ( input->GetPointData()->GetScalars() )
    {
    newScalars = vtkScalars::New();
    newScalars->Allocate(2500);
    }
  newVectors = vtkVectors::New();
  newVectors->Allocate(2500);
  newNormals = vtkNormals::New();
  newNormals->Allocate(2500);
  newStrips = vtkCellArray::New();
  newStrips->Allocate(newStrips->EstimateSize(3*this->NumberOfStreamers,
                                              VTK_CELL_SIZE));
  //
  // Loop over all hyperstreamlines generating points
  //
  for (ptId=0; ptId < this->NumberOfStreamers; ptId++)
    {
    if ( (numIntPts=this->Streamers[ptId].GetNumberOfPoints()) < 2 )
      {
      continue;
      }
    sPrev = this->Streamers[ptId].GetHyperPoint(0);
    sPtr = this->Streamers[ptId].GetHyperPoint(1);

    // compute scale factor
    i = (sPrev->W[ix] > sPrev->W[iy] ? ix : iy);
    if ( sPrev->W[i] == 0.0 )
      {
      sFactor = 1.0;
      }
    else
      {
      sFactor = this->Radius / sPrev->W[i];
      }

    if ( numIntPts == 2 && sPtr->CellId < 0 )
      {
      continue;
      }

    dOffset = sPrev->D;

    for ( npts=0, i=1; i < numIntPts && sPtr->CellId >= 0;
    i++, sPrev=sPtr, sPtr=this->Streamers[ptId].GetHyperPoint(i) )
      {
  //
  // Bracket steps and construct tube points
  //
      while ( dOffset >= sPrev->D && dOffset < sPtr->D )
        {
        r = (dOffset - sPrev->D) / (sPtr->D - sPrev->D);

        for (j=0; j<3; j++) //compute point in center of tube
          {
          x[j] = sPrev->X[j] + r * (sPtr->X[j] - sPrev->X[j]);
          v[j] = sPrev->V[j][iv] + r * (sPtr->V[j][iv] - sPrev->V[j][iv]);
          r1[j] = sPrev->V[j][ix] + r * (sPtr->V[j][ix] - sPrev->V[j][ix]);
          r2[j] = sPrev->V[j][iy] + r * (sPtr->V[j][iy] - sPrev->V[j][iy]);
          w[j] = sPrev->W[j] + r * (sPtr->W[j] - sPrev->W[j]);
          }

        // construct points around tube
        for (k=0; k < this->NumberOfSides; k++)
          {
          for (j=0; j<3; j++) 
            {
            normal[j] = w[ix]*r1[j]*cos((double)k*theta) + 
                        w[iy]*r2[j]*sin((double)k*theta);
            xT[j] = x[j] + sFactor * normal[j];
            }
          id = newPts->InsertNextPoint(xT);
          newVectors->InsertVector(id,v);
          vtkMath::Normalize(normal);
          newNormals->InsertNormal(id,normal);
          }

        if ( newScalars ) //add scalars around tube
          {
          s = sPrev->S + r * (sPtr->S - sPrev->S);
          for (k=0; k<this->NumberOfSides; k++)
	    {
            newScalars->InsertNextScalar(s);
	    }
          }

        npts++;
        dOffset += stepLength;

        } //while
      } //for this hyperstreamline

    //
    // Generate the strips for this hyperstreamline
    //
    for (k=0; k<this->NumberOfSides; k++)
      {
      i1 = (k+1) % this->NumberOfSides;
      newStrips->InsertNextCell(npts*2);
      for (i=0; i < npts; i++) 
        {
        //make sure strip definition consistent with normals
        if (this->Streamers[ptId].Direction > 0.0)
	  {
	  i2 = i*this->NumberOfSides;
	  }
        else
	  {
	  i2 = (npts - i - 1) * this->NumberOfSides;
	  }
        newStrips->InsertCellPoint(ptOffset+i2+k);
        newStrips->InsertCellPoint(ptOffset+i2+i1);
        }
      }//for all tube sides

    ptOffset += this->NumberOfSides*npts;

    } //for all hyperstreamlines
  //
  // Update ourselves
  //
  output->SetPoints(newPts);
  newPts->Delete();

  output->SetStrips(newStrips);
  newStrips->Delete();

  if ( newScalars )
    {
    outPD->SetScalars(newScalars);
    newScalars->Delete();
    }

  outPD->SetNormals(newNormals);
  newNormals->Delete();

  outPD->SetVectors(newVectors);
  newVectors->Delete();

  output->Squeeze();
}

void vtkHyperStreamline::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkDataSetToPolyDataFilter::PrintSelf(os,indent);

  if ( this->StartFrom == VTK_START_FROM_POSITION )
    {
    os << indent << "Starting Position: (" << this->StartPosition[0] << ","
       << this->StartPosition[1] << ", " << this->StartPosition[2] << ")\n";
    }
  else 
    {
    os << indent << "Starting Location:\n\tCell: " << this->StartCell 
       << "\n\tSubId: " << this->StartSubId << "\n\tP.Coordinates: ("
       << this->StartPCoords[0] << ", " 
       << this->StartPCoords[1] << ", " 
       << this->StartPCoords[2] << ")\n";
    }

  os << indent << "Maximum Propagation Distance: " 
     << this->MaximumPropagationDistance << "\n";

  if ( this->IntegrationDirection == VTK_INTEGRATE_FORWARD )
    {
    os << indent << "Integration Direction: FORWARD\n";
    }
  else if ( this->IntegrationDirection == VTK_INTEGRATE_BACKWARD )
    {
    os << indent << "Integration Direction: BACKWARD\n";
    }
  else
    {
    os << indent << "Integration Direction: FORWARD & BACKWARD\n";
    }

  os << indent << "Integration Step Length: " << this->IntegrationStepLength << "\n";
  os << indent << "Step Length: " << this->StepLength << "\n";

  os << indent << "Terminal Eigenvalue: " << this->TerminalEigenvalue << "\n";

  os << indent << "Radius: " << this->Radius << "\n";
  os << indent << "Number Of Sides: " << this->NumberOfSides << "\n";
  os << indent << "Logarithmic Scaling: " << (this->LogScaling ? "On\n" : "Off\n");
  
  if ( this->IntegrationEigenvector == 0 )
    {
    os << indent << "Integrate Along Major Eigenvector\n";
    }
  else if ( this->IntegrationEigenvector == 1 )
    {
    os << indent << "Integrate Along Medium Eigenvector\n";
    }
  else
    {
    os << indent << "Integrate Along Minor Eigenvector\n";
    }
}


