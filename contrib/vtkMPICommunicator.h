/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMPICommunicator.h
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
// .NAME vtkMPICommunicator - Class for creating user defined MPI communicators.

// .SECTION Description
// This class can be used to create user defined MPI communicators.
// The actual creation (with MPI_Comm_create) occurs in Initialize
// which takes as arguments a super-communicator and a group of
// process ids. The new communicator is created by including the
// processes contained in the group. The global communicator 
// (equivalent to MPI_COMM_WORLD) can be obtained using the class
// method GetWorldCommunicator. It is important to note that 
// this communicator should not be used on the processes not contained
// in the group. For example, if the group contains processes 0 and 1,
// controller->SetCommunicator(communicator) would cause an MPI error
// on any other process.

// .SECTION See Also
// vtkMPIController vtkMPIGroup

#ifndef __vtkMPICommunicator_h
#define __vtkMPICommunicator_h

#include "vtkMPIGroup.h"
#include "mpi.h"

class vtkMPIController;

class VTK_EXPORT vtkMPICommunicator : public vtkObject
{

public:

  vtkTypeMacro( vtkMPICommunicator,vtkObject);
  
  // Description:
  // Creates an empty communicator.
  static vtkMPICommunicator* New();

  // Description:
  // Returns the singleton which behaves as the global
  // communicator (MPI_COMM_WORLD)
  static vtkMPICommunicator* GetWorldCommunicator();
  
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Used to initialize (i.e. create the underlying MPI_Comm)
  // the communicator. Note that group is also stored in
  // an instance variable (the reference count of group
  // is increased by 1). group->Delete() can be safely invoked after 
  // this.
  int Initialize(vtkMPICommunicator* mpiComm, vtkMPIGroup* group);

//BTX

  friend vtkMPIController;

//ETX

protected:

  vtkSetObjectMacro(Group, vtkMPIGroup);

  static vtkMPICommunicator* WorldCommunicator;

  MPI_Comm* Handle;
  vtkMPIGroup* Group;

  int Initialized;

  vtkMPICommunicator();
  ~vtkMPICommunicator();
  vtkMPICommunicator(const vtkMPICommunicator&) {};
  void operator=(const vtkMPICommunicator&) {};

};

#endif




