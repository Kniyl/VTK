/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestPBGLGraphSQLReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/* 
 * Copyright (C) 2008 The Trustees of Indiana University.
 * Use, modification and distribution is subject to the Boost Software
 * License, Version 1.0. (See http://www.boost.org/LICENSE_1_0.txt)
 */

#include <mpi.h>

#include "vtkDataObject.h"
#include "vtkDataSetAttributes.h"
#include "vtkEdgeListIterator.h"
#include "vtkGraph.h"
#include "vtkInformation.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkPBGLDistributedGraphHelper.h"
#include "vtkPBGLGraphSQLReader.h"
#include "vtkSmartPointer.h"
#include "vtkSQLiteDatabase.h"
#include "vtkSQLQuery.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkVertexListIterator.h"
#include "vtkVariant.h"
#include "vtkVariantArray.h"

#include <vtksys/ios/sstream>
#include <boost/graph/distributed/mpi_process_group.hpp>

#define myassert(Cond)                                  \
  if (!(Cond))                                          \
    {                                                   \
      cerr << "error (" __FILE__ ":" << dec << __LINE__ \
           << ") assertion \"" #Cond "\" failed."       \
           << endl;                                     \
    MPI_Abort(MPI_COMM_WORLD, -1);                      \
    }

void TestPSQLGraphReader()
{
  vtksys_ios::ostringstream oss;
  // Make a database containing a cycle.
  int vertices = 11;
  vtkSmartPointer<vtkSQLiteDatabase> db =
    vtkSmartPointer<vtkSQLiteDatabase>::New();
  db->SetDatabaseFileName(":memory:");
  bool ok = db->Open();
  if (!ok)
    {
    cerr << "Could not open database!" << endl;
    cerr << db->GetLastErrorText() << endl;
    return;
    }
  vtkSmartPointer<vtkSQLQuery> query;
  query.TakeReference(db->GetQueryInstance());
  query->SetQuery("create table vertices (id INTEGER)");
  query->Execute();
  for (int i = 0; i < vertices; ++i)
    {
    oss.str("");
    oss << "insert into vertices values(" << i << ")" << endl;
    query->SetQuery(oss.str().c_str());
    query->Execute();
    }
  query->SetQuery("create table edges (source INTEGER, target INTEGER)");
  query->Execute();
  for (int i = 0; i < vertices; ++i)
    {
    oss.str("");
    oss << "insert into edges values(" << i << ", "
      << (i+1)%vertices << ")" << endl;
    query->SetQuery(oss.str().c_str());
    query->Execute();
    }

  vtkSmartPointer<vtkPBGLGraphSQLReader> reader =
    vtkSmartPointer<vtkPBGLGraphSQLReader>::New();
  reader->SetDatabase(db);
  reader->SetVertexTable("vertices");
  reader->SetEdgeTable("edges");
  reader->SetVertexIdField("id");
  reader->SetSourceField("source");
  reader->SetTargetField("target");
  vtkStreamingDemandDrivenPipeline* exec =
    vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  vtkSmartPointer<vtkPBGLDistributedGraphHelper> helper =
    vtkSmartPointer<vtkPBGLDistributedGraphHelper>::New();
  int total = num_processes(helper->GetProcessGroup());
  int rank = process_id(helper->GetProcessGroup());
  reader->UpdateInformation();
  exec->SetUpdateNumberOfPieces(exec->GetOutputInformation(0), total);
  exec->SetUpdatePiece(exec->GetOutputInformation(0), rank);
  reader->Update();
  vtkGraph* output = reader->GetOutput();
  vtkSmartPointer<vtkEdgeListIterator> it =
    vtkSmartPointer<vtkEdgeListIterator>::New();
  output->GetEdges(it);
  while (it->HasNext())
    {
    vtkEdgeType e = it->Next();
    cerr << "PROCESS " << rank << ": " << hex << e.Id << " (" << e.Source << "," << e.Target << ")" << endl;
    }
  vtkSmartPointer<vtkVertexListIterator> vit =
    vtkSmartPointer<vtkVertexListIterator>::New();
  output->GetVertices(vit);
  while (vit->HasNext())
    {
    vtkIdType v = vit->Next();
    cerr << "PROCESS " << rank << ": " << hex << v << endl;
    }
}

//----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  TestPSQLGraphReader();
  cerr << "finalizing." << endl;
  MPI_Finalize();
  cerr << "done." << endl;
  return 0;
}
