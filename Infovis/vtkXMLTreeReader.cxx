/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXMLTreeReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------*/

#include "vtkXMLTreeReader.h"

#include "vtkBitArray.h"
#include "vtkInformation.h"
#include "vtkIdTypeArray.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkStringArray.h"
#include "vtkTree.h"

#include "vtk_libxml2.h"
#include VTKLIBXML2_HEADER(parser.h)
#include VTKLIBXML2_HEADER(tree.h)

vtkCxxRevisionMacro(vtkXMLTreeReader, "1.8");
vtkStandardNewMacro(vtkXMLTreeReader);

const char * vtkXMLTreeReader::TagNameField = ".tagname";
const char * vtkXMLTreeReader::CharDataField = ".chardata";

vtkXMLTreeReader::vtkXMLTreeReader()
{
  this->FileName = 0;
  this->XMLString = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->ReadCharData = 0;
  this->MaskArrays = 0;
}

vtkXMLTreeReader::~vtkXMLTreeReader()
{
  this->SetFileName(0);
  this->SetXMLString(0);
}

void vtkXMLTreeReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " 
     << (this->FileName ? this->FileName : "(none)") << endl;
  os << indent << "ReadCharData: "
     << (this->ReadCharData ? "on" : "off") << endl;
  os << indent << "MaskArrays: "
     << (this->MaskArrays ? "on" : "off") << endl;
  os << indent << "XMLString: " 
     << (this->XMLString ? this->XMLString : "(none)") << endl;
}

void vtkXMLTreeReaderProcessElement(vtkMutableDirectedGraph *tree,
   vtkIdType parent, xmlNode *node, int readCharData, int maskArrays)
{
  vtkDataSetAttributes *data = tree->GetVertexData();
  vtkStringArray *nameArr = vtkStringArray::SafeDownCast(data->GetAbstractArray(vtkXMLTreeReader::TagNameField));
  vtkStdString content;
  for (xmlNode *curNode = node; curNode; curNode = curNode->next)
    {
    //cerr << "type=" << curNode->type << ",name=" << curNode->name << endl;
    if (curNode->content)
      {
      const char *curContent = reinterpret_cast<const char*>(curNode->content);
      content += curContent;
      //cerr << "content=" << curNode->content << endl;
      }

    if (curNode->type != XML_ELEMENT_NODE)
      {
      continue;
      }

    vtkIdType vertex = -1;
    vertex = tree->AddVertex();
    if (parent != -1)
      {
      tree->AddEdge(parent, vertex);
      }

    // Append the node tag and character data to the vtkPointData
    nameArr->InsertValue(vertex, reinterpret_cast<const char*>(curNode->name));
    
    // Add pedigree ids
    vtkIdTypeArray *idArr = vtkIdTypeArray::SafeDownCast(data->GetAbstractArray("PedigreeVertexId"));
    idArr->InsertValue(vertex, vertex);
    
    // Append the element attributes to the vtkPointData
    for (xmlAttr *curAttr = curNode->properties; curAttr; curAttr = curAttr->next)
      {
      const char *name = reinterpret_cast<const char*>(curAttr->name);
      int len = static_cast<int>(strlen(name));
      char *validName = new char[len+8];
      strcpy(validName, ".valid.");
      strcat(validName, name);
      vtkStringArray *stringArr = vtkStringArray::SafeDownCast(data->GetAbstractArray(name));
      vtkBitArray *bitArr = 0;
      if (maskArrays)
        {
        bitArr = vtkBitArray::SafeDownCast(data->GetAbstractArray(validName));
        }
      if (!stringArr)
        {
        stringArr = vtkStringArray::New();
        stringArr->SetName(name);
        data->AddArray(stringArr);
        stringArr->Delete();
        if (maskArrays)
          {
          bitArr = vtkBitArray::New();
          bitArr->SetName(validName);
          data->AddArray(bitArr);
          bitArr->Delete();
          }
        }
      const char *value = reinterpret_cast<const char*>(curAttr->children->content);
      stringArr->InsertValue(vertex, value);
      if (maskArrays)
        {
        for (vtkIdType i = bitArr->GetNumberOfTuples(); i < vertex; i++)
          {
          bitArr->InsertNextValue(false);
          }
        bitArr->InsertNextValue(true);
        }
      //cerr << "attname=" << name << ",value=" << value << endl;
      delete [] validName;
      }

    // Process this node's children
    vtkXMLTreeReaderProcessElement(tree, vertex, curNode->children, readCharData, maskArrays);
    }

  if (readCharData && parent >= 0)
    {
    vtkStringArray *charArr = vtkStringArray::SafeDownCast(data->GetAbstractArray(vtkXMLTreeReader::CharDataField));
    charArr->InsertValue(parent, content);
    }
}

int vtkXMLTreeReader::RequestData(
  vtkInformation*, 
  vtkInformationVector**, 
  vtkInformationVector *outputVector)
{
  if (!this->FileName && !this->XMLString)
    {
    vtkErrorMacro("A FileName or XMLString must be specified");
    return 0;
    }

  xmlDoc *doc = NULL;
  if (this->FileName)
    {
    // Parse the file and get the DOM
    doc = xmlReadFile(this->FileName, NULL, 0);
    }
  else if (this->XMLString)
    {
    // Parse from memory and get the DOM
    doc = xmlReadMemory(this->XMLString, static_cast<int>(strlen(this->XMLString)), "noname.xml", NULL, 0);
    }

  // Store the XML hierarchy into a vtkMutableDirectedGraph,
  // later to be placed in a vtkTree.
  vtkSmartPointer<vtkMutableDirectedGraph> builder = 
    vtkSmartPointer<vtkMutableDirectedGraph>::New();

  vtkDataSetAttributes *data = builder->GetVertexData();
  vtkStringArray *nameArr = vtkStringArray::New();
  nameArr->SetName(vtkXMLTreeReader::TagNameField);
  data->AddArray(nameArr);
  nameArr->Delete();

  if (this->ReadCharData)
    {
    vtkStringArray *charArr = vtkStringArray::New();
    charArr->SetName(vtkXMLTreeReader::CharDataField);
    data->AddArray(charArr);
    charArr->Delete();
    }
  
  // Add pedigree id array
  vtkIdTypeArray *idArr = vtkIdTypeArray::New();
  idArr->SetName("PedigreeVertexId");
  data->AddArray(idArr);
  idArr->Delete();

  // Get the root element node
  xmlNode *rootElement = xmlDocGetRootElement(doc);
  vtkXMLTreeReaderProcessElement(builder, -1, rootElement, this->ReadCharData, this->MaskArrays);

  // Make all the arrays the same size
  for (int i = 0; i < data->GetNumberOfArrays(); i++)
    {
    vtkStringArray *stringArr = vtkStringArray::SafeDownCast(data->GetAbstractArray(i));
    if (stringArr && (stringArr->GetNumberOfTuples() < builder->GetNumberOfVertices()))
      {
      stringArr->InsertValue(builder->GetNumberOfVertices() - 1, vtkStdString(""));
      }
    }

  // Move the XML hierarchy into a vtkTree
  vtkTree *output = vtkTree::GetData(outputVector);
  if (!output->CheckedShallowCopy(builder))
    {
    vtkErrorMacro(<<"Structure is not a valid tree!");
    return 0;
    }

  return 1;
}

