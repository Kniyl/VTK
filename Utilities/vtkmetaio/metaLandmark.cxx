/*=========================================================================

  Program:   MetaIO
  Module:    metaLandmark.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "metaLandmark.h"

#include <stdio.h>
#include <ctype.h>
#include <string>

#if (METAIO_USE_NAMESPACE)
namespace METAIO_NAMESPACE {
#endif

//
// MedImage Constructors
//
MetaLandmark::
MetaLandmark()
:MetaObject()
{
  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark()" << METAIO_STREAM::endl;
  m_NPoints = 0;
  Clear();
}

//
MetaLandmark::
MetaLandmark(const char *_headerName)
:MetaObject()
{
  if(META_DEBUG)  METAIO_STREAM::cout << "MetaLandmark()" << METAIO_STREAM::endl;
  m_NPoints = 0;
  Clear();
  Read(_headerName);
}

//
MetaLandmark::
MetaLandmark(const MetaLandmark *_tube)
:MetaObject()
{
  if(META_DEBUG)  METAIO_STREAM::cout << "MetaLandmark()" << METAIO_STREAM::endl;
  m_NPoints = 0;
  Clear();
  CopyInfo(_tube);
}



//
MetaLandmark::
MetaLandmark(unsigned int dim)
:MetaObject(dim)
{
  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark()" << METAIO_STREAM::endl;
  m_NPoints = 0;
  Clear();
}

//
MetaLandmark::
~MetaLandmark()
{
  Clear();
  M_Destroy();
}

//
void MetaLandmark::
PrintInfo() const
{
  MetaObject::PrintInfo();
  METAIO_STREAM::cout << "PointDim = " << m_PointDim << METAIO_STREAM::endl;
  METAIO_STREAM::cout << "NPoints = " << m_NPoints << METAIO_STREAM::endl;
  char str[255];
  MET_TypeToString(m_ElementType, str);
  METAIO_STREAM::cout << "ElementType = " << str << METAIO_STREAM::endl;
}

void MetaLandmark::
CopyInfo(const MetaObject * _object)
{
  MetaObject::CopyInfo(_object);
}

    

void MetaLandmark::
PointDim(const char* pointDim)
{
  strcpy(m_PointDim,pointDim);
}
    
const char* MetaLandmark::
PointDim(void) const
{
  return m_PointDim;
}

void MetaLandmark::
NPoints(int npnt)
{
  m_NPoints = npnt;
}

int MetaLandmark::
NPoints(void) const
{
  return m_NPoints;
}


/** Clear tube information */
void MetaLandmark::
Clear(void)
{
  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark: Clear" << METAIO_STREAM::endl;
  MetaObject::Clear();
  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark: Clear: m_NPoints" << METAIO_STREAM::endl;
  // Delete the list of pointers to tubes.
  PointListType::iterator it = m_PointList.begin();
  while(it != m_PointList.end())
  {
    LandmarkPnt* pnt = *it;
    it++;
    delete pnt;
  }  
  m_PointList.clear();
  m_NPoints = 0;
  strcpy(m_PointDim, "x y z red green blue alpha");
  m_ElementType = MET_FLOAT;
}
        
/** Destroy tube information */
void MetaLandmark::
M_Destroy(void)
{
  MetaObject::M_Destroy();
}

/** Set Read fields */
void MetaLandmark::
M_SetupReadFields(void)
{
  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark: M_SetupReadFields" << METAIO_STREAM::endl;

  MetaObject::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "PointDim", MET_STRING, true);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "NPoints", MET_INT, true);
  m_Fields.push_back(mF);
 
  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "ElementType", MET_STRING, true);
  mF->required = true;
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Points", MET_NONE, true);
  mF->terminateRead = true;
  m_Fields.push_back(mF);

}

MET_ValueEnumType MetaLandmark::
ElementType(void) const
{
  return m_ElementType;
}

void MetaLandmark::
ElementType(MET_ValueEnumType _elementType)
{
  m_ElementType = _elementType;
}


void MetaLandmark::
M_SetupWriteFields(void)
{
  strcpy(m_ObjectTypeName,"Landmark");
  MetaObject::M_SetupWriteFields();

  MET_FieldRecordType * mF;

  char s[255];
  mF = new MET_FieldRecordType;
  MET_TypeToString(m_ElementType, s);
  MET_InitWriteField(mF, "ElementType", MET_STRING, strlen(s), s);
  m_Fields.push_back(mF);

  if(strlen(m_PointDim)>0)
  {
    mF = new MET_FieldRecordType;
    MET_InitWriteField(mF, "PointDim", MET_STRING,
                           strlen(m_PointDim),m_PointDim);
    m_Fields.push_back(mF);
  }

  m_NPoints = m_PointList.size();
  mF = new MET_FieldRecordType;
  MET_InitWriteField(mF, "NPoints", MET_INT,m_NPoints);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitWriteField(mF, "Points", MET_NONE);
  m_Fields.push_back(mF);

}



bool MetaLandmark::
M_Read(void)
{
  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark: M_Read: Loading Header" << METAIO_STREAM::endl;

  if(!MetaObject::M_Read())
  {
    METAIO_STREAM::cout << "MetaLandmark: M_Read: Error parsing file" << METAIO_STREAM::endl;
    return false;
  }

  if(META_DEBUG) METAIO_STREAM::cout << "MetaLandmark: M_Read: Parsing Header" << METAIO_STREAM::endl;
 
  MET_FieldRecordType * mF;
 
  mF = MET_GetFieldRecord("NPoints", &m_Fields);
  if(mF->defined)
  {
    m_NPoints= (int)mF->value[0];
  }

  mF = MET_GetFieldRecord("ElementType", &m_Fields);
  if(mF->defined)
  {
    MET_StringToType((char *)(mF->value), &m_ElementType);
  }


  mF = MET_GetFieldRecord("PointDim", &m_Fields);
  if(mF->defined)
  {
    strcpy(m_PointDim,(char *)(mF->value));
  }

  int* posDim= new int[m_NDims];
  int i;
  for(i= 0; i < m_NDims; i++)
  {
    posDim[i] = -1;
  }

  int pntDim;
  char** pntVal = NULL;
  MET_StringToWordArray(m_PointDim, &pntDim, &pntVal); 
 
    
  int j;
  for(j = 0; j < pntDim; j++) 
  {
    if(!strcmp(pntVal[j], "x") || !strcmp(pntVal[j], "X"))
    {
      posDim[0] = j;
    }
    if(!strcmp(pntVal[j], "y") || !strcmp(pntVal[j], "Y"))
    {
      posDim[1] = j;
    }
    if(!strcmp(pntVal[j], "z") || !strcmp(pntVal[j], "Z"))
    {
      posDim[2] = j;
    }

  }

  for(i=0;i<pntDim;i++)
    {
      delete [] pntVal[i];
    }
  delete [] pntVal;

  float v[16];
  
  if(m_BinaryData)
  {
    int elementSize;
    MET_SizeOfType(m_ElementType, &elementSize);
    int readSize = m_NPoints*(m_NDims+4)*elementSize;
    
    char* _data = new char[readSize];
    m_ReadStream->read((char *)_data, readSize);

    int gc = m_ReadStream->gcount();
    if(gc != readSize)
    {
      METAIO_STREAM::cout << "MetaLandmark: m_Read: data not read completely" 
                << METAIO_STREAM::endl;
      METAIO_STREAM::cout << "   ideal = " << readSize << " : actual = " << gc << METAIO_STREAM::endl;
      return false;
    }

    i=0;
    int d;
    unsigned int k;
    for(j=0; j<m_NPoints; j++) 
    {
      LandmarkPnt* pnt = new LandmarkPnt(m_NDims);

      for(d=0; d<m_NDims; d++)
      {
        char* num = new char[sizeof(float)];
        for(k=0;k<sizeof(float);k++)
          {
          num[k] = _data[i+k];
          }
        float td = (float)((float*)num)[0];
        MET_SwapByteIfSystemMSB(&td,MET_FLOAT);
        i+=sizeof(float); 
        pnt->m_X[d] = (float)td;
        delete [] num;
      }

      for(d=0; d<4; d++)
      {
        char* num = new char[sizeof(float)];
        for(k=0;k<sizeof(float);k++)
          {
          num[k] = _data[i+k];
          }
        float td = (float)((float*)num)[0];
        MET_SwapByteIfSystemMSB(&td,MET_FLOAT);
        i+=sizeof(float); 
        pnt->m_Color[d] = (float)td;
        delete [] num;
      }

      m_PointList.push_back(pnt);
    }
    delete [] _data;
  }
  else
  {
    for(j=0; j<m_NPoints; j++) 
    {
      LandmarkPnt* pnt = new LandmarkPnt(m_NDims);
      
      for(int k=0; k<pntDim; k++)
      {
        *m_ReadStream >> v[k];
        m_ReadStream->get(); 
      }

      int d;
      for(d=0; d<m_NDims; d++)
      {
        pnt->m_X[d] = v[posDim[d]];
      }

      for(d=0; d<4; d++)
      {
        pnt->m_Color[d] = v[d+m_NDims];
      }
     
      m_PointList.push_back(pnt);
    }

      
    char c = ' ';
    while( (c!='\n') && (!m_ReadStream->eof()))
    {
      c = m_ReadStream->get();// to avoid unrecognize charactere
    }
  }
  
  delete [] posDim;
  return true;
}


bool MetaLandmark::
M_Write(void)
{

  if(!MetaObject::M_Write())
  {
    METAIO_STREAM::cout << "MetaLandmark: M_Read: Error parsing file" << METAIO_STREAM::endl;
    return false;
  }

  /** Then copy all points */
  if(m_BinaryData)
  {
    PointListType::const_iterator it = m_PointList.begin();
    int elementSize;
    MET_SizeOfType(m_ElementType, &elementSize);

    char* data = new char[(m_NDims+4)*m_NPoints*elementSize];
    int i=0;
    int d;
    while(it != m_PointList.end())
    {
      for(d = 0; d < m_NDims; d++)
      {
        float x = (*it)->m_X[d];
        MET_SwapByteIfSystemMSB(&x,MET_FLOAT);
        MET_DoubleToValue((double)x,m_ElementType,data,i++);
      }

      for(d = 0; d < 4; d++)
      {
        float c = (*it)->m_Color[d];
        MET_SwapByteIfSystemMSB(&c,MET_FLOAT);
        MET_DoubleToValue((double)c,m_ElementType,data,i++);
      }
      it++;
    }  
    m_WriteStream->write((char *)data,(m_NDims+4)*m_NPoints*elementSize);
    m_WriteStream->write("\n",1);
    delete [] data;
  }
  else
  {
    PointListType::const_iterator it = m_PointList.begin();
  
    int d;
    while(it != m_PointList.end())
    {
      for(d = 0; d < m_NDims; d++)
      {
        *m_WriteStream << (*it)->m_X[d] << " ";
      }

      for(d = 0; d < 4; d++)
      {
        *m_WriteStream << (*it)->m_Color[d] << " ";
      }

      *m_WriteStream << METAIO_STREAM::endl;
      it++;
    }
  }

  return true;

}

#if (METAIO_USE_NAMESPACE)
};
#endif

