/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : Id  */
/*  Date : $Date$ */
/*  Version : $Revision$ */
/*                                                                 */
/*  Author:                                                        */
/*     Jerry A. Clarke                                             */
/*     clarke@arl.army.mil                                         */
/*     US Army Research Laboratory                                 */
/*     Aberdeen Proving Ground, MD                                 */
/*                                                                 */
/*     Copyright @ 2002 US Army Research Laboratory                */
/*     All Rights Reserved                                         */
/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/
#ifndef _vtkXdmfDataArray_h
#define _vtkXdmfDataArray_h

#include "vtkIOXdmf2Module.h" // For export macro
#include <vtkDataArray.h> // For?
#include <vtkFloatArray.h> // For?
#include <XdmfArray.h> //For?

class VTKIOXDMF2_EXPORT vtkXdmfDataArray : public vtkObject
{
public:
  static vtkXdmfDataArray *New();
  vtkTypeMacro(vtkXdmfDataArray,vtkObject);

  vtkDataArray *FromArray( void ) {
    return( this->FromXdmfArray() );
    }
  char *ToArray( void ) {
    return( this->ToXdmfArray() );
    }
  vtkDataArray *FromXdmfArray( char *ArrayName = NULL, int CopyShape = 1,
   int rank = 1, int Components = 1 , int MakeCopy = 1);
  char *ToXdmfArray( vtkDataArray *DataArray = NULL, int CopyShape = 1 );

  void SetArray( char *TagName ) {
    this->Array = TagNameToArray( TagName );
    if( this->Array ) {
      this->FromXdmfArray();
      }
    }

  char *GetArray( void ) {
    if ( this->Array != NULL ) {
      return( this->Array->GetTagName() );
    }
    return( NULL );
    }

  void SetVtkArray( vtkDataArray *array) {
    this->vtkArray = array;
    this->ToXdmfArray( array );
    }

  vtkDataArray *GetVtkArray( void ) {
    return( this->vtkArray );
    }

protected:
  vtkXdmfDataArray();

private:
  vtkDataArray  *vtkArray;
  XdmfArray  *Array;
  vtkXdmfDataArray(const vtkXdmfDataArray&); // Not implemented
  void operator=(const vtkXdmfDataArray&); // Not implemented
};

#endif /* _vtkXdmfDataArray_h */
