/*=========================================================================*/
// .NAME vtkRayBounder
// .SECTION Description
// 
// 

// .SECTION see also
// 

#ifndef __vtkRayBounder_h
#define __vtkRayBounder_h

#include "vtkObject.h"
#include "vtkRenderer.h"

class VTK_EXPORT vtkRayBounder : public vtkObject
{
public:
  const char *GetClassName() {return "vtkRayBounder";};

  virtual float *GetRayBounds( vtkRenderer *ren )=0;

protected:
};

#endif

