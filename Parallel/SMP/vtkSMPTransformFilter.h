#ifndef __vtkSMPTransformFilter_h__
#define __vtkSMPTransformFilter_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkTransformFilter.h"

class vtkInformation;
class vtkInformationVector;

class VTKPARALLELSMP_EXPORT vtkSMPTransformFilter : public vtkTransformFilter
{
public :
  vtkTypeMacro(vtkSMPTransformFilter,vtkTransformFilter);
  static vtkSMPTransformFilter *New();
  void PrintSelf(ostream& os, vtkIndent indent);

protected :
  vtkSMPTransformFilter();
  ~vtkSMPTransformFilter();

  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *);

private :
  vtkSMPTransformFilter(const vtkSMPTransformFilter&);  // Not implemented.
  void operator=(const vtkSMPTransformFilter&);  // Not implemented.

};

#endif //__vtkSMPTransformFilter_h__
