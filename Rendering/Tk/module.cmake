vtk_module(vtkRenderingTk
  GROUPS
    Tk
  DEPENDS
    vtkRenderingOpenGL
    vtkInteractionStyle
    vtkInteractionImage
  COMPILE_DEPENDS
    vtkTclTk
  TEST_DEPENDS
    vtkRenderingVolume
  EXCLUDE_FROM_WRAPPING
  )
