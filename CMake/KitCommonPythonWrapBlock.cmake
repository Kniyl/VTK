# Create custom commands to generate the python wrappers for this kit.
VTK_WRAP_PYTHON3(vtk${KIT}Python KitPython_SRCS "${Kit_SRCS}")

# Create a shared library containing the python wrappers.  Executables
# can link to this but it is not directly loaded dynamically as a
# module.
ADD_LIBRARY(vtk${KIT}PythonD ${KitPython_SRCS} ${Kit_PYTHON_EXTRA_SRCS})
TARGET_LINK_LIBRARIES(vtk${KIT}PythonD vtk${KIT} ${KIT_PYTHON_LIBS})
IF(NOT VTK_INSTALL_NO_LIBRARIES)
  INSTALL_TARGETS(${VTK_INSTALL_LIB_DIR} RUNTIME_DIRECTORY ${VTK_INSTALL_BIN_DIR} vtk${KIT}PythonD)
ENDIF(NOT VTK_INSTALL_NO_LIBRARIES)
SET(KIT_LIBRARY_TARGETS ${KIT_LIBRARY_TARGETS} vtk${KIT}PythonD)

# On some UNIX platforms the python library is static and therefore
# should not be linked into the shared library.  Instead the symbols
# are exported from the python executable so that they can be used by
# shared libraries that are linked or loaded.  On Windows and OSX we
# want to link to the python libray to resolve its symbols
# immediately.
IF(WIN32 OR APPLE)
  TARGET_LINK_LIBRARIES (vtk${KIT}PythonD ${VTK_PYTHON_LIBRARIES})
ENDIF(WIN32 OR APPLE)

# Add dependencies that may have been generated by VTK_WRAP_PYTHON3 to
# the python wrapper library.  This is needed for the
# pre-custom-command hack in Visual Studio 6.
IF(KIT_PYTHON_DEPS)
  ADD_DEPENDENCIES(vtk${KIT}PythonD ${KIT_PYTHON_DEPS})
ENDIF(KIT_PYTHON_DEPS)

# Create a python module that can be loaded dynamically.  It links to
# the shared library containing the wrappers for this kit.
ADD_LIBRARY(vtk${KIT}Python MODULE vtk${KIT}PythonInit.cxx)
TARGET_LINK_LIBRARIES(vtk${KIT}Python vtk${KIT}PythonD)

IF(WIN32 OR APPLE)
  TARGET_LINK_LIBRARIES (vtk${KIT}Python ${VTK_PYTHON_LIBRARIES})
ENDIF(WIN32 OR APPLE)
