
# - This module locates the Matlab Engine and Matlab mex.
# Defines the following:
# MATLAB_ROOT_DIR - Matlab installation directory
# MATLAB_INCLUDE_DIR - Path to Matlab include directory
# MATLAB_LIB_DIR - Path to Matlab library directory
# MATLAB_MEX_EXECUTABLE - Path to mex compiler executable
# MATLAB_EXECUTABLE - Path to the matlab executable

# Determine platform

IF(UNIX)
  IF(APPLE)
   SET(MATLAB_ROOT_SEARCH "/Applications")
   IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(MATLAB_DIR_PREFIX "maci")
   ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(MATLAB_DIR_PREFIX "maci64")
   ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
  ELSE(APPLE)
   SET(MATLAB_ROOT_SEARCH "/usr/local/")
   IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(MATLAB_DIR_PREFIX "glnx86")
   ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(MATLAB_DIR_PREFIX "glnxa64")
   ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
  ENDIF(APPLE)
ELSE(UNIX)
   SET(MATLAB_ROOT_SEARCH "C:/PROGRA~1/")
   IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(MATLAB_DIR_PREFIX "win32")
   ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(MATLAB_DIR_PREFIX "win64")
   ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
ENDIF(UNIX)

FIND_PROGRAM(MATLAB_EXECUTABLE matlab PATHS ${MATLAB_ROOT_SEARCH})

IF(NOT MATLAB_EXECUTABLE)
  MESSAGE( FATAL_ERROR "Matlab program not found, please specify MATLAB_EXECUTABLE" )
ENDIF(NOT MATLAB_EXECUTABLE)

GET_FILENAME_COMPONENT(MATLAB_ROOT_DIR ${MATLAB_EXECUTABLE} PATH)

IF("${MATLAB_ROOT_DIR}" MATCHES "/${MATLAB_DIR_PREFIX}$")
  GET_FILENAME_COMPONENT(MATLAB_ROOT_DIR ${MATLAB_ROOT_DIR} PATH)
ENDIF("${MATLAB_ROOT_DIR}" MATCHES "/${MATLAB_DIR_PREFIX}$")

IF("${MATLAB_ROOT_DIR}" MATCHES "/bin$")
  GET_FILENAME_COMPONENT(MATLAB_ROOT_DIR ${MATLAB_ROOT_DIR} PATH)
ENDIF("${MATLAB_ROOT_DIR}" MATCHES "/bin$")

SET(MATLAB_INCLUDE_DIR ${MATLAB_ROOT_DIR}/extern/include CACHE FILEPATH DOCSTRING)

IF(NOT WIN32)
  SET(MATLAB_LIB_DIR ${MATLAB_ROOT_DIR}/bin/${MATLAB_DIR_PREFIX} CACHE FILEPATH DOCSTRING)
ELSE(NOT WIN32)
  SET(MATLAB_LIB_DIR ${MATLAB_ROOT_DIR}/extern/lib/${MATLAB_DIR_PREFIX}/microsoft CACHE FILEPATH DOCSTRING)
ENDIF(NOT WIN32)

# Find mex executable 

IF(NOT WIN32)
  FIND_PROGRAM(MATLAB_MEX_EXECUTABLE mex ${MATLAB_ROOT_DIR} DOC "Path to Matlab mex compiler")  
ELSE(NOT WIN32)
  FIND_FILE(MATLAB_MEX_EXECUTABLE mex.bat ${MATLAB_ROOT_DIR} DOC "Path to Matlab mex compiler")
ENDIF(NOT WIN32)

SET(MATLAB_MEX_USE_FILE ${VTK_SOURCE_DIR}/CMake/UseMatlabMex.cmake)

