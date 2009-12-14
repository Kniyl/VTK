
#
# - This module locates an installed R distribution.  
#
# Defines the following:
#
# R_ROOT_DIR - R installation directory
# R_INCLUDE_DIR - Path to R include directory
# R_LIB_DIR - Path to R library directory
# R_COMMAND - Path to R command
#

FIND_FILE(R_HEADER "R.h" PATHS "c:/R-2.9.0/include" "/Applications"
                               "/Frameworks"
                               "/usr/local"
                               "/usr/include/R")

IF(NOT R_HEADER)
 MESSAGE( FATAL_ERROR "R.h header file not found, please specify path in R_HEADER" )
ENDIF(NOT R_HEADER)

# Determine platform

IF(UNIX)
  IF(APPLE)
   GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_HEADER} PATH)
   GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_ROOT_DIR} PATH)
   IF("${R_ROOT_DIR}" MATCHES "/Resources$")
    GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_ROOT_DIR} PATH)
   ENDIF("${R_ROOT_DIR}" MATCHES "/Resources$")
   SET(R_COMMAND ${R_ROOT_DIR}/Resources/bin/R CACHE FILEPATH DOCSTRING)
   SET(R_INCLUDE_DIR ${R_ROOT_DIR}/Resources/include CACHE FILEPATH DOCSTRING)
   SET(R_LIB_DIR ${R_ROOT_DIR}/Resources/lib CACHE FILEPATH DOCSTRING)
  ELSE(APPLE)
   IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(BIN_DIR_PREFIX "lib")
   ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(BIN_DIR_PREFIX "lib64")
   ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
   GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_HEADER} PATH)
   GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_ROOT_DIR} PATH)
   GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_ROOT_DIR} PATH)
   SET(R_COMMAND ${R_ROOT_DIR}/bin/R CACHE FILEPATH DOCSTRING)
   SET(R_INCLUDE_DIR ${R_ROOT_DIR}/include/R CACHE FILEPATH DOCSTRING)
   SET(R_LIB_DIR ${R_ROOT_DIR}/${BIN_DIR_PREFIX}/R/lib CACHE FILEPATH DOCSTRING)
  ENDIF(APPLE)
ELSE(UNIX)
   GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_HEADER} PATH)
   IF("${R_ROOT_DIR}" MATCHES "/include$")
     GET_FILENAME_COMPONENT(R_ROOT_DIR ${R_ROOT_DIR} PATH)
   ENDIF("${R_ROOT_DIR}" MATCHES "/include$")
   SET(R_COMMAND ${R_ROOT_DIR}/bin/R CACHE FILEPATH DOCSTRING)
   SET(R_INCLUDE_DIR "${R_ROOT_DIR}/include" CACHE FILEPATH DOCSTRING)
   SET(R_LIB_DIR ${R_ROOT_DIR}/bin CACHE FILEPATH DOCSTRING)
ENDIF(UNIX)

SET(R_USE_FILE ${VTK_SOURCE_DIR}/CMake/UseR.cmake)

