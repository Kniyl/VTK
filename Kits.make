#
# Makefile for Visualization Toolkit sources. 
# 
#------------------------------------------------------------------------------
#

SHELL = /bin/sh
.SUFFIXES: .cxx .java .class

#------------------------------------------------------------------------------

CC_FLAGS = ${CPPFLAGS} ${USER_CFLAGS} ${CFLAGS} ${USE_TOOLKIT_FLAGS} \
	   ${GRAPHICS_API_FLAGS} ${CONTROLLER_API_FLAGS} \
	   ${JAVA_INCLUDES} ${PYTHON_INCLUDES}

CXX_FLAGS = ${CPPFLAGS} ${USER_CXXFLAGS} ${CXXFLAGS} -I${srcdir} \
	${KIT_FLAGS} -I. ${USE_TOOLKIT_FLAGS} \
	${GRAPHICS_API_FLAGS} ${CONTROLLER_API_FLAGS} \
	 -I${srcdir}/../common ${TK_INCLUDE} ${TCL_INCLUDE} \
	${JAVA_INCLUDES} ${PYTHON_INCLUDES}

all: ${VTK_LIB_FILE} ${BUILD_TCL} ${BUILD_JAVA} ${BUILD_PYTHON}

.c.o:
	${CC} ${CC_FLAGS} -c $< -o $@
.cxx.o:
	${CXX} ${CXX_FLAGS} -c $< -o $@

#------------------------------------------------------------------------------
../targets:	
	cd ..; ${MAKE} targets

depend: ../targets
	../targets ${srcdir}/.. extra ${srcdir} ${KIT_EXTRA_DEPENDS} concrete $(CONCRETE) abstract $(ABSTRACT) concrete_h $(CONCRETE_H) abstract_h $(ABSTRACT_H)


targets.make: ../targets Makefile
	../targets ${srcdir}/.. extra ${srcdir} ${KIT_EXTRA_DEPENDS} concrete $(CONCRETE) abstract $(ABSTRACT) concrete_h $(CONCRETE_H) abstract_h $(ABSTRACT_H)

#------------------------------------------------------------------------------
# rules for the normal library
#
libVTK${ME}.a: ${SRC_OBJ} ${KIT_OBJ}
	${AR} cr libVTK${ME}.a ${KIT_OBJ}
	${RANLIB} libVTK$(ME).a


libVTK$(ME)$(SHLIB_SUFFIX): ${KIT_OBJ}
	rm -f libVTK$(ME)$(SHLIB_SUFFIX)
	$(CXX) ${CXX_FLAGS} ${VTK_SHLIB_BUILD_FLAGS} -o \
	libVTK$(ME)$(SHLIB_SUFFIX) \
	   ${KIT_OBJ} ${SHLIB_LD_LIBS}

#------------------------------------------------------------------------------
# rules for the tcl library
#
build_tcl: ${TCL_LIB_FILE}

tcl/${ME}Init.cxx: ../wrap/vtkWrapTclInit ${KIT_NEWS} Makefile
	../wrap/vtkWrapTclInit VTK${ME}Tcl ${KIT_NEWS} > tcl/${ME}Init.cxx

libVTK${ME}Tcl.a: tcl/${ME}Init.o ${KIT_LIBS} ${KIT_TCL_OBJ} 
	${AR} cr libVTK${ME}Tcl.a tcl/${ME}Init.o ${KIT_LIBS} ${KIT_TCL_OBJ}
	${RANLIB} libVTK$(ME)Tcl.a

libVTK$(ME)Tcl$(SHLIB_SUFFIX): tcl/${ME}Init.o ${KIT_LIBS} ${KIT_TCL_OBJ}
	rm -f libVTK$(ME)Tcl$(SHLIB_SUFFIX)
	$(CXX) ${CXX_FLAGS} ${VTK_SHLIB_BUILD_FLAGS} -o \
	libVTK$(ME)Tcl$(SHLIB_SUFFIX) \
	tcl/${ME}Init.o ${KIT_LIBS} ${KIT_TCL_OBJ}

#------------------------------------------------------------------------------
# rules for the java library
#
build_java: ${JAVA_CLASSES} ${JAVA_CODE} ${JAVA_CODE_ADD} ${JAVA_O_ADD} ${JAVA_WRAP} libVTK${ME}Java${SHLIB_SUFFIX}

.java.class:
	${JAVAC} -d ${JAVA_CLASS_HOME} -classpath ${JAVA_CLASS_HOME} $< 

libVTK$(ME)Java$(SHLIB_SUFFIX): ${JAVA_O_ADD} ${JAVA_WRAP}
	rm -f libVTK$(ME)Java$(SHLIB_SUFFIX)
	$(CXX) ${CXX_FLAGS} ${VTK_SHLIB_BUILD_FLAGS} \
	-o libVTK$(ME)Java$(SHLIB_SUFFIX) \
	${JAVA_O_ADD} ${JAVA_WRAP} 

#------------------------------------------------------------------------------
# rules for the python library
#
build_python: ${PYTHON_O_ADD} ${PYTHON_WRAP} libVTK${ME}Python${SHLIB_SUFFIX} 

python/${ME}Init.cxx: ../wrap/vtkWrapPythonInit ${KIT_NEWS} Makefile
	../wrap/vtkWrapPythonInit libVTK${ME}Python ${KIT_NEWS} > python/${ME}Init.cxx

libVTK$(ME)Python$(SHLIB_SUFFIX): python/${ME}Init.o ${KIT_OBJ} \
	${PYTHON_O_ADD} ${PYTHON_WRAP}
	rm -f libVTK$(ME)Python$(SHLIB_SUFFIX)
	$(CXX) ${CXX_FLAGS} ${VTK_SHLIB_BUILD_FLAGS} \
	-o libVTK$(ME)Python$(SHLIB_SUFFIX) python/${ME}Init.o \
	  ${KIT_OBJ} ${PYTHON_O_ADD} ${PYTHON_WRAP} \
	-L. ${XLDFLAGS} ${PYTHON_LIBS} \
	${XLIBS} -lXext -lXt ${X_PRE_LIBS} -lX11 ${X_EXTRA_LIBS} ${DL_LIBS} \
	${THREAD_LIBS}


#------------------------------------------------------------------------------
clean: ${CLEAN_TCL} $(CLEAN_JAVA) $(CLEAN_PYTHON)
	-rm -f *.o *.a *.so *.sl *~ Makefile

clean_tcl:
	-cd tcl; rm -f *

clean_java:
	-cd java; rm -f *

clean_python:
	-cd python; rm -f *

#------------------------------------------------------------------------------
install_python: install build_python
	@echo "Installing libVTK${ME}Python${SHLIB_SUFFIX}"
	${INSTALL} -m 755 libVTK${ME}Python${SHLIB_SUFFIX} $(LIB_INSTALL_DIR)/libVTK${ME}Python${SHLIB_SUFFIX}

install_java: install build_java
	@echo "Installing libVTK${ME}Java${SHLIB_SUFFIX}"
	${INSTALL} -m 755 libVTK${ME}Java${SHLIB_SUFFIX} $(LIB_INSTALL_DIR)/libVTK${ME}Java${SHLIB_SUFFIX}

install_tcl: install ${TCL_LIB_FILE}
	@echo "Installing ${TCL_LIB_FILE}"
	${INSTALL} -m 755 $(TCL_LIB_FILE) $(LIB_INSTALL_DIR)/$(TCL_LIB_FILE)

install: ${VTK_LIB_FILE} 
	@echo "Installing ${VTK_LIB_FILE}"
	${INSTALL} -m 755 $(VTK_LIB_FILE) $(LIB_INSTALL_DIR)/$(VTK_LIB_FILE)
