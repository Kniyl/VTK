/* 
 * tkAppInit.c --
 *
 *	Provides a default version of the Tcl_AppInit procedure for
 *	use in wish and similar Tk-based applications.
 *
 * Copyright (c) 1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include "tk.h"

#ifdef USE_TIX
#include "tix.h"
#endif

#include <iostream.h>

/*
 *----------------------------------------------------------------------
 *
 * main --
 *
 *	This is the main program for the application.
 *
 * Results:
 *	None: Tk_Main never returns here, so this procedure never
 *	returns either.
 *
 * Side effects:
 *	Whatever the application does.
 *
 *----------------------------------------------------------------------
 */
#if (TK_MAJOR_VERSION == 3)

EXTERN int main _ANSI_ARGS_((int     argc,
                             char  **argv));
int (*tclDummyMainPtr)() = (int (*)()) main;

#if defined(DOMAIN) && defined(SING)
EXTERN "C" int matherr _ANSI_ARGS_((struct exception *));
int (*tclDummyMathPtr)() = (int (*)()) matherr;
#endif

#else
/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl.
 */
extern "C" int matherr();
int *tclDummyMathPtr = (int *) matherr;


int
main(int argc, char **argv)
{
  ios::sync_with_stdio();
  Tk_Main(argc, argv, Tcl_AppInit);
  return 0;			/* Needed only to prevent compiler warning. */
}

#endif

/*
 *----------------------------------------------------------------------
 *
 * Tcl_AppInit --
 *
 *	This procedure performs application-specific initialization.
 *	Most applications, especially those that incorporate additional
 *	packages, will have their own version of this procedure.
 *
 * Results:
 *	Returns a standard Tcl completion code, and leaves an error
 *	message in interp->result if an error occurs.
 *
 * Side effects:
 *	Depends on the startup script.
 *
 *----------------------------------------------------------------------
 */

extern "C" int Vtkcommontcl_Init(Tcl_Interp *interp);

#ifdef VTK_USE_GRAPHICS
extern "C" int Vtkgraphicstcl_Init(Tcl_Interp *interp);
#ifdef VTK_USE_TKWIDGET
extern "C" int Vtktkrenderwidget_Init(Tcl_Interp *interp);
#endif
#endif

#ifdef VTK_USE_IMAGING
extern "C" int Vtkimagingtcl_Init(Tcl_Interp *interp);
#ifdef VTK_USE_TKWIDGET
extern "C" int Vtktkimageviewerwidget_Init(Tcl_Interp *interp);
extern "C" int Vtktkimagewindowwidget_Init(Tcl_Interp *interp);
#endif
#endif

#ifdef VTK_USE_PATENTED
extern "C" int Vtkpatentedtcl_Init(Tcl_Interp *interp);
#endif

#ifdef VTK_USE_CONTRIB
extern "C" int Vtkcontribtcl_Init(Tcl_Interp *interp);
#endif

#ifdef VTK_USE_LOCAL
extern "C" int Vtklocaltcl_Init(Tcl_Interp *interp);
#endif

int Tcl_AppInit(Tcl_Interp *interp)
{
  Tk_Window main;
  
  if (Tcl_Init(interp) == TCL_ERROR) {
  return TCL_ERROR;
  }
  if (Tk_Init(interp) == TCL_ERROR) {
  return TCL_ERROR;
  }
#ifdef USE_TIX
  if (Tix_Init(interp) == TCL_ERROR) {
  return TCL_ERROR;
  }
#endif

  /* init the core vtk stuff */
  if (Vtkcommontcl_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
    
#ifdef VTK_USE_GRAPHICS
  if (Vtkgraphicstcl_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#ifdef VTK_USE_TKWIDGET
  if (Vtktkrenderwidget_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#endif
#endif

#ifdef VTK_USE_IMAGING
  if (Vtkimagingtcl_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#ifdef VTK_USE_TKWIDGET
  if (Vtktkimagewindowwidget_Init(interp) == TCL_ERROR)
    {
    return TCL_ERROR;
    }
  if (Vtktkimageviewerwidget_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#endif
#endif

#ifdef VTK_USE_PATENTED
  if (Vtkpatentedtcl_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#endif

#ifdef VTK_USE_CONTRIB
  if (Vtkcontribtcl_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#endif

#ifdef VTK_USE_LOCAL
  if (Vtklocaltcl_Init(interp) == TCL_ERROR) 
    {
    return TCL_ERROR;
    }
#endif

  /*
   * Specify a user-specific startup file to invoke if the application
   * is run interactively.  Typically the startup file is "~/.apprc"
   * where "app" is the name of the application.  If this line is deleted
   * then no user-specific startup file will be run under any conditions.
   */
  
#if (((TK_MAJOR_VERSION == 4)&&(TK_MINOR_VERSION >= 1))||((TK_MAJOR_VERSION == 8)&&(TK_MINOR_VERSION >= 0)))
    Tcl_SetVar(interp, "tcl_rcFileName", "~/.vtkrc", TCL_GLOBAL_ONLY);
#else
    tcl_RcFileName = "~/.vtkrc";
#endif
    return TCL_OK;
}






