/* this is a CMake loadable command to wrap vtk objects into Python */

#include "cmCPluginAPI.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct 
{
  char *LibraryName;
  int NumberWrapped;
  void **SourceFiles;
  char **HeaderFiles;
} cmVTKWrapPythonData;

/* this roputine creates the init file */
static void CreateInitFile(cmLoadedCommandInfo *info,
                           void *mf, const char *kitName, 
                           int numClasses, const char **classes) 
{
  /* we have to make sure that the name is the correct case */
  int i;
  char *tempOutputFile;  
  char *outFileName = 
    (char *)malloc(strlen(info->CAPI->GetCurrentOutputDirectory(mf)) + 
                   strlen(kitName) + 10);
  FILE *fout;
  
  sprintf(outFileName,"%s/%sInit.cxx",
          info->CAPI->GetCurrentOutputDirectory(mf), kitName);
  
  tempOutputFile = (char *)malloc(strlen(outFileName) + 5);
  sprintf(tempOutputFile,"%s.tmp",outFileName);
  fout = fopen(tempOutputFile,"w");
  if (!fout)
    {
    return;
    }
  
  fprintf(fout,"// Generated by cmVTKWrapPythonCommand2 in VTK/CMake\n\n");
  fprintf(fout,"#include \"vtkSystemIncludes.h\"\n");
  fprintf(fout,"#include <string.h>\n");
  fprintf(fout,"#include \"Python.h\"\n\n");
  fprintf(fout,"// Handle compiler warning messages, etc.\n"
          "#if defined( _MSC_VER ) && !defined(VTK_DISPLAY_WIN32_WARNINGS)\n"
          "#pragma warning ( disable : 4706 )\n"
          "#endif // Windows Warnings \n\n");
  
  for (i = 0; i < numClasses; i++)
    {
#ifdef _WIN32
    fprintf(fout,"extern  \"C\" {__declspec( dllexport) PyObject *PyVTKClass_%sNew(char *); }\n",classes[i]);
#else
    fprintf(fout,"extern  \"C\" {PyObject *PyVTKClass_%sNew(char *); }\n",classes[i]);
#endif
    }
  
  fprintf(fout,"\nstatic PyMethodDef Py%s_ClassMethods[] = {\n",
          kitName);
  fprintf(fout,"{NULL, NULL, 0, NULL}};\n\n");
  
#ifdef _WIN32
  fprintf(fout,"extern  \"C\" {__declspec( dllexport) void init%s();}\n\n",kitName);
  fprintf(fout,"void init%s()\n{\n",kitName);
#else
  fprintf(fout,"extern  \"C\" {void initlib%s();}\n\n",kitName);
  fprintf(fout,"void initlib%s()\n{\n",kitName);
#endif
  
  
  /* module init function */
  fprintf(fout,"  PyObject *m, *d, *c;\n\n");
#ifdef _WIN32
  fprintf(fout,"  static const char modulename[] = \"%s\";\n",kitName);
#else
  fprintf(fout,"  static const char modulename[] = \"lib%s\";\n",kitName);
#endif
  fprintf(fout,"  m = Py_InitModule((char*)modulename, Py%s_ClassMethods);\n",
    kitName);
  
  fprintf(fout,"  d = PyModule_GetDict(m);\n");
  fprintf(fout,"  if (!d) Py_FatalError((char*)\"can't get dictionary for module %s!\");\n\n",
          kitName);
  
  for (i = 0; i < numClasses; i++)
    {
    fprintf(fout,"  if ((c = PyVTKClass_%sNew((char*)modulename)))\n",
            classes[i]);
    fprintf(fout,"    if (-1 == PyDict_SetItemString(d, (char*)\"%s\", c))\n",
      classes[i]);
    fprintf(fout,"      Py_FatalError((char*)\"can't add class %s to dictionary!\");\n\n",
            classes[i]);
    }
  fprintf(fout,"}\n\n");
  fclose(fout);

  /* copy the file if different */
  info->CAPI->CopyFileIfDifferent(tempOutputFile, outFileName);
  info->CAPI->RemoveFile(tempOutputFile);
}

/* do almost everything in the initial pass */
static int InitialPass(void *inf, void *mf, int argc, char *argv[])
{
  cmLoadedCommandInfo *info = (cmLoadedCommandInfo *)inf;
  int i;
  int newArgc;
  char **newArgv;
  int numClasses = 0;
  char **classes = 0;
  int numWrapped = 0;
  cmVTKWrapPythonData *cdata = 
    (cmVTKWrapPythonData *)malloc(sizeof(cmVTKWrapPythonData));
  const char *cdir = info->CAPI->GetCurrentDirectory(mf);
  const char *def = 0;
  char *sourceListValue = 0;
  char *newName;
  void *cfile = 0;

  if(argc < 3 )
    {
    info->CAPI->SetError(info, "called with incorrect number of arguments");
    return 0;
    }
  
  info->CAPI->ExpandSourceListArguments(mf, argc, argv, 
                                        &newArgc, &newArgv, 2);
  
  /* Now check and see if the value has been stored in the cache */
  /* already, if so use that value and don't look for the program */
  if(!info->CAPI->IsOn(mf,"VTK_WRAP_PYTHON"))
    {
    info->CAPI->FreeArguments(newArgc, newArgv);
    return 1;
    }

  /* keep the library name */
  classes = (char **)malloc(sizeof(char *)*newArgc);
  cdata->LibraryName = strdup(newArgv[0]);
  cdata->SourceFiles = (void **)malloc(sizeof(void *)*newArgc);
  cdata->HeaderFiles = (char **)malloc(sizeof(char *)*newArgc);

  /* was the list already populated */
  def = info->CAPI->GetDefinition(mf, newArgv[1]);
  sourceListValue = 
    (char *)malloc(info->CAPI->GetTotalArgumentSize(newArgc,newArgv)+
                   newArgc*14 + (def ? strlen(def) : 0) + 10);
  if (def)
    {
    sprintf(sourceListValue,"%s;%sInit.cxx",def,newArgv[0]);
    }
  else
    {
      /* Don't include the Init.cxx file in the library on OSX */
      /* It is linked into a MODULE separate from the rest of the dylib */
#if defined(__APPLE__)
    sprintf(sourceListValue,"");
#else
    sprintf(sourceListValue,"%sInit.cxx",newArgv[0]);
#endif
    }

  /* get the classes for this lib */
  for(i = 2; i < newArgc; ++i)
    {   
    void *curr = info->CAPI->GetSource(mf,newArgv[i]);
    
    /* if we should wrap the class */
    if (!curr || 
        !info->CAPI->SourceFileGetPropertyAsBool(curr,"WRAP_EXCLUDE"))
      {
      void *file = info->CAPI->CreateSourceFile();
      char *srcName;
      char *hname;
      char *pathName;
      srcName = info->CAPI->GetFilenameWithoutExtension(newArgv[i]);
      pathName = info->CAPI->GetFilenamePath(newArgv[i]);
      if (curr)
        {
        int abst = info->CAPI->SourceFileGetPropertyAsBool(curr,"ABSTRACT");
        info->CAPI->SourceFileSetProperty(file,"ABSTRACT",
                                          (abst ? "1" : "0"));
        }
      classes[numClasses] = strdup(srcName);
      numClasses++;
      newName = (char *)malloc(strlen(srcName)+7);
      sprintf(newName,"%sPython",srcName);
      info->CAPI->SourceFileSetName2(file, newName, 
                                     info->CAPI->GetCurrentOutputDirectory(mf),
                                     "cxx",0);

      if (strlen(pathName) > 1)
        {
        hname = (char *)malloc(strlen(pathName) + strlen(srcName) + 4);
        sprintf(hname,"%s/%s.h",pathName,srcName);
        }
      else
        {
        hname = (char *)malloc(strlen(cdir) + strlen(srcName) + 4);
        sprintf(hname,"%s/%s.h",cdir,srcName);
        }
      /* add starting depends */
      info->CAPI->SourceFileAddDepend(file,hname);
      info->CAPI->AddSource(mf,file);
      cdata->SourceFiles[numWrapped] = file;
      cdata->HeaderFiles[numWrapped] = hname;
      numWrapped++;
      if(sourceListValue[0])
        {
        /* This is not the first value, add a separator. */
        strcat(sourceListValue,";");
        }
      strcat(sourceListValue,newName);
      strcat(sourceListValue,".cxx");        
      free(newName);
      info->CAPI->Free(srcName);
      info->CAPI->Free(pathName);
      }
    }
  
  /* add the init file */
  cfile = info->CAPI->CreateSourceFile();
  info->CAPI->SourceFileSetProperty(cfile,"ABSTRACT","0");
  newName = (char *)malloc(strlen(newArgv[0]) + 5);
  sprintf(newName,"%sInit",newArgv[0]);
  CreateInitFile(info,mf,newArgv[0],numClasses,classes);
  info->CAPI->SourceFileSetName2(cfile, newName, 
                                 info->CAPI->GetCurrentOutputDirectory(mf),
                                 "cxx",0);
  free(newName);
  info->CAPI->AddSource(mf,cfile);

  cdata->NumberWrapped = numWrapped;
  info->CAPI->SetClientData(info,cdata);

  info->CAPI->AddDefinition(mf, newArgv[1], sourceListValue);
  info->CAPI->FreeArguments(newArgc, newArgv);
  free(sourceListValue);
  return 1;
}
  
  
static void FinalPass(void *inf, void *mf) 
{
  cmLoadedCommandInfo *info = (cmLoadedCommandInfo *)inf;
  /* get our client data from initial pass */
  cmVTKWrapPythonData *cdata = 
    (cmVTKWrapPythonData *)info->CAPI->GetClientData(info);
  
  /* first we add the rules for all the .h to Python.cxx files */
  const char *wpython = "${VTK_WRAP_PYTHON_EXE}";
  const char *hints = info->CAPI->GetDefinition(mf,"VTK_WRAP_HINTS");
  const char *args[4];
  const char *depends[2];
  int i;
  int numDepends, numArgs;
  const char *cdir = info->CAPI->GetCurrentDirectory(mf);
  
  /* wrap all the .h files */
  depends[0] = wpython;
  numDepends = 1;
  if (hints)
    {
    depends[1] = hints;
    numDepends++;
    }
  for(i = 0; i < cdata->NumberWrapped; i++)
    {
    char *res;
    const char *srcName = info->CAPI->SourceFileGetSourceName(cdata->SourceFiles[i]);
    args[0] = cdata->HeaderFiles[i];
    numArgs = 1;
    if (hints)
      {
      args[1] = hints;
      numArgs++;
      }
    args[numArgs] = 
      (info->CAPI->SourceFileGetPropertyAsBool(cdata->SourceFiles[i],"ABSTRACT") ?"0" :"1");
    numArgs++;
    res = (char *)malloc(strlen(info->CAPI->GetCurrentOutputDirectory(mf)) + 
                         strlen(srcName) + 6);
    sprintf(res,"%s/%s.cxx",info->CAPI->GetCurrentOutputDirectory(mf),srcName);
    args[numArgs] = res;
    numArgs++;
    info->CAPI->AddCustomCommand(mf, args[0],
                       wpython, numArgs, args, numDepends, depends, 
                       1, &res, cdata->LibraryName);
    free(res);
    }
}

static void Destructor(void *inf) 
{
  int i;
  cmLoadedCommandInfo *info = (cmLoadedCommandInfo *)inf;
  /* get our client data from initial pass */
  cmVTKWrapPythonData *cdata = 
    (cmVTKWrapPythonData *)info->CAPI->GetClientData(info);
  if (cdata)
    {
    for (i = 0; i < cdata->NumberWrapped; ++i)
      {              
      info->CAPI->DestroySourceFile(cdata->SourceFiles[i]);
      free(cdata->HeaderFiles[i]);
      }
    free(cdata->SourceFiles);
    free(cdata->HeaderFiles);
    free(cdata->LibraryName);
    free(cdata);
    }
}

static const char* GetTerseDocumentation() 
{
  return "Create Python Wrappers.";
}

static const char* GetFullDocumentation()
{
  return
    "VTK_WRAP_PYTHON(resultingLibraryName SourceListName SourceLists ...)";
}

void CM_PLUGIN_EXPORT VTK_WRAP_PYTHON2Init(cmLoadedCommandInfo *info)
{
  info->InitialPass = InitialPass;
  info->FinalPass = FinalPass;
  info->Destructor = Destructor;
  info->m_Inherited = 0;
  info->GetTerseDocumentation = GetTerseDocumentation;
  info->GetFullDocumentation = GetFullDocumentation;  
  info->Name = "VTK_WRAP_PYTHON2";
}




