#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_DEPENDS 2000

struct DEPENDS_STRUCT
{
	int indices[100];  // aloows up to 100 includes in a single files
	int numIndices;
	char name[256];
};

static DEPENDS_STRUCT *DependsStructArray[MAX_DEPENDS];
static int NumInDepends = 0;

static int  num;

// 1000 leaves a LOT of room to grow (8/97)
static int dependIndices[1000];

//static char depends[400][256];
//static char names[400][256];


void GetDepends(int index)
{
  int i, j;

  for (i = 0; i < DependsStructArray[index]->numIndices; i++)
    {
    // see if entry is alreay in the list
    for (j = 0; j < num; j++)
      {
      if ( DependsStructArray[index]->indices[i] == dependIndices[j] )
        break;
      }
    if (j != num ) // already been added
      continue;

    dependIndices[ num++ ] = DependsStructArray[index]->indices[i];
    GetDepends(DependsStructArray[index]->indices[i]);
    }
}


extern void OutputUNIXDepends(char *file, FILE *fp)
{
  int i;

  num = 0;

  // find this entry in DependsStructArray
  for (i = 0; i < NumInDepends; i++)
    {
    if ( !strcmp(file,DependsStructArray[i]->name) )
      break;
    }

  if ( i == NumInDepends )
    {
    fprintf(stderr,"Bad Error!!! Could not find file %s to generate depends on.!!!!\nConsider checking the Makefile and Makefile.in for a bad file name.\nIf special file (not in Makefile.in) may need to add to SetupDepends() \nin targets.cxx\n",file);
    exit(-1);
    }

  GetDepends(i);

  // now output the results
  for (i = 0; i < num; i++)
  {
    fprintf(fp," \\\n  %s",DependsStructArray[ dependIndices[i] ]->name);
  }
  fprintf(fp,"\n");
}




void AddToDepends(const char *filename)
{
  DEPENDS_STRUCT *dependsEntry;

// allocate new entry
  dependsEntry = new DEPENDS_STRUCT;

  dependsEntry->numIndices = 0;
  strcpy( dependsEntry->name, filename );

  if ( NumInDepends >= MAX_DEPENDS )
    {
    fprintf(stderr,
       "ERROR:  To many depends files... recompile with larger MAX_DEPENDS!!!");
    exit(-1);
    }

  DependsStructArray[ NumInDepends++ ] = dependsEntry;
}


int GetFullPath(char *name, const char *vtkHome, char *fullPath)
{
  struct stat statBuff;

  // does the file exist ?
  // search for it in the vtk src code
  sprintf(fullPath,"%s/common/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/graphics/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/imaging/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/contrib/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/patented/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/working/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/gemsvolume/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/gemsio/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/gemsip/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  // if control reaches here then it hasn't been found yet
  sprintf(fullPath,"%s/geae/%s",vtkHome,name);
  if (!stat(fullPath,&statBuff))
    return 1;

  return 0;
}



void GetIncludes(DEPENDS_STRUCT *dependsEntry, const char *vtkHome )
{
  ifstream *IS;
  char line[256];
  int j, k;
  char name[256], fullPath[512];
  struct stat statBuff;

  // does the file exist ?
  // search for it in the vtk src code
  if (stat(dependsEntry->name,&statBuff))
    {
    fprintf(stderr,"ERROR:  file %s not found... Continuing anyway!", dependsEntry->name);
    return;
    }

  IS = new ifstream(dependsEntry->name);

  // search for includes
  while (!IS->eof())
    {
    IS->getline(line,255);
    // do we have an include
    if (!strncmp(line,"#include",8))
      {
      // is it a quoted include
      for (j = 8; j < (int)strlen(line); j++)
        {
        if (line[j] == '<') j = 1000;
        else
          {
          if (line[j] == '"')
            {
            // we found a quoted include, process it
            // make sure it is a vtk include file
            if (!strncmp(line +j +1,"vtk",3))
              {
              // extract the class name
              // there should always be an end quote or this will die
              for (k = 0; line[j+k+1] != '"'; k++) 
                {
                name[k] = line[j+k+1];
                }
              name[k] = '\0';

              // Get the full name
              if (!GetFullPath(name, vtkHome, fullPath))
                {
                fprintf(stderr,"ERROR:  Dependency %s not found!!!", name);
                exit(-1);
                }

              // get the index in depends
              for (k = 0; k < NumInDepends; k++)
                {
                if ( !strcmp(fullPath,DependsStructArray[k]->name) )
                  {
                  dependsEntry->indices[ dependsEntry->numIndices++ ] = k;
                  break;
                  }
                }

              // if not found, add it to the end
              if ( k == NumInDepends )
                {
                AddToDepends(fullPath);
                dependsEntry->indices[ dependsEntry->numIndices++ ] = k;
                }

              break; // break for (j = 8... loop
              }
            else j = 1000;
            }
          } // end if line[j] == '<' and else
        } // end for j = 8 ... strlen(line)
      } // end if (!strncmp(line,"#include"))
    } // end while

  IS->close();
  delete IS;  

}



void BuildDepends(const char *vtkHome)
{
  int i;

  for (i = 0; i < NumInDepends; i++)
    GetIncludes(DependsStructArray[i],vtkHome);
}
