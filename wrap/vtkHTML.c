#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vtkParse.h"

/* print the parsed structures */
void vtkParseOutput(FILE *fp, FileInfo *data)
{
  int i,closed;
  
  fprintf(fp,"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n");
  fprintf(fp,"<!-- manual pages generated by vtkHTML -->\n\n");
  fprintf(fp,"<HTML><HEADER>\n<TITLE>Documentation for %s</TITLE></HEADER>\n",
	  data->ClassName);
  
  fprintf(fp,"<body bgcolor=\"#FFFFFF\">\n");
  
  fprintf(fp,"<p>%s</p>\n\n",  
	  data->NameComment ? data->NameComment : "No description provided.");
  
  if (data->NumberOfSuperClasses > 0)
    {
    fprintf(fp,"<p>Super Class: <A HREF=\"%s.html\">%s</a></p>\n", 
	    data->SuperClasses[0], data->SuperClasses[0]);
    }
  
  fprintf(fp,"<p>Description:</p>\n\n<blockquote>\n");
  fprintf(fp,"    <p>%s</p>\n</blockquote>\n\n",  
	  data->Description ? data->Description : "None provided.");
  
  fprintf(fp,"<p>&nbsp;</p>\n\n");
  
  if (data->Caveats)
    {
    fprintf(fp,"<p>Caveats:</p>\n\n<blockquote>\n");
    fprintf(fp,"    <p>%s</p>\n</blockquote>\n\n", data->Caveats);
    fprintf(fp,"<p>&nbsp;</p>\n\n");
    }
  
  if (data->SeeAlso)
    {
    char *dup, *tok;
    
    fprintf(fp,"<p>See Also:</p>\n\n<blockquote>\n");
    fprintf(fp,"    <p>");
    /* change simple see also into html references */
    dup = strdup(data->SeeAlso);
    tok = strtok(dup," ");
    while (tok)
      {
      if (!strncmp("vtk",tok,3))
        {
        fprintf(fp,"<a href=\"%s.html\">%s</a> ",tok, tok);
        }
      else
        {
        fprintf(fp,"%s ",tok);
        }
      tok = strtok(NULL," ");
      }
    fprintf(fp,"</p>\n</blockquote>\n<p>&nbsp;</p>\n\n");
    free(dup);
    }
  
  fprintf(fp,"<p>Methods:</p>\n\n<ul>\n\n");
  for (i = 0; i < data->NumberOfFunctions; i++)
    {
    if (data->Functions[i].Name)
      {
      if (data->Functions[i].Comment)
        {
        fprintf(fp,"    <li><a href=\"#method%i\">%s</a></li>", i, 
		data->Functions[i].Signature);
	}
      else
        {
        fprintf(fp,"    <li>%s</li>", data->Functions[i].Signature);
        }
      }
    }
  fprintf(fp,"</ul>\n\n");
  fprintf(fp,"<p>&nbsp;</p>\n\n");
  
  fprintf(fp,"<p>Detailed Method Descriptions:</p><hr>\n");
  closed = 1;
  for (i = 0; i < data->NumberOfFunctions; i++)
    {
    if (data->Functions[i].Name && data->Functions[i].Comment)
      {
      if (closed) 
        {
        fprintf(fp,"<ul>\n");
        closed = 0;
        }
      fprintf(fp,"    <li><a name=\"method%i\">%s</a></li>\n", i, 
	      data->Functions[i].Signature);
      /* if we have a comment */
      if (data->Functions[i].Comment)
	{
	/* and this is the last method */
	if (i == (data->NumberOfFunctions - 1) ||
	    /* or it isn't the last */
	    (i < (data->NumberOfFunctions - 1) &&
	     /* and the next comment either doesn't exist */
	     (!data->Functions[i+1].Comment ||
	      /* or it does exist but is different */
	      (data->Functions[i+1].Comment &&
	       strcmp(data->Functions[i].Comment,
		      data->Functions[i+1].Comment)))))
	  {
	  closed = 1;
	  fprintf(fp,"<blockquote>");
	  fprintf(fp,"   <p>%s</p>\n</blockquote></ul><hr>\n\n",
		  data->Functions[i].Comment);
	  }
	}
      }
    }
  fprintf(fp,"<p>&nbsp;</p>\n</body>\n</html>\n");
}

