/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVRMLExporter.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1998 Ken Martin, Will Schroeder, Bill Lorensen.

This software is copyrighted by Ken Martin, Will Schroeder and Bill Lorensen.
The following terms apply to all files associated with the software unless
explicitly disclaimed in individual files. This copyright specifically does
not apply to the related textbook "The Visualization Toolkit" ISBN
013199837-4 published by Prentice Hall which is covered by its own copyright.

The authors hereby grant permission to use, copy, and distribute this
software and its documentation for any purpose, provided that existing
copyright notices are retained in all copies and that this notice is included
verbatim in any distributions. Additionally, the authors grant permission to
modify this software and its documentation for any purpose, provided that
such modifications are not distributed without the explicit consent of the
authors and that existing copyright notices are retained in all copies. Some
of the algorithms implemented by this software are patented, observe all
applicable patent law.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN
"AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


=========================================================================*/
#include "vtkVRMLExporter.h"
#include "vtkGeometryFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkVRMLExporter* vtkVRMLExporter::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkVRMLExporter");
  if(ret)
    {
    return (vtkVRMLExporter*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkVRMLExporter;
}




vtkVRMLExporter::vtkVRMLExporter()
{
  this->Speed = 4.0;
  this->FileName = NULL;
  this->FilePointer = NULL;
}

vtkVRMLExporter::~vtkVRMLExporter()
{
  if ( this->FileName )
    {
    delete [] this->FileName;
    }
}

void vtkVRMLExporter::SetFilePointer(FILE *fp)
{
  if (fp != this->FilePointer)
    {
    this->Modified();
    this->FilePointer = fp;
    }
}

  
void vtkVRMLExporter::WriteData()
{
  vtkRenderer *ren;
  vtkActorCollection *ac;
  vtkActor *anActor, *aPart;
  vtkLightCollection *lc;
  vtkLight *aLight;
  vtkCamera *cam;
  float *tempf;
  FILE *fp;
  
  // make sure the user specified a FileName or FilePointer
  if (!this->FilePointer && (this->FileName == NULL))
    {
    vtkErrorMacro(<< "Please specify FileName to use");
    return;
    }

  // first make sure there is only one renderer in this rendering window
  if (this->RenderWindow->GetRenderers()->GetNumberOfItems() > 1)
    {
    vtkErrorMacro(<< "VRML files only support one renderer per window.");
    return;
    }

  // get the renderer
  this->RenderWindow->GetRenderers()->InitTraversal();
  ren = this->RenderWindow->GetRenderers()->GetNextItem();
  
  // make sure it has at least one actor
  if (ren->GetActors()->GetNumberOfItems() < 1)
    {
    vtkErrorMacro(<< "no actors found for writing VRML file.");
    return;
    }
    
  // try opening the files
  if (!this->FilePointer)
    {
    fp = fopen(this->FileName,"w");
    if (!fp)
      {
      vtkErrorMacro(<< "unable to open VRML file " << this->FileName);
      return;
      }
    }
  else
    {
    fp = this->FilePointer;
    }
  
  //
  //  Write header
  //
  vtkDebugMacro("Writing VRML file");
  fprintf(fp,"#VRML V2.0 utf8\n");
  fprintf(fp,"# VRML file written by the visualization toolkit\n\n");

  // do the camera
  cam = ren->GetActiveCamera();
  fprintf(fp,"    Viewpoint\n      {\n      fieldOfView %f\n",
	  cam->GetViewAngle()*3.1415926/180.0);
  fprintf(fp,"      position %f %f %f\n",cam->GetPosition()[0],
	  cam->GetPosition()[1], cam->GetPosition()[2]);
  fprintf(fp,"      description \"Default View\"\n");
  tempf = cam->GetOrientationWXYZ();
  fprintf(fp,"      orientation %g %g %g %g\n      }\n", tempf[1], tempf[2], 
	  tempf[3], tempf[0]*3.1415926/180.0);

  // do the lights first the ambient then the others
  fprintf(fp,"    NavigationInfo {\n      type [\"EXAMINE\",\"FLY\"]\n      speed %f\n", this->Speed);
  if (ren->GetLights()->GetNumberOfItems() == 0)
    {
    fprintf(fp,"      headlight TRUE}\n\n");
    }
  else
    {
    fprintf(fp,"      headlight FALSE}\n\n");
    }
  fprintf(fp,"    DirectionalLight { ambientIntensity 1 intensity 0 # ambient light\n");
  fprintf(fp,"      color %f %f %f }\n\n", ren->GetAmbient()[0],
	  ren->GetAmbient()[1], ren->GetAmbient()[2]);
  
  // make sure we have a default light
  // if we dont then use a headlight
  lc = ren->GetLights();
  for (lc->InitTraversal(); (aLight = lc->GetNextItem()); )
    {
    this->WriteALight(aLight, fp);
    }

  // do the actors now
  ac = ren->GetActors();
  for (ac->InitTraversal(); (anActor = ac->GetNextActor()); )
    {
    for (anActor->InitPartTraversal();(aPart=anActor->GetNextPart()); )
      {
      this->WriteAnActor(aPart, fp);
      }
    }

  if (!this->FilePointer)
    {
    fclose(fp);
    }
}

void vtkVRMLExporter::WriteALight(vtkLight *aLight, FILE *fp)
{
  float *pos, *focus, *color;
  float dir[3];
  
  pos = aLight->GetPosition();
  focus = aLight->GetFocalPoint();
  color = aLight->GetColor();

  dir[0] = focus[0] - pos[0];
  dir[1] = focus[1] - pos[1];
  dir[2] = focus[2] - pos[2];
  vtkMath::Normalize(dir);
    
  if (aLight->GetPositional())
    {
    float *attn;
    
    if (aLight->GetConeAngle() >= 180.0)
      {
      fprintf(fp,"    PointLight {\n");
      }
    else
      { 
      fprintf(fp,"    SpotLight {\n");
      fprintf(fp,"      direction %f %f %f\n",dir[0], dir[1], dir[2]);
      fprintf(fp,"      cutOffAngle %f\n", aLight->GetConeAngle());
      }
    fprintf(fp,"      location %f %f %f\n", pos[0], pos[1], pos[2]);
    attn = aLight->GetAttenuationValues();
    fprintf(fp,"      attenuation %f %f %f\n", attn[0], attn[1], attn[2]);
    }
  else
    {
    fprintf(fp,"    DirectionalLight {\n");
    fprintf(fp,"      direction %f %f %f\n",dir[0], dir[1], dir[2]);
    }

  fprintf(fp,"      color %f %f %f\n", color[0], color[1], color[2]);
  fprintf(fp,"      intensity %f\n", aLight->GetIntensity());
  if (aLight->GetSwitch())
    {
    fprintf(fp,"      on TRUE\n      }\n");
    }
  else
    {
    fprintf(fp,"      on FALSE\n      }\n");
    }
}

void vtkVRMLExporter::WriteAnActor(vtkActor *anActor, FILE *fp)
{
  vtkDataSet *ds;
  vtkPolyData *pd;
  vtkGeometryFilter *gf = NULL;
  vtkPointData *pntData;
  vtkPoints *points = NULL;
  vtkNormals *normals = NULL;
  vtkTCoords *tcoords = NULL;
  int i, i1, i2;
  vtkProperty *prop;
  float *tempf;
  vtkCellArray *cells;
  int npts, *indx;
  float tempf2;
  int pointDataWritten = 0;
  vtkPolyDataMapper *pm;
  vtkScalars *colors;
  float *p;
  unsigned char *c;
  vtkTransform *trans;
  int totalValues;
  
  // see if the actor has a mapper. it could be an assembly
  if (anActor->GetMapper() == NULL)
    {
    return;
    }

  // first stuff out the transform
  trans = vtkTransform::New();
  trans->SetMatrix(*(anActor->vtkProp3D::GetMatrixPointer()));
  
  fprintf(fp,"    Transform {\n");
  tempf = trans->GetPosition();
  fprintf(fp,"      translation %g %g %g\n", tempf[0], tempf[1], tempf[2]);
  tempf = trans->GetOrientationWXYZ();
  fprintf(fp,"      rotation %g %g %g %g\n", tempf[1], tempf[2], 
	  tempf[3], tempf[0]*3.1415926/180.0);
  tempf = trans->GetScale();
  fprintf(fp,"      scale %g %g %g\n", tempf[0], tempf[1], tempf[2]);
  fprintf(fp,"      children [\n");
  trans->Delete();
  
  // get the mappers input and matrix
  ds = anActor->GetMapper()->GetInput();
  
  // we really want polydata
  if ( ds->GetDataObjectType() != VTK_POLY_DATA )
    {
    gf = vtkGeometryFilter::New();
    gf->SetInput(ds);
    gf->Update();
    pd = gf->GetOutput();
    }
  else
    {
    pd = (vtkPolyData *)ds;
    }

  pm = vtkPolyDataMapper::New();
  pm->SetInput(pd);
  pm->SetScalarRange(anActor->GetMapper()->GetScalarRange());
  pm->SetScalarVisibility(anActor->GetMapper()->GetScalarVisibility());
  pm->SetLookupTable(anActor->GetMapper()->GetLookupTable());
  pm->SetScalarMode(anActor->GetMapper()->GetScalarMode());

  points = pd->GetPoints();
  pntData = pd->GetPointData();
  normals = pntData->GetNormals();
  tcoords = pntData->GetTCoords();
  colors  = pm->GetColors();
  
  fprintf(fp,"        Shape {\n");
  
  // write out the material properties to the mat file
  fprintf(fp,"          appearance Appearance {\n");
  fprintf(fp,"            material Material {\n");
  prop = anActor->GetProperty();
  fprintf(fp,"              ambientIntensity %g\n", prop->GetAmbient());
  // if we don't have colors and we have only lines & points
  // use emissive to color them
  if (!(normals || colors || pd->GetNumberOfPolys() || 
	pd->GetNumberOfStrips()))
    {
    tempf2 = prop->GetAmbient();
    tempf = prop->GetAmbientColor();
    fprintf(fp,"              emissiveColor %g %g %g\n",
	    tempf[0]*tempf2, tempf[1]*tempf2, tempf[2]*tempf2);
    }
  tempf2 = prop->GetDiffuse();
  tempf = prop->GetDiffuseColor();
  fprintf(fp,"              diffuseColor %g %g %g\n",
	  tempf[0]*tempf2, tempf[1]*tempf2, tempf[2]*tempf2);
  tempf2 = prop->GetSpecular();
  tempf = prop->GetSpecularColor();
  fprintf(fp,"              specularColor %g %g %g\n",
	  tempf[0]*tempf2, tempf[1]*tempf2, tempf[2]*tempf2);
  fprintf(fp,"              shininess %g\n",prop->GetSpecularPower()/128.0);
  fprintf(fp,"              transparency %g\n",1.0 - prop->GetOpacity());
  fprintf(fp,"              }\n"); // close matrial

  // is there a texture map
  if (anActor->GetTexture())
    {
    vtkTexture *aTexture = anActor->GetTexture();
    int *size, xsize, ysize, bpp;
    vtkScalars *scalars;
    vtkScalars *mappedScalars;
    unsigned char *txtrData;
    
    // make sure it is updated and then get some info
    if (aTexture->GetInput() == NULL)
      {
      vtkErrorMacro(<< "texture has no input!\n");
      return;
      }
    aTexture->GetInput()->Update();
    size = aTexture->GetInput()->GetDimensions();
    scalars = (aTexture->GetInput()->GetPointData())->GetScalars();

    // make sure scalars are non null
    if (!scalars) 
      {
      vtkErrorMacro(<< "No scalar values found for texture input!\n");
      return;
      }

    // make sure using unsigned char data of color scalars type
    if (aTexture->GetMapColorScalarsThroughLookupTable () ||
        (scalars->GetDataType() != VTK_UNSIGNED_CHAR) )
      {
      mappedScalars = aTexture->GetMappedScalars ();
      }
    else
      {
      mappedScalars = scalars;
      }

    // we only support 2d texture maps right now
    // so one of the three sizes must be 1, but it 
    // could be any of them, so lets find it
    if (size[0] == 1)
      {
      xsize = size[1]; ysize = size[2];
      }
    else
      {
      xsize = size[0];
      if (size[1] == 1)
	{
	ysize = size[2];
	}
      else
	{
	ysize = size[1];
	if (size[2] != 1)
	  {
	  vtkErrorMacro(<< "3D texture maps currently are not supported!\n");
	  return;
	  }
	}
      }

    fprintf(fp,"            texture PixelTexture {\n");
    bpp = mappedScalars->GetNumberOfComponents();
    fprintf(fp,"              image %i %i %i\n", xsize, ysize, bpp);
    txtrData = ((vtkUnsignedCharArray *)mappedScalars->GetData())->GetPointer(0);
    totalValues = xsize*ysize;
    for (i = 0; i < totalValues; i++)
      {
      fprintf(fp,"0x%.2x",*txtrData);
      txtrData++;
      if (bpp > 1) 
	{
	fprintf(fp,"%.2x",*txtrData);
	txtrData++;
	}
      if (bpp > 2) 
	{
	fprintf(fp,"%.2x",*txtrData);
	txtrData++;
	}
      if (bpp > 3) 
	{
	fprintf(fp,"%.2x",*txtrData);
	txtrData++;
	}
      if (i%8 == 0)
	{
	fprintf(fp,"\n");
	}
      else
	{
	fprintf(fp," ");
	}
      }
    if (!(aTexture->GetRepeat()))
      {
      fprintf(fp,"              repeatS FALSE\n");
      fprintf(fp,"              repeatT FALSE\n");
      }
    fprintf(fp,"              }\n"); // close texture
    }
  fprintf(fp,"            }\n"); // close appearance

  // write out polys if any
  if (pd->GetNumberOfPolys() > 0)
    {
    fprintf(fp,"          geometry IndexedFaceSet {\n");
    // two sided lighting ? for now assume it is on
    fprintf(fp,"            solid FALSE\n");
    if (!pointDataWritten)
      {
      this->WritePointData(points, normals, tcoords, colors, fp);
      pointDataWritten = 1;
      }
    else
      {
      fprintf(fp,"            coord  USE VTKcoordinates\n");
      if (normals)
	{
	fprintf(fp,"            normal  USE VTKnormals\n");
	}
      if (tcoords)
	{
	fprintf(fp,"            texCoord  USE VTKtcoords\n");
	}
      if (colors)
	{
	fprintf(fp,"            color  USE VTKcolors\n");
	}
      }
    
    fprintf(fp,"            coordIndex  [\n");
    
    cells = pd->GetPolys();
    for (cells->InitTraversal(); cells->GetNextCell(npts,indx); )
      {
      fprintf(fp,"              ");
      for (i = 0; i < npts; i++)
	{
	fprintf(fp,"%i, ",indx[i]);
	}
      fprintf(fp,"-1,\n");
      }
    fprintf(fp,"            ]\n");
    fprintf(fp,"          }\n");
    }

  // write out tstrips if any
  if (pd->GetNumberOfStrips() > 0)
    {
    fprintf(fp,"          geometry IndexedFaceSet {\n");
    if (!pointDataWritten)
      {
      this->WritePointData(points, normals, tcoords, colors, fp);
      pointDataWritten = 1;
      }
    else
      {
      fprintf(fp,"            coord  USE VTKcoordinates\n");
      if (normals)
	{
	fprintf(fp,"            normal  USE VTKnormals\n");
	}
      if (tcoords)
	{
	fprintf(fp,"            texCoord  USE VTKtcoords\n");
	}
      if (colors)
	{
	fprintf(fp,"            color  USE VTKcolors\n");
	}
      }
    fprintf(fp,"            coordIndex  [\n");
    cells = pd->GetStrips();
    for (cells->InitTraversal(); cells->GetNextCell(npts,indx); )
      {
      for (i = 2; i < npts; i++)
	{
	if (i%2)
	  {
	  i1 = i - 1;
	  i2 = i - 2;
	  }
	else
	  {
	  i1 = i - 2;
	  i2 = i - 1;
	  }
	fprintf(fp,"              %i, %i, %i, -1,\n",indx[i1], 
		indx[i2], indx[i]);
	}
      }
    fprintf(fp,"            ]\n");
    fprintf(fp,"          }\n");
    }
  
  // write out lines if any
  if (pd->GetNumberOfLines() > 0)
    {
    fprintf(fp,"          geometry IndexedLineSet {\n");
    if (!pointDataWritten)
      {
      this->WritePointData(points, NULL, NULL, colors, fp);
      pointDataWritten = 1;
      }
    else
      {
      fprintf(fp,"            coord  USE VTKcoordinates\n");
      if (colors)
	{
	fprintf(fp,"            color  USE VTKcolors\n");
	}
      }
    
    fprintf(fp,"            coordIndex  [\n");
    
    cells = pd->GetLines();
    for (cells->InitTraversal(); cells->GetNextCell(npts,indx); )
      {
      fprintf(fp,"              ");
      for (i = 0; i < npts; i++)
	{
	fprintf(fp,"%i, ",indx[i]);
	}
      fprintf(fp,"-1,\n");
      }
    fprintf(fp,"            ]\n");
    fprintf(fp,"          }\n");
    }

  // write out verts if any
  if (pd->GetNumberOfVerts() > 0)
    {
    fprintf(fp,"          geometry PointSet {\n");
    cells = pd->GetVerts();
    fprintf(fp,"            coord Coordinate {");
    fprintf(fp,"              point [");
    for (cells->InitTraversal(); cells->GetNextCell(npts,indx); )
      {
      fprintf(fp,"              ");
      for (i = 0; i < npts; i++)
	{
	p = points->GetPoint(indx[i]);
	fprintf (fp,"              %g %g %g,\n", p[0], p[1], p[2]);
	}
      }
    fprintf(fp,"              ]\n");
    fprintf(fp,"            }\n");
    if (colors)
	{
	fprintf(fp,"            color Color {");
	fprintf(fp,"              color [");
	for (cells->InitTraversal(); cells->GetNextCell(npts,indx); )
	  {
	  fprintf(fp,"              ");
	  for (i = 0; i < npts; i++)
	    {
	    c = colors->GetColor(indx[i]);
	    fprintf (fp,"           %g %g %g,\n", c[0]/255.0, c[1]/255.0, 
		     c[2]/255.0);
	    }
	  }
	fprintf(fp,"              ]\n");
	fprintf(fp,"            }\n");
	}
  
    fprintf(fp,"          }\n");
    }

  fprintf(fp,"        }\n"); // close the  Shape
  fprintf(fp,"      ]\n"); // close the original transforms children
  fprintf(fp,"    }\n"); // close the original transform
  
  if (gf)
    {
    gf->Delete();
    }
  pm->Delete();
}

void vtkVRMLExporter::WritePointData(vtkPoints *points, vtkNormals *normals,
				     vtkTCoords *tcoords, 
				     vtkScalars *colors, FILE *fp)
{
  float *p;
  int i;
  unsigned char *c;
  
  // write out the points
  fprintf(fp,"            coord DEF VTKcoordinates Coordinate {\n");
  fprintf(fp,"              point [\n");
  for (i = 0; i < points->GetNumberOfPoints(); i++)
    {
    p = points->GetPoint(i);
    fprintf (fp,"              %g %g %g,\n", p[0], p[1], p[2]);
    }
  fprintf(fp,"              ]\n");
  fprintf(fp,"            }\n");
  
  // write out the point data
  if (normals)
    {
    fprintf(fp,"            normal DEF VTKnormals Normal {\n");
    fprintf(fp,"              vector [\n");
    for (i = 0; i < normals->GetNumberOfNormals(); i++)
      {
      p = normals->GetNormal(i);
      fprintf (fp,"           %g %g %g,\n", p[0], p[1], p[2]);
      }
    fprintf(fp,"            ]\n");
    fprintf(fp,"          }\n");
    }

  // write out the point data
  if (tcoords)
    {
    fprintf(fp,"            texCoord DEF VTKtcoords TextureCoordinate {\n");
    fprintf(fp,"              point [\n");
    for (i = 0; i < tcoords->GetNumberOfTCoords(); i++)
      {
      p = tcoords->GetTCoord(i);
      fprintf (fp,"           %g %g,\n", p[0], p[1]);
      }
    fprintf(fp,"            ]\n");
    fprintf(fp,"          }\n");
    }

  // write out the point data
  if (colors)
    {
    fprintf(fp,"            color DEF VTKcolors Color {\n");
    fprintf(fp,"              color [\n");
    for (i = 0; i < colors->GetNumberOfScalars(); i++)
      {
      c = colors->GetColor(i);
      fprintf (fp,"           %g %g %g,\n", c[0]/255.0, c[1]/255.0, 
	       c[2]/255.0);
      }
    fprintf(fp,"            ]\n");
    fprintf(fp,"          }\n");
    }
}


void vtkVRMLExporter::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkExporter::PrintSelf(os,indent);
 
  os << indent << "File Name: " << this->FileName << "\n";
  os << indent << "Speed: " << this->Speed << "\n";
}

