/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCaptionActor2D.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCaptionActor2D - draw text label associated with a point
// .SECTION Description
// vtkCaptionActor2D is a hybrid 2D/3D actor that is used to associate text 
// with a point (the AttachmentPoint) in the scene. The caption can be 
// drawn with a rectangular border and a leader connecting 
// the caption to the attachment point. Optionally, the leader can be 
// glyphed at its endpoint to create arrow heads or other indicators.
// 
// To use the caption actor, you normally specify the Position and Position2
// coordinates (these are inherited from the vtkActor2D superclass). (Note
// that Position2 can be set using vtkActor2D's SetWidth() and SetHeight()
// methods.)  Position and Position2 define the size of the caption, and a
// third point, the AttachmentPoint, defines a point that the caption is
// associated with.  You must also define the caption text, 
// whether you want a border around the caption, and whether you want a
// leader from the caption to the attachment point. The font attributes of 
// the text can be set through the vtkTextProperty associated to this actor.
// You also indicate whether you want
// the leader to be 2D or 3D. (2D leaders are always drawn over the
// underlying geometry. 3D leaders may be occluded by the geometry.) The
// leader may also be terminated by an optional glyph (e.g., arrow).
//
// The trickiest part about using this class is setting Position, Position2,
// and AttachmentPoint correctly. These instance variables are
// vtkCoordinates, and can be set up in various ways. In default usage, the
// AttachmentPoint is defined in the world coordinate system, Position is the
// lower-left corner of the caption and relative to AttachmentPoint (defined
// in display coordaintes, i.e., pixels), and Position2 is relative to
// Position and is the upper-right corner (also in display
// coordinates). However, the user has full control over the coordinates, and
// can do things like place the caption in a fixed position in the renderer,
// with the leader moving with the AttachmentPoint.

// .SECTION See Also
// vtkLegendBoxActor vtkTextMapper vtkTextActor vtkTextProperty
// vtkCoordinate

#ifndef __vtkCaptionActor2D_h
#define __vtkCaptionActor2D_h

#include "vtkActor2D.h"

class vtkActor;
class vtkAppendPolyData;
class vtkGlyph2D;
class vtkGlyph3D;
class vtkPolyData;
class vtkPolyDataMapper2D;
class vtkPolyDataMapper;
class vtkTextActor;
class vtkTextMapper;
class vtkTextProperty;

class VTK_HYBRID_EXPORT vtkCaptionActor2D : public vtkActor2D
{
public:
  vtkTypeRevisionMacro(vtkCaptionActor2D,vtkActor2D);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkCaptionActor2D *New();

  // Description:
  // Define the text to be placed in the caption. The text can be multiple
  // lines (separated by "\n").
  vtkSetStringMacro(Caption);
  vtkGetStringMacro(Caption);

  // Description:
  // Set/Get the attachment point for the caption. By default, the attachment
  // point is defined in world coordinates, but this can be changed using
  // vtkCoordinate methods.
  vtkWorldCoordinateMacro(AttachmentPoint);

  // Description:
  // Enable/disable the placement of a border around the text.
  vtkSetMacro(Border,int);
  vtkGetMacro(Border,int);
  vtkBooleanMacro(Border,int);

  // Description:
  // Enable/disable drawing a "line" from the caption to the 
  // attachment point.
  vtkSetMacro(Leader,int);
  vtkGetMacro(Leader,int);
  vtkBooleanMacro(Leader,int);

  // Description:
  // Indicate whether the leader is 2D (no hidden line) or 3D (z-buffered).
  vtkSetMacro(ThreeDimensionalLeader,int);
  vtkGetMacro(ThreeDimensionalLeader,int);
  vtkBooleanMacro(ThreeDimensionalLeader,int);

  // Description:
  // Specify a glyph to be used as the leader "head". This could be something
  // like an arrow or sphere. If not specified, no glyph is drawn.
  virtual void SetLeaderGlyph(vtkPolyData*);
  vtkGetObjectMacro(LeaderGlyph,vtkPolyData);

  // Description:
  // Specify the relative size of the leader head. This is expressed as a
  // fraction of the size (diagonal length) of the renderer. The leader
  // head is automatically scaled so that window resize, zooming or other 
  // camera motion results in proportional changes in size to the leader
  // glyph.
  vtkSetClampMacro(LeaderGlyphSize,float,0.0,0.1);
  vtkGetMacro(LeaderGlyphSize,float);

  // Description:
  // Specify the maximum size of the leader head (if any) in pixels. This 
  // is used in conjunction with LeaderGlyphSize to cap the maximum size of
  // the LeaderGlyph.
  vtkSetClampMacro(MaximumLeaderGlyphSize,int,1,1000);
  vtkGetMacro(MaximumLeaderGlyphSize,int);

  // Description:
  // Set/Get the padding between the caption and the border. The value
  // is specified in pixels.
  vtkSetClampMacro(Padding, int, 0, 50);
  vtkGetMacro(Padding, int);

  // Description:
  // Set/Get the text property.
  virtual void SetCaptionTextProperty(vtkTextProperty *p);
  vtkGetObjectMacro(CaptionTextProperty,vtkTextProperty);

  // Description:
  // Set/Get the font family. Three font types are allowed: Arial (VTK_ARIAL),
  // Courier (VTK_COURIER), and Times (VTK_TIMES).
  // Warning: these functions remain for backward compatibility. Use the
  // vtkTextProperty through the (Set/Get)CaptionTextProperty() methods.
  virtual void SetFontFamily(int val);
  virtual int GetFontFamily();
  void SetFontFamilyToArial()   { this->SetFontFamily(VTK_ARIAL);  };
  void SetFontFamilyToCourier() { this->SetFontFamily(VTK_COURIER);};
  void SetFontFamilyToTimes()   { this->SetFontFamily(VTK_TIMES);  };

  // Description:
  // Enable/disable text bolding.
  // Warning: these functions remain for backward compatibility. Use the
  // vtkTextProperty through the (Set/Get)CaptionTextProperty() methods.
  virtual void SetBold(int val);
  virtual int GetBold();
  vtkBooleanMacro(Bold, int);

  // Description:
  // Enable/disable text italic.
  // Warning: these functions remain for backward compatibility. Use the
  // vtkTextProperty through the (Set/Get)CaptionTextProperty() methods.
  virtual void SetItalic(int val);
  virtual int GetItalic();
  vtkBooleanMacro(Italic, int);

  // Description:
  // Enable/disable text shadows.
  // Warning: these functions remain for backward compatibility. Use the
  // vtkTextProperty through the (Set/Get)CaptionTextProperty() methods.
  virtual void SetShadow(int val);
  virtual int GetShadow();
  vtkBooleanMacro(Shadow, int);
    
  // Description:
  // Set/Get the horizontal justification to left (default), centered,
  // or right.
  // Warning: these functions remain for backward compatibility. Use the
  // vtkTextProperty through the (Set/Get)CaptionTextProperty() methods.
  virtual void SetJustification(int val);
  virtual int GetJustification();
  void SetJustificationToLeft()     
    { this->SetJustification(VTK_TEXT_LEFT);};
  void SetJustificationToCentered() 
    { this->SetJustification(VTK_TEXT_CENTERED);};
  void SetJustificationToRight()    
    { this->SetJustification(VTK_TEXT_RIGHT);};
    
  // Description:
  // Set/Get the vertical justification to bottom (default), middle,
  // or top.
  // Warning: these functions remain for backward compatibility. Use the
  // vtkTextProperty through the (Set/Get)CaptionTextProperty() methods.
  virtual void SetVerticalJustification(int val);
  virtual int GetVerticalJustification();
  void SetVerticalJustificationToBottom() 
    {this->SetVerticalJustification(VTK_TEXT_BOTTOM);};
  void SetVerticalJustificationToCentered() 
    {this->SetVerticalJustification(VTK_TEXT_CENTERED);};
  void SetVerticalJustificationToTop() 
    {this->SetVerticalJustification(VTK_TEXT_TOP);};
    
  // Description:
  // Shallow copy of this scaled text actor. Overloads the virtual
  // vtkProp method.
  void ShallowCopy(vtkProp *prop);

//BTX
  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS.
  // Release any graphics resources that are being consumed by this actor.
  // The parameter window could be used to determine which graphic
  // resources to release.
  virtual void ReleaseGraphicsResources(vtkWindow *);

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS.
  // Draw the legend box to the screen.
  int RenderOpaqueGeometry(vtkViewport* viewport);
  int RenderTranslucentGeometry(vtkViewport* ) {return 0;}
  int RenderOverlay(vtkViewport* viewport);
//ETX

protected:
  vtkCaptionActor2D();
  ~vtkCaptionActor2D();

  vtkCoordinate *AttachmentPointCoordinate;

  char  *Caption;
  int   Border;
  int   Leader;
  int   ThreeDimensionalLeader;
  float LeaderGlyphSize;
  int   MaximumLeaderGlyphSize;

  vtkPolyData *LeaderGlyph; //what to put on the end of the leader
  
  int   Padding;

private:
  vtkTextActor        *CaptionActor;
  vtkTextProperty     *CaptionTextProperty;

  vtkPolyData         *BorderPolyData;
  vtkPolyDataMapper2D *BorderMapper;
  vtkActor2D          *BorderActor;

  vtkPolyData         *HeadPolyData;    // single attachment point for glyphing
  vtkGlyph3D          *HeadGlyph;       // for 3D leader
  vtkPolyData         *LeaderPolyData;  // line represents the leader
  vtkAppendPolyData   *AppendLeader;    // append head and leader
  
  // for 2D leader
  vtkCoordinate       *MapperCoordinate2D;
  vtkPolyDataMapper2D *LeaderMapper2D;
  vtkActor2D          *LeaderActor2D;

  // for 3D leader
  vtkPolyDataMapper   *LeaderMapper3D;
  vtkActor            *LeaderActor3D;

private:
  vtkCaptionActor2D(const vtkCaptionActor2D&);  // Not implemented.
  void operator=(const vtkCaptionActor2D&);  // Not implemented.
};


#endif



