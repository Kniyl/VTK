#include  "FTPixmapGlyph.h"
#include  "FTGLgl.h"
#ifdef FTGL_DEBUG
  #include "mmgr.h"
#endif

#ifndef GetCurrentColorFunctionName
#define GetCurrentColorFunctionName GetCurrentColorOpenGL
#endif

void FTPixmapGlyph::GetCurrentColorFunctionName(float colour[4],
                                                const FTGLRenderContext *context)
{
  glGetFloatv( GL_CURRENT_COLOR, colour);
}

#ifndef RenderFunctionName
#define RenderFunctionName RenderOpenGL
#endif

#define ToString(arg) ToString0(arg)
#define ToString0(arg) #arg

void FTPixmapGlyph::RenderFunctionName(const FT_Vector& pen,
                                       const FTGLRenderContext *context)
{
  // Move the glyph origin
  glBitmap( 0, 0, 0.0, 0.0, (float)(pen.x + pos.x), (float)(pen.y - pos.y), (const GLubyte *)0);
  
  printf("FTPixmapGlyph::"ToString(RenderFunctionName)"\n");

  glDrawPixels( destWidth, destHeight, GL_RGBA, GL_UNSIGNED_BYTE, (const GLvoid*)data);

  // Restore the glyph origin
  glBitmap( 0, 0, 0.0, 0.0, (float)(-pen.x - pos.x), (float)(-pen.y + pos.y), (const GLubyte *)0);
}
