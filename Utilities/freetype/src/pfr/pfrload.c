/***************************************************************************/
/*                                                                         */
/*  pfrload.c                                                              */
/*                                                                         */
/*    FreeType PFR loader (body).                                          */
/*                                                                         */
/*  Copyright 2002 by                                                      */
/*  David Turner, Robert Wilhelm, and Werner Lemberg.                      */
/*                                                                         */
/*  This file is part of the FreeType project, and may only be used,       */
/*  modified, and distributed under the terms of the FreeType project      */
/*  license, LICENSE.TXT.  By continuing to use, modify, or distribute     */
/*  this file you indicate that you have read the license and              */
/*  understand and accept it fully.                                        */
/*                                                                         */
/***************************************************************************/


#include "pfrload.h"
#include FT_INTERNAL_DEBUG_H
#include FT_INTERNAL_STREAM_H

#include "pfrerror.h"

#undef  FT_COMPONENT
#define FT_COMPONENT  trace_pfr


  /*************************************************************************/
  /*************************************************************************/
  /*****                                                               *****/
  /*****                          EXTRA ITEMS                          *****/
  /*****                                                               *****/
  /*************************************************************************/
  /*************************************************************************/


  FT_LOCAL_DEF( FT_Error )
  pfr_extra_items_skip( FT_Byte*  *pp,
                        FT_Byte*   limit )
  {
    return pfr_extra_items_parse( pp, limit, NULL, NULL );
  }


  FT_LOCAL_DEF( FT_Error )
  pfr_extra_items_parse( FT_Byte*       *pp,
                         FT_Byte*        limit,
                         PFR_ExtraItem   item_list,
                         FT_Pointer      item_data )
  {
    FT_Error  error = 0;
    FT_Byte*  p     = *pp;
    FT_UInt   num_items, item_type, item_size;
    

    PFR_CHECK( 1 );
    num_items = PFR_NEXT_BYTE( p );

    for ( ; num_items > 0; num_items-- )
    {
      PFR_CHECK( 2 );
      item_size = PFR_NEXT_BYTE( p );
      item_type = PFR_NEXT_BYTE( p );
      
      PFR_CHECK( item_size );

      if ( item_list )
      {      
        PFR_ExtraItem  extra = item_list;


        for ( extra = item_list; extra->parser != NULL; extra++ )
        {
          if ( extra->type == item_type )
          {
            error = extra->parser( p, p + item_size, item_data );
            if ( error ) goto Exit;
            
            break;
          }
        }
      }
      
      p += item_size;
    }

  Exit:
    *pp = p;
    return error;
    
  Too_Short:
    FT_ERROR(( "pfr_extra_items_parse: invalid extra items table\n" ));
    error = PFR_Err_Invalid_Table;
    goto Exit;
  }


  /*************************************************************************/
  /*************************************************************************/
  /*****                                                               *****/
  /*****                          PFR HEADER                           *****/
  /*****                                                               *****/
  /*************************************************************************/
  /*************************************************************************/

   static const FT_Frame_Field  pfr_header_fields[] =
   {
#undef  FT_STRUCTURE
#define FT_STRUCTURE  PFR_HeaderRec

     FT_FRAME_START( 58 ),
       FT_FRAME_ULONG ( signature ),
       FT_FRAME_USHORT( version ),
       FT_FRAME_USHORT( signature2 ),
       FT_FRAME_USHORT( header_size ),

       FT_FRAME_USHORT( log_dir_size ),
       FT_FRAME_USHORT( log_dir_offset ),

       FT_FRAME_USHORT( log_font_max_size ),
       FT_FRAME_UOFF3 ( log_font_section_size ),
       FT_FRAME_UOFF3 ( log_font_section_offset ),

       FT_FRAME_USHORT( phy_font_max_size ),
       FT_FRAME_UOFF3 ( phy_font_section_size ),
       FT_FRAME_UOFF3 ( phy_font_section_offset ),

       FT_FRAME_USHORT( gps_max_size ),
       FT_FRAME_UOFF3 ( gps_section_size ),
       FT_FRAME_UOFF3 ( gps_section_offset ),

       FT_FRAME_BYTE  ( max_blue_values ),
       FT_FRAME_BYTE  ( max_x_orus ),
       FT_FRAME_BYTE  ( max_y_orus ),

       FT_FRAME_BYTE  ( phy_font_max_size_high ),
       FT_FRAME_BYTE  ( color_flags ),

       FT_FRAME_UOFF3 ( bct_max_size ),
       FT_FRAME_UOFF3 ( bct_set_max_size ),
       FT_FRAME_UOFF3 ( phy_bct_set_max_size ),

       FT_FRAME_USHORT( num_phy_fonts ),
       FT_FRAME_BYTE  ( max_vert_stem_snap ),
       FT_FRAME_BYTE  ( max_horz_stem_snap ),
       FT_FRAME_USHORT( max_chars ),
     FT_FRAME_END
   };


  FT_LOCAL_DEF( FT_Error )
  pfr_header_load( PFR_Header  header,
                   FT_Stream   stream )
  {
    FT_Error   error;


    /* read header directly */
    if ( !FT_STREAM_SEEK( 0 )                                &&
         !FT_STREAM_READ_FIELDS( pfr_header_fields, header ) )
    {
      /* make a few adjustments to the header */
      header->phy_font_max_size +=
        (FT_UInt32)header->phy_font_max_size_high << 16;
    }

    return error;
  }


  FT_LOCAL_DEF( FT_Bool )
  pfr_header_check( PFR_Header  header )
  {
    FT_Bool  result = 1;


    /* check signature and header size */
    if ( header->signature  != 0x50465230L ||   /* "PFR0" */
         header->version     > 4           ||
         header->header_size < 58          ||
         header->signature2 != 0x0d0a      )    /* CR/LF  */
    {
      result = 0;
    }
    return  result;
  }


  /***********************************************************************/
  /***********************************************************************/
  /*****                                                             *****/
  /*****                    PFR LOGICAL FONTS                        *****/
  /*****                                                             *****/
  /***********************************************************************/
  /***********************************************************************/


  FT_LOCAL_DEF( FT_Error )
  pfr_log_font_count( FT_Stream  stream,
                      FT_UInt32  section_offset,
                      FT_UInt   *acount )
  {
    FT_Error   error;
    FT_UInt    count;
    FT_UInt    result = 0;
    

    if ( FT_STREAM_SEEK( section_offset ) || FT_READ_USHORT( count ) )
      goto Exit;
    
    result = count;
    
  Exit:
    *acount = result;
    return error;
  }                      


  FT_LOCAL_DEF( FT_Error )
  pfr_log_font_load( PFR_LogFont  log_font,
                     FT_Stream    stream,
                     FT_UInt      index,
                     FT_UInt32    section_offset,
                     FT_Bool      size_increment )
  {
    FT_UInt    num_log_fonts;
    FT_UInt    flags;
    FT_UInt32  offset;
    FT_UInt32  size;
    FT_Error   error;


    if ( FT_STREAM_SEEK( section_offset ) ||
         FT_READ_USHORT( num_log_fonts )  )
      goto Exit;

    if ( index >= num_log_fonts )
      return PFR_Err_Invalid_Argument;

    if ( FT_STREAM_SKIP( index * 5 ) ||
         FT_READ_USHORT( size )      ||
         FT_READ_UOFF3 ( offset )    )
      goto Exit;

    /* save logical font size and offset */
    log_font->size   = size;
    log_font->offset = offset;

    /* now, check the rest of the table before loading it */
    {
      FT_Byte*  p;
      FT_Byte*  limit;
      FT_UInt   local;
      

      if ( FT_STREAM_SEEK( offset ) || FT_FRAME_ENTER( size ) )
        goto Exit;

      p     = stream->cursor;
      limit = p + size;
      
      PFR_CHECK(13);
      
      log_font->matrix[0] = PFR_NEXT_LONG( p );
      log_font->matrix[1] = PFR_NEXT_LONG( p );
      log_font->matrix[2] = PFR_NEXT_LONG( p );
      log_font->matrix[3] = PFR_NEXT_LONG( p );
      
      flags = PFR_NEXT_BYTE( p );
      
      local = 0;
      if ( flags & PFR_LOG_STROKE )
      {
        local++;
        if ( flags & PFR_LOG_2BYTE_STROKE )
          local++;
          
        if ( (flags & PFR_LINE_JOIN_MASK) == PFR_LINE_JOIN_MITER )
          local += 3;
      }
      if ( flags & PFR_LOG_BOLD )
      {
        local++;
        if ( flags & PFR_LOG_2BYTE_BOLD )
          local++;
      }

      PFR_CHECK( local );

      if ( flags & PFR_LOG_STROKE )
      {
        log_font->stroke_thickness = ( flags & PFR_LOG_2BYTE_STROKE )
                                     ? PFR_NEXT_SHORT( p )
                                     : PFR_NEXT_BYTE( p );
  
        if ( ( flags & PFR_LINE_JOIN_MASK ) == PFR_LINE_JOIN_MITER )
          log_font->miter_limit = PFR_NEXT_LONG( p );
      }
  
      if ( flags & PFR_LOG_BOLD )
      {
        log_font->bold_thickness = ( flags & PFR_LOG_2BYTE_BOLD )
                                   ? PFR_NEXT_SHORT( p )
                                   : PFR_NEXT_BYTE( p );
      }
  
      if ( flags & PFR_LOG_EXTRA_ITEMS )
      {
        error = pfr_extra_items_skip( &p, limit );
        if (error) goto Fail;
      }

      PFR_CHECK(5);
      log_font->phys_size   = PFR_NEXT_USHORT( p );
      log_font->phys_offset = PFR_NEXT_ULONG( p );
      if ( size_increment )
      {
        PFR_CHECK( 1 );
        log_font->phys_size += (FT_UInt32)PFR_NEXT_BYTE( p ) << 16;
      }
    }

  Fail:
    FT_FRAME_EXIT();

  Exit:
    return error;
  
  Too_Short:
    FT_ERROR(( "pfr_log_font_load: invalid logical font table\n" ));
    error = PFR_Err_Invalid_Table;
    goto Fail;
  }


  /***********************************************************************/
  /***********************************************************************/
  /*****                                                             *****/
  /*****                    PFR PHYSICAL FONTS                       *****/
  /*****                                                             *****/
  /***********************************************************************/
  /***********************************************************************/


  /* load bitmap strikes lists */
  FT_CALLBACK_DEF( FT_Error )
  pfr_extra_item_load_bitmap_info( FT_Byte*     p,
                                   FT_Byte*     limit,
                                   PFR_PhyFont  phy_font )
  {
    FT_Memory   memory = phy_font->memory;
    PFR_Strike  strike;
    FT_UInt     flags0;
    FT_UInt     n, count, size1;
    FT_Error    error = 0;

 
    PFR_CHECK( 5 );

    p += 3;  /* skip bctSize */
    flags0 = PFR_NEXT_BYTE( p );
    count  = PFR_NEXT_BYTE( p );

    /* re-allocate when needed */
    if ( phy_font->num_strikes + count > phy_font->max_strikes )
    {
      FT_UInt  new_max = (phy_font->num_strikes + count + 3) & -4;
      
      if ( FT_RENEW_ARRAY( phy_font->strikes,
                           phy_font->num_strikes,
                           new_max ) )
        goto Exit;
        
      phy_font->max_strikes = new_max;
    }

    size1 = 1 + 1 + 1 + 2 + 2 + 1;
    if ( flags0 & PFR_STRIKE_2BYTE_XPPM )
      size1++;

    if ( flags0 & PFR_STRIKE_2BYTE_YPPM )
      size1++;

    if ( flags0 & PFR_STRIKE_3BYTE_SIZE )
      size1++;

    if ( flags0 & PFR_STRIKE_3BYTE_OFFSET )
      size1++;

    if ( flags0 & PFR_STRIKE_2BYTE_COUNT )
      size1++;

    strike = phy_font->strikes + phy_font->num_strikes;

    PFR_CHECK( count * size1 );

    for ( n = 0; n < count; n++, strike++ )
    {
      strike->x_ppm       = ( flags0 & PFR_STRIKE_2BYTE_XPPM )
                            ? PFR_NEXT_USHORT( p )
                            : PFR_NEXT_BYTE( p );

      strike->y_ppm       = ( flags0 & PFR_STRIKE_2BYTE_YPPM )
                            ? PFR_NEXT_USHORT( p )
                            : PFR_NEXT_BYTE( p );

      strike->flags       = PFR_NEXT_BYTE( p );

      strike->bct_size    = ( flags0 & PFR_STRIKE_3BYTE_SIZE )
                            ? PFR_NEXT_ULONG( p )
                            : PFR_NEXT_USHORT( p );

      strike->bct_offset  = ( flags0 & PFR_STRIKE_3BYTE_OFFSET )
                            ? PFR_NEXT_ULONG( p )
                            : PFR_NEXT_USHORT( p );

      strike->num_bitmaps = ( flags0 & PFR_STRIKE_2BYTE_COUNT )
                            ? PFR_NEXT_USHORT( p )
                            : PFR_NEXT_BYTE( p );
    }

    phy_font->num_strikes += count;
    
  Exit:
    return error;
    
  Too_Short:
    error = PFR_Err_Invalid_Table;
    FT_ERROR(( "pfr_extra_item_load_bitmap_info: invalid bitmap info table\n" ));
    goto Exit;
  }                                  


  /* load font ID, i.e. name */
  FT_CALLBACK_DEF( FT_Error )
  pfr_extra_item_load_font_id( FT_Byte*     p,
                               FT_Byte*     limit,
                               PFR_PhyFont  phy_font )
  {
    FT_Error   error  = 0;
    FT_Memory  memory = phy_font->memory;
    FT_UInt    len    = (FT_UInt)( limit - p );
    

    if ( phy_font->font_id != NULL )
      goto Exit;
    
    if ( FT_ALLOC( phy_font->font_id, len+1 ) )
      goto Exit;

    /* copy font ID name, and terminate it for safety */
    FT_MEM_COPY( phy_font->font_id, p, len );
    phy_font->font_id[len] = 0;
    
  Exit:
    return error;
  }


  /* load stem snap tables */
  FT_CALLBACK_DEF( FT_Error )
  pfr_extra_item_load_stem_snaps( FT_Byte*     p,
                                  FT_Byte*     limit,
                                  PFR_PhyFont  phy_font )
  {
    FT_UInt    count, num_vert, num_horz;
    FT_Int*    snaps;
    FT_Error   error  = 0;
    FT_Memory  memory = phy_font->memory;


    if ( phy_font->vertical.stem_snaps != NULL )
      goto Exit;
    
    PFR_CHECK( 1 );
    count = PFR_NEXT_BYTE( p );
    
    num_vert = count & 15;
    num_horz = count >> 4;
    count    = num_vert + num_horz;
    
    PFR_CHECK( count * 2 );
    
    if ( FT_NEW_ARRAY( snaps, count ) )
      goto Exit;
    
    phy_font->vertical.stem_snaps = snaps;
    phy_font->horizontal.stem_snaps = snaps + num_vert;
    
    for ( ; count > 0; count--, snaps++ )
      *snaps = FT_NEXT_SHORT( p );
    
  Exit:
    return error;
  
  Too_Short:
    error = PFR_Err_Invalid_Table;
    FT_ERROR(( "pfr_exta_item_load_stem_snaps: invalid stem snaps table\n" ));
    goto Exit;
  }                                 


  /* load kerning pair data */
  FT_CALLBACK_DEF( FT_Error )
  pfr_extra_item_load_kerning_pairs( FT_Byte*     p,
                                     FT_Byte*     limit,
                                     PFR_PhyFont  phy_font )
  {
    FT_Int        count;
    FT_UShort     base_adj;
    FT_UInt       flags;
    FT_UInt       num_pairs;
    PFR_KernPair  pairs;
    FT_Error      error  = 0;
    FT_Memory     memory = phy_font->memory;


    /* XXX: there may be multiple extra items for kerning */
    if ( phy_font->kern_pairs != NULL )
      goto Exit;

    FT_TRACE2(( "pfr_extra_item_load_kerning_pairs()\n" ));
    
    PFR_CHECK( 4 );

    num_pairs = PFR_NEXT_BYTE( p );
    base_adj  = PFR_NEXT_SHORT( p );
    flags     = PFR_NEXT_BYTE( p );
    
#ifndef PFR_CONFIG_NO_CHECKS
    count = 3;
    
    if ( flags & PFR_KERN_2BYTE_CHAR )
      count += 2;
    
    if ( flags & PFR_KERN_2BYTE_ADJ )
      count += 1;

    PFR_CHECK( num_pairs * count );
#endif

    if ( FT_NEW_ARRAY( pairs, num_pairs ) )
      goto Exit;
    
    phy_font->num_kern_pairs = num_pairs;
    phy_font->kern_pairs     = pairs;
    
    for (count = num_pairs ; count > 0; count--, pairs++ )
    {
      if ( flags & PFR_KERN_2BYTE_CHAR )
      {
        pairs->glyph1 = PFR_NEXT_USHORT( p );
        pairs->glyph2 = PFR_NEXT_USHORT( p );
      }
      else
      {
        pairs->glyph1 = PFR_NEXT_BYTE( p );
        pairs->glyph2 = PFR_NEXT_BYTE( p );
      }
     
      if ( flags & PFR_KERN_2BYTE_ADJ )
        pairs->kerning.x = base_adj + PFR_NEXT_SHORT( p );
      else
        pairs->kerning.x = base_adj + PFR_NEXT_INT8( p );

      pairs->kerning.y = 0;

      FT_TRACE2(( "kerning %d <-> %d : %ld\n",
                   pairs->glyph1, pairs->glyph2, pairs->kerning.x ));
    }
      
  Exit:
    return error;
  
  Too_Short:
    error = PFR_Err_Invalid_Table;
    FT_ERROR(( "pfr_extra_item_load_kerning_pairs: "
               "invalid kerning pairs table\n" ));
    goto Exit;
  }                                 


  static const PFR_ExtraItemRec  pfr_phy_font_extra_items[] =
  {
    { 1, (PFR_ExtraItem_ParseFunc) pfr_extra_item_load_bitmap_info },
    { 2, (PFR_ExtraItem_ParseFunc) pfr_extra_item_load_font_id },
    { 3, (PFR_ExtraItem_ParseFunc) pfr_extra_item_load_stem_snaps },
    { 4, (PFR_ExtraItem_ParseFunc) pfr_extra_item_load_kerning_pairs },
    { 0, NULL }
  };


  FT_LOCAL_DEF( void )
  pfr_phy_font_done( PFR_PhyFont  phy_font,
                     FT_Memory    memory )
  {
    if ( phy_font->font_id )
      FT_FREE( phy_font->font_id );

    FT_FREE( phy_font->vertical.stem_snaps );
    phy_font->vertical.num_stem_snaps = 0;

    phy_font->horizontal.stem_snaps     = NULL;
    phy_font->horizontal.num_stem_snaps = 0;

    FT_FREE( phy_font->strikes );
    phy_font->num_strikes = 0;
    phy_font->max_strikes = 0;
  
    FT_FREE( phy_font->chars );
    phy_font->num_chars    = 0;
    phy_font->chars_offset = 0;

    FT_FREE( phy_font->blue_values );
    phy_font->num_blue_values = 0;

    FT_FREE( phy_font->kern_pairs );
    phy_font->num_kern_pairs = 0;
  }
  

  FT_LOCAL_DEF( FT_Error )
  pfr_phy_font_load( PFR_PhyFont  phy_font,
                     FT_Stream    stream,
                     FT_UInt32    offset,
                     FT_UInt32    size )
  {
    FT_Error   error;
    FT_Memory  memory = stream->memory;
    FT_UInt    flags, num_aux;
    FT_Byte*   p;
    FT_Byte*   limit;


    phy_font->memory = memory;
    phy_font->offset = offset;

    if ( FT_STREAM_SEEK( offset ) || FT_FRAME_ENTER( size ) )
      goto Exit;

    p     = stream->cursor;
    limit = p + size;

    PFR_CHECK( 15 );
    phy_font->font_ref_number    = PFR_NEXT_USHORT( p );
    phy_font->outline_resolution = PFR_NEXT_USHORT( p );
    phy_font->metrics_resolution = PFR_NEXT_USHORT( p );
    phy_font->bbox.xMin          = PFR_NEXT_SHORT( p );
    phy_font->bbox.yMin          = PFR_NEXT_SHORT( p );
    phy_font->bbox.xMax          = PFR_NEXT_SHORT( p );
    phy_font->bbox.yMax          = PFR_NEXT_SHORT( p );
    phy_font->flags      = flags = PFR_NEXT_BYTE( p );

    /* get the standard advance for non-proprotional fonts */
    if ( !(flags & PFR_PHY_PROPORTIONAL) )
    {
      PFR_CHECK( 2 );
      phy_font->standard_advance = PFR_NEXT_SHORT( p );
    }

    /* load the extra items when present */
    if ( flags & PFR_PHY_EXTRA_ITEMS )
    {
      error =  pfr_extra_items_parse( &p, limit,
                                      pfr_phy_font_extra_items, phy_font );
        
      if ( error )
        goto Fail;
    }
    
    /* skip the aux bytes */
    PFR_CHECK( 3 );
    num_aux = PFR_NEXT_ULONG( p );
    
    PFR_CHECK( num_aux );
    p += num_aux;

    /* read the blue values */
    {
      FT_UInt  n, count;

      PFR_CHECK( 1 );
      phy_font->num_blue_values = count = PFR_NEXT_BYTE( p );

      PFR_CHECK( count * 2 );
      
      if ( FT_NEW_ARRAY( phy_font->blue_values, count ) )
        goto Fail;

      for ( n = 0; n < count; n++ )
        phy_font->blue_values[n] = PFR_NEXT_SHORT( p );
    }

    PFR_CHECK( 8 );
    phy_font->blue_fuzz  = PFR_NEXT_BYTE( p );
    phy_font->blue_scale = PFR_NEXT_BYTE( p );

    phy_font->vertical.standard   = PFR_NEXT_USHORT( p );
    phy_font->horizontal.standard = PFR_NEXT_USHORT( p );

    /* read the character descriptors */
    {
      FT_UInt  n, count, Size;


      phy_font->num_chars    = count = PFR_NEXT_USHORT( p );
      phy_font->chars_offset = offset + ( p - stream->cursor );

      if ( FT_NEW_ARRAY( phy_font->chars, count ) )
        goto Fail;

      Size = 1 + 1 + 2;
      if ( flags & PFR_PHY_2BYTE_CHARCODE )
        Size += 1;
      
      if ( flags & PFR_PHY_PROPORTIONAL )
        Size += 2;
      
      if ( flags & PFR_PHY_ASCII_CODE )
        Size += 1;
        
      if ( flags & PFR_PHY_2BYTE_GPS_SIZE )
        Size += 1;
        
      if ( flags & PFR_PHY_3BYTE_GPS_OFFSET )
        Size += 1;
      
      PFR_CHECK( count * Size );

      for ( n = 0; n < count; n++ )
      {
        PFR_Char  cur = &phy_font->chars[n];


        cur->char_code = ( flags & PFR_PHY_2BYTE_CHARCODE )
                         ? PFR_NEXT_USHORT( p )
                         : PFR_NEXT_BYTE( p );

        cur->advance   = ( flags & PFR_PHY_PROPORTIONAL )
                         ? PFR_NEXT_SHORT( p )
                         : (FT_Int)phy_font->standard_advance;

#if 0
        cur->ascii     = ( flags & PFR_PHY_ASCII_CODE )
                         ? PFR_NEXT_BYTE( p )
                         : 0;
#else
        if ( flags & PFR_PHY_ASCII_CODE )
          p += 1;
#endif
        cur->gps_size  = ( flags & PFR_PHY_2BYTE_GPS_SIZE )
                         ? PFR_NEXT_USHORT( p )
                         : PFR_NEXT_BYTE( p );

        cur->gps_offset = ( flags & PFR_PHY_3BYTE_GPS_OFFSET )
                          ? PFR_NEXT_ULONG( p )
                          : PFR_NEXT_USHORT( p );
      }
    }

    /* that's it !! */
  Fail:
    FT_FRAME_EXIT();

  Exit:
    return error;
  
  Too_Short:
    error = PFR_Err_Invalid_Table;
    FT_ERROR(( "pfr_phy_font_load: invalid physical font table\n" ));
    goto Fail;
  }


/* END */
