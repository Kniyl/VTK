/*

VTK_FREETYPE_CHANGE this file is new for VTK.

This header file mangles all symbols exported from the freetype library.
It is included in all files while building the freetype library.  Due to
namespace pollution, no freetype headers should be included in .h files in
VTK.

The following command was used to obtain the symbol list:

nm libvtkfreetype.a |grep " [TR] "

This is the way to recreate the whole list:

nm bin/libvtkfreetype.so |grep " [TR] " | awk '{ print "#define "$3" vtk_freetype_"$3 }'

*/

#ifndef vtk_freetype_mangle_h
#define vtk_freetype_mangle_h

#define _FT_Activate_Size vtk_freetype__FT_Activate_Size
#define _FT_Add_Module vtk_freetype__FT_Add_Module
#define _FT_Alloc vtk_freetype__FT_Alloc
#define _FT_Angle_Diff vtk_freetype__FT_Angle_Diff
#define _FT_Atan2 vtk_freetype__FT_Atan2
#define _FT_Attach_File vtk_freetype__FT_Attach_File
#define _FT_Attach_Stream vtk_freetype__FT_Attach_Stream
#define _FT_CMap_Done vtk_freetype__FT_CMap_Done
#define _FT_CMap_New vtk_freetype__FT_CMap_New
#define _FT_CeilFix vtk_freetype__FT_CeilFix
#define _FT_Cos vtk_freetype__FT_Cos
#define _FT_DivFix vtk_freetype__FT_DivFix
#define _FT_Done_Face vtk_freetype__FT_Done_Face
#define _FT_Done_GlyphSlot vtk_freetype__FT_Done_GlyphSlot
#define _FT_Done_Library vtk_freetype__FT_Done_Library
#define _FT_Done_Size vtk_freetype__FT_Done_Size
#define _FT_FloorFix vtk_freetype__FT_FloorFix
#define _FT_Free vtk_freetype__FT_Free
#define _FT_GetFilePath_From_Mac_ATS_Name vtk_freetype__FT_GetFilePath_From_Mac_ATS_Name
#define _FT_GetFile_From_Mac_ATS_Name vtk_freetype__FT_GetFile_From_Mac_ATS_Name
#define _FT_GetFile_From_Mac_Name vtk_freetype__FT_GetFile_From_Mac_Name
#define _FT_Get_CMap_Language_ID vtk_freetype__FT_Get_CMap_Language_ID
#define _FT_Get_Char_Index vtk_freetype__FT_Get_Char_Index
#define _FT_Get_Charmap_Index vtk_freetype__FT_Get_Charmap_Index
#define _FT_Get_First_Char vtk_freetype__FT_Get_First_Char
#define _FT_Get_Glyph_Name vtk_freetype__FT_Get_Glyph_Name
#define _FT_Get_Kerning vtk_freetype__FT_Get_Kerning
#define _FT_Get_Module vtk_freetype__FT_Get_Module
#define _FT_Get_Module_Interface vtk_freetype__FT_Get_Module_Interface
#define _FT_Get_Name_Index vtk_freetype__FT_Get_Name_Index
#define _FT_Get_Next_Char vtk_freetype__FT_Get_Next_Char
#define _FT_Get_Postscript_Name vtk_freetype__FT_Get_Postscript_Name
#define _FT_Get_Renderer vtk_freetype__FT_Get_Renderer
#define _FT_Get_Sfnt_Name vtk_freetype__FT_Get_Sfnt_Name
#define _FT_Get_Sfnt_Name_Count vtk_freetype__FT_Get_Sfnt_Name_Count
#define _FT_Get_Sfnt_Table vtk_freetype__FT_Get_Sfnt_Table
#define _FT_Get_SubGlyph_Info vtk_freetype__FT_Get_SubGlyph_Info
#define _FT_Get_Track_Kerning vtk_freetype__FT_Get_Track_Kerning
#define _FT_Get_TrueType_Engine_Type vtk_freetype__FT_Get_TrueType_Engine_Type
#define _FT_GlyphLoader_Add vtk_freetype__FT_GlyphLoader_Add
#define _FT_GlyphLoader_CheckPoints vtk_freetype__FT_GlyphLoader_CheckPoints
#define _FT_GlyphLoader_CheckSubGlyphs vtk_freetype__FT_GlyphLoader_CheckSubGlyphs
#define _FT_GlyphLoader_CopyPoints vtk_freetype__FT_GlyphLoader_CopyPoints
#define _FT_GlyphLoader_CreateExtra vtk_freetype__FT_GlyphLoader_CreateExtra
#define _FT_GlyphLoader_Done vtk_freetype__FT_GlyphLoader_Done
#define _FT_GlyphLoader_New vtk_freetype__FT_GlyphLoader_New
#define _FT_GlyphLoader_Prepare vtk_freetype__FT_GlyphLoader_Prepare
#define _FT_GlyphLoader_Reset vtk_freetype__FT_GlyphLoader_Reset
#define _FT_GlyphLoader_Rewind vtk_freetype__FT_GlyphLoader_Rewind
#define _FT_Library_Version vtk_freetype__FT_Library_Version
#define _FT_List_Add vtk_freetype__FT_List_Add
#define _FT_List_Finalize vtk_freetype__FT_List_Finalize
#define _FT_List_Find vtk_freetype__FT_List_Find
#define _FT_List_Insert vtk_freetype__FT_List_Insert
#define _FT_List_Iterate vtk_freetype__FT_List_Iterate
#define _FT_List_Remove vtk_freetype__FT_List_Remove
#define _FT_List_Up vtk_freetype__FT_List_Up
#define _FT_Load_Char vtk_freetype__FT_Load_Char
#define _FT_Load_Glyph vtk_freetype__FT_Load_Glyph
#define _FT_Load_Sfnt_Table vtk_freetype__FT_Load_Sfnt_Table
#define _FT_Lookup_Renderer vtk_freetype__FT_Lookup_Renderer
#define _FT_Match_Size vtk_freetype__FT_Match_Size
#define _FT_MulDiv vtk_freetype__FT_MulDiv
#define _FT_MulDiv_No_Round vtk_freetype__FT_MulDiv_No_Round
#define _FT_MulFix vtk_freetype__FT_MulFix
#define _FT_New_Face vtk_freetype__FT_New_Face
#define _FT_New_Face_From_FOND vtk_freetype__FT_New_Face_From_FOND
#define _FT_New_Face_From_FSRef vtk_freetype__FT_New_Face_From_FSRef
#define _FT_New_Face_From_FSSpec vtk_freetype__FT_New_Face_From_FSSpec
#define _FT_New_GlyphSlot vtk_freetype__FT_New_GlyphSlot
#define _FT_New_Library vtk_freetype__FT_New_Library
#define _FT_New_Memory_Face vtk_freetype__FT_New_Memory_Face
#define _FT_New_Size vtk_freetype__FT_New_Size
#define _FT_Open_Face vtk_freetype__FT_Open_Face
#define _FT_Outline_Check vtk_freetype__FT_Outline_Check
#define _FT_Outline_Copy vtk_freetype__FT_Outline_Copy
#define _FT_Outline_Decompose vtk_freetype__FT_Outline_Decompose
#define _FT_Outline_Done vtk_freetype__FT_Outline_Done
#define _FT_Outline_Done_Internal vtk_freetype__FT_Outline_Done_Internal
#define _FT_Outline_Embolden vtk_freetype__FT_Outline_Embolden
#define _FT_Outline_Get_Bitmap vtk_freetype__FT_Outline_Get_Bitmap
#define _FT_Outline_Get_CBox vtk_freetype__FT_Outline_Get_CBox
#define _FT_Outline_Get_Orientation vtk_freetype__FT_Outline_Get_Orientation
#define _FT_Outline_New vtk_freetype__FT_Outline_New
#define _FT_Outline_New_Internal vtk_freetype__FT_Outline_New_Internal
#define _FT_Outline_Render vtk_freetype__FT_Outline_Render
#define _FT_Outline_Reverse vtk_freetype__FT_Outline_Reverse
#define _FT_Outline_Transform vtk_freetype__FT_Outline_Transform
#define _FT_Outline_Translate vtk_freetype__FT_Outline_Translate
#define _FT_QAlloc vtk_freetype__FT_QAlloc
#define _FT_QRealloc vtk_freetype__FT_QRealloc
#define _FT_Raccess_Get_DataOffsets vtk_freetype__FT_Raccess_Get_DataOffsets
#define _FT_Raccess_Get_HeaderInfo vtk_freetype__FT_Raccess_Get_HeaderInfo
#define _FT_Raccess_Guess vtk_freetype__FT_Raccess_Guess
#define _FT_Realloc vtk_freetype__FT_Realloc
#define _FT_Remove_Module vtk_freetype__FT_Remove_Module
#define _FT_Render_Glyph vtk_freetype__FT_Render_Glyph
#define _FT_Render_Glyph_Internal vtk_freetype__FT_Render_Glyph_Internal
#define _FT_Request_Metrics vtk_freetype__FT_Request_Metrics
#define _FT_Request_Size vtk_freetype__FT_Request_Size
#define _FT_RoundFix vtk_freetype__FT_RoundFix
#define _FT_Select_Charmap vtk_freetype__FT_Select_Charmap
#define _FT_Select_Metrics vtk_freetype__FT_Select_Metrics
#define _FT_Select_Size vtk_freetype__FT_Select_Size
#define _FT_Set_Char_Size vtk_freetype__FT_Set_Char_Size
#define _FT_Set_Charmap vtk_freetype__FT_Set_Charmap
#define _FT_Set_Debug_Hook vtk_freetype__FT_Set_Debug_Hook
#define _FT_Set_Pixel_Sizes vtk_freetype__FT_Set_Pixel_Sizes
#define _FT_Set_Renderer vtk_freetype__FT_Set_Renderer
#define _FT_Set_Transform vtk_freetype__FT_Set_Transform
#define _FT_Sfnt_Table_Info vtk_freetype__FT_Sfnt_Table_Info
#define _FT_Sin vtk_freetype__FT_Sin
#define _FT_Sqrt32 vtk_freetype__FT_Sqrt32
#define _FT_SqrtFixed vtk_freetype__FT_SqrtFixed
#define _FT_Stream_Close vtk_freetype__FT_Stream_Close
#define _FT_Stream_EnterFrame vtk_freetype__FT_Stream_EnterFrame
#define _FT_Stream_ExitFrame vtk_freetype__FT_Stream_ExitFrame
#define _FT_Stream_ExtractFrame vtk_freetype__FT_Stream_ExtractFrame
#define _FT_Stream_Free vtk_freetype__FT_Stream_Free
#define _FT_Stream_GetChar vtk_freetype__FT_Stream_GetChar
#define _FT_Stream_GetLong vtk_freetype__FT_Stream_GetLong
#define _FT_Stream_GetLongLE vtk_freetype__FT_Stream_GetLongLE
#define _FT_Stream_GetOffset vtk_freetype__FT_Stream_GetOffset
#define _FT_Stream_GetShort vtk_freetype__FT_Stream_GetShort
#define _FT_Stream_GetShortLE vtk_freetype__FT_Stream_GetShortLE
#define _FT_Stream_New vtk_freetype__FT_Stream_New
#define _FT_Stream_OpenMemory vtk_freetype__FT_Stream_OpenMemory
#define _FT_Stream_Pos vtk_freetype__FT_Stream_Pos
#define _FT_Stream_Read vtk_freetype__FT_Stream_Read
#define _FT_Stream_ReadAt vtk_freetype__FT_Stream_ReadAt
#define _FT_Stream_ReadChar vtk_freetype__FT_Stream_ReadChar
#define _FT_Stream_ReadFields vtk_freetype__FT_Stream_ReadFields
#define _FT_Stream_ReadLong vtk_freetype__FT_Stream_ReadLong
#define _FT_Stream_ReadLongLE vtk_freetype__FT_Stream_ReadLongLE
#define _FT_Stream_ReadOffset vtk_freetype__FT_Stream_ReadOffset
#define _FT_Stream_ReadShort vtk_freetype__FT_Stream_ReadShort
#define _FT_Stream_ReadShortLE vtk_freetype__FT_Stream_ReadShortLE
#define _FT_Stream_ReleaseFrame vtk_freetype__FT_Stream_ReleaseFrame
#define _FT_Stream_Seek vtk_freetype__FT_Stream_Seek
#define _FT_Stream_Skip vtk_freetype__FT_Stream_Skip
#define _FT_Stream_TryRead vtk_freetype__FT_Stream_TryRead
#define _FT_Tan vtk_freetype__FT_Tan
#define _FT_Vector_From_Polar vtk_freetype__FT_Vector_From_Polar
#define _FT_Vector_Length vtk_freetype__FT_Vector_Length
#define _FT_Vector_Polarize vtk_freetype__FT_Vector_Polarize
#define _FT_Vector_Rotate vtk_freetype__FT_Vector_Rotate
#define _FT_Vector_Transform vtk_freetype__FT_Vector_Transform
#define _FT_Vector_Unit vtk_freetype__FT_Vector_Unit
#define _ft_corner_is_flat vtk_freetype__ft_corner_is_flat
#define _ft_corner_orientation vtk_freetype__ft_corner_orientation
#define _ft_glyphslot_alloc_bitmap vtk_freetype__ft_glyphslot_alloc_bitmap
#define _ft_glyphslot_free_bitmap vtk_freetype__ft_glyphslot_free_bitmap
#define _ft_glyphslot_set_bitmap vtk_freetype__ft_glyphslot_set_bitmap
#define _ft_highpow2 vtk_freetype__ft_highpow2
#define _ft_mem_alloc vtk_freetype__ft_mem_alloc
#define _ft_mem_dup vtk_freetype__ft_mem_dup
#define _ft_mem_free vtk_freetype__ft_mem_free
#define _ft_mem_qalloc vtk_freetype__ft_mem_qalloc
#define _ft_mem_qrealloc vtk_freetype__ft_mem_qrealloc
#define _ft_mem_realloc vtk_freetype__ft_mem_realloc
#define _ft_mem_strcpyn vtk_freetype__ft_mem_strcpyn
#define _ft_mem_strdup vtk_freetype__ft_mem_strdup
#define _ft_module_get_service vtk_freetype__ft_module_get_service
#define _ft_service_list_lookup vtk_freetype__ft_service_list_lookup
#define _ft_stub_set_char_sizes vtk_freetype__ft_stub_set_char_sizes
#define _ft_stub_set_pixel_sizes vtk_freetype__ft_stub_set_pixel_sizes
#define _ft_synthesize_vertical_metrics vtk_freetype__ft_synthesize_vertical_metrics
#define _ft_validator_error vtk_freetype__ft_validator_error
#define _ft_validator_init vtk_freetype__ft_validator_init
#define _ft_validator_run vtk_freetype__ft_validator_run
#define _FT_Outline_Get_BBox vtk_freetype__FT_Outline_Get_BBox
#define _FT_Bitmap_Convert vtk_freetype__FT_Bitmap_Convert
#define _FT_Bitmap_Copy vtk_freetype__FT_Bitmap_Copy
#define _FT_Bitmap_Done vtk_freetype__FT_Bitmap_Done
#define _FT_Bitmap_Embolden vtk_freetype__FT_Bitmap_Embolden
#define _FT_Bitmap_New vtk_freetype__FT_Bitmap_New
#define _FT_Done_Glyph vtk_freetype__FT_Done_Glyph
#define _FT_Get_Glyph vtk_freetype__FT_Get_Glyph
#define _FT_Glyph_Copy vtk_freetype__FT_Glyph_Copy
#define _FT_Glyph_Get_CBox vtk_freetype__FT_Glyph_Get_CBox
#define _FT_Glyph_To_Bitmap vtk_freetype__FT_Glyph_To_Bitmap
#define _FT_Glyph_Transform vtk_freetype__FT_Glyph_Transform
#define _FT_Matrix_Invert vtk_freetype__FT_Matrix_Invert
#define _FT_Matrix_Multiply vtk_freetype__FT_Matrix_Multiply
#define _FT_Add_Default_Modules vtk_freetype__FT_Add_Default_Modules
#define _FT_Done_FreeType vtk_freetype__FT_Done_FreeType
#define _FT_Init_FreeType vtk_freetype__FT_Init_FreeType
#define _FT_Get_MM_Var vtk_freetype__FT_Get_MM_Var
#define _FT_Get_Multi_Master vtk_freetype__FT_Get_Multi_Master
#define _FT_Set_MM_Blend_Coordinates vtk_freetype__FT_Set_MM_Blend_Coordinates
#define _FT_Set_MM_Design_Coordinates vtk_freetype__FT_Set_MM_Design_Coordinates
#define _FT_Set_Var_Blend_Coordinates vtk_freetype__FT_Set_Var_Blend_Coordinates
#define _FT_Set_Var_Design_Coordinates vtk_freetype__FT_Set_Var_Design_Coordinates
#define _FTC_CMapCache_Lookup vtk_freetype__FTC_CMapCache_Lookup
#define _FTC_CMapCache_New vtk_freetype__FTC_CMapCache_New
#define _FTC_ImageCache_Lookup vtk_freetype__FTC_ImageCache_Lookup
#define _FTC_ImageCache_New vtk_freetype__FTC_ImageCache_New
#define _FTC_Image_Cache_Lookup vtk_freetype__FTC_Image_Cache_Lookup
#define _FTC_Image_Cache_New vtk_freetype__FTC_Image_Cache_New
#define _FTC_Manager_Done vtk_freetype__FTC_Manager_Done
#define _FTC_Manager_LookupFace vtk_freetype__FTC_Manager_LookupFace
#define _FTC_Manager_LookupSize vtk_freetype__FTC_Manager_LookupSize
#define _FTC_Manager_Lookup_Face vtk_freetype__FTC_Manager_Lookup_Face
#define _FTC_Manager_Lookup_Size vtk_freetype__FTC_Manager_Lookup_Size
#define _FTC_Manager_New vtk_freetype__FTC_Manager_New
#define _FTC_Manager_RemoveFaceID vtk_freetype__FTC_Manager_RemoveFaceID
#define _FTC_Manager_Reset vtk_freetype__FTC_Manager_Reset
#define _FTC_Node_Unref vtk_freetype__FTC_Node_Unref
#define _FTC_SBitCache_Lookup vtk_freetype__FTC_SBitCache_Lookup
#define _FTC_SBitCache_New vtk_freetype__FTC_SBitCache_New
#define _FTC_SBit_Cache_Lookup vtk_freetype__FTC_SBit_Cache_Lookup
#define _FTC_SBit_Cache_New vtk_freetype__FTC_SBit_Cache_New
#define _ftc_node_destroy vtk_freetype__ftc_node_destroy
#define _FT_Stream_OpenGzip vtk_freetype__FT_Stream_OpenGzip
#define _FT_Stream_OpenLZW vtk_freetype__FT_Stream_OpenLZW
#define _ps_hints_apply vtk_freetype__ps_hints_apply
#define _TT_New_Context vtk_freetype__TT_New_Context
#define _TT_RunIns vtk_freetype__TT_RunIns
#define _FT_Trace_Get_Count vtk_freetype__FT_Trace_Get_Count
#define _FT_Trace_Get_Name vtk_freetype__FT_Trace_Get_Name
#define _ft_debug_init vtk_freetype__ft_debug_init
#define _FT_Done_Memory vtk_freetype__FT_Done_Memory
#define _FT_New_Memory vtk_freetype__FT_New_Memory
#define _FT_Stream_Open vtk_freetype__FT_Stream_Open

#endif
