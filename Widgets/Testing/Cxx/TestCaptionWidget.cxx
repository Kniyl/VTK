/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestCaptionWidget.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This example tests the vtkCaptionWidget.

// First include the required header files for the VTK classes we are using.
#include "vtkSmartPointer.h"
#include "vtkCaptionWidget.h"
#include "vtkCaptionRepresentation.h"
#include "vtkCaptionActor2D.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkTextActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCommand.h"
#include "vtkInteractorEventRecorder.h"
#include "vtkTextProperty.h"

const char eventLog[] =
"# StreamVersion 1\n"
"EnterEvent 56 3 0 0 0 0 0\n"
"MouseMoveEvent 46 5 0 0 0 0 0\n"
"MouseMoveEvent 45 33 0 0 0 0 0\n"
"MiddleButtonPressEvent 45 33 0 0 0 0 0\n"
"MouseMoveEvent 46 33 0 0 0 0 0\n"
"RenderEvent 46 33 0 0 0 0 0\n"
"MouseMoveEvent 47 33 0 0 0 0 0\n"
"RenderEvent 47 33 0 0 0 0 0\n"
"MouseMoveEvent 48 33 0 0 0 0 0\n"
"RenderEvent 48 33 0 0 0 0 0\n"
"MouseMoveEvent 49 33 0 0 0 0 0\n"
"RenderEvent 49 33 0 0 0 0 0\n"
"MouseMoveEvent 52 33 0 0 0 0 0\n"
"RenderEvent 52 33 0 0 0 0 0\n"
"MouseMoveEvent 53 33 0 0 0 0 0\n"
"RenderEvent 53 33 0 0 0 0 0\n"
"MouseMoveEvent 55 33 0 0 0 0 0\n"
"RenderEvent 55 33 0 0 0 0 0\n"
"MouseMoveEvent 56 33 0 0 0 0 0\n"
"RenderEvent 56 33 0 0 0 0 0\n"
"MouseMoveEvent 58 33 0 0 0 0 0\n"
"RenderEvent 58 33 0 0 0 0 0\n"
"MouseMoveEvent 59 33 0 0 0 0 0\n"
"RenderEvent 59 33 0 0 0 0 0\n"
"MouseMoveEvent 62 33 0 0 0 0 0\n"
"RenderEvent 62 33 0 0 0 0 0\n"
"MouseMoveEvent 67 33 0 0 0 0 0\n"
"RenderEvent 67 33 0 0 0 0 0\n"
"MouseMoveEvent 71 33 0 0 0 0 0\n"
"RenderEvent 71 33 0 0 0 0 0\n"
"MouseMoveEvent 79 33 0 0 0 0 0\n"
"RenderEvent 79 33 0 0 0 0 0\n"
"MouseMoveEvent 89 33 0 0 0 0 0\n"
"RenderEvent 89 33 0 0 0 0 0\n"
"MouseMoveEvent 100 33 0 0 0 0 0\n"
"RenderEvent 100 33 0 0 0 0 0\n"
"MouseMoveEvent 116 33 0 0 0 0 0\n"
"RenderEvent 116 33 0 0 0 0 0\n"
"MouseMoveEvent 127 33 0 0 0 0 0\n"
"RenderEvent 127 33 0 0 0 0 0\n"
"MouseMoveEvent 138 33 0 0 0 0 0\n"
"RenderEvent 138 33 0 0 0 0 0\n"
"MouseMoveEvent 142 33 0 0 0 0 0\n"
"RenderEvent 142 33 0 0 0 0 0\n"
"MouseMoveEvent 145 33 0 0 0 0 0\n"
"RenderEvent 145 33 0 0 0 0 0\n"
"MouseMoveEvent 146 33 0 0 0 0 0\n"
"RenderEvent 146 33 0 0 0 0 0\n"
"MouseMoveEvent 150 33 0 0 0 0 0\n"
"RenderEvent 150 33 0 0 0 0 0\n"
"MouseMoveEvent 154 33 0 0 0 0 0\n"
"RenderEvent 154 33 0 0 0 0 0\n"
"MouseMoveEvent 156 33 0 0 0 0 0\n"
"RenderEvent 156 33 0 0 0 0 0\n"
"MouseMoveEvent 161 33 0 0 0 0 0\n"
"RenderEvent 161 33 0 0 0 0 0\n"
"MouseMoveEvent 163 33 0 0 0 0 0\n"
"RenderEvent 163 33 0 0 0 0 0\n"
"MouseMoveEvent 167 33 0 0 0 0 0\n"
"RenderEvent 167 33 0 0 0 0 0\n"
"MouseMoveEvent 169 33 0 0 0 0 0\n"
"RenderEvent 169 33 0 0 0 0 0\n"
"MouseMoveEvent 175 33 0 0 0 0 0\n"
"RenderEvent 175 33 0 0 0 0 0\n"
"MouseMoveEvent 180 33 0 0 0 0 0\n"
"RenderEvent 180 33 0 0 0 0 0\n"
"MouseMoveEvent 183 33 0 0 0 0 0\n"
"RenderEvent 183 33 0 0 0 0 0\n"
"MouseMoveEvent 193 33 0 0 0 0 0\n"
"RenderEvent 193 33 0 0 0 0 0\n"
"MouseMoveEvent 196 33 0 0 0 0 0\n"
"RenderEvent 196 33 0 0 0 0 0\n"
"MouseMoveEvent 200 33 0 0 0 0 0\n"
"RenderEvent 200 33 0 0 0 0 0\n"
"MouseMoveEvent 202 33 0 0 0 0 0\n"
"RenderEvent 202 33 0 0 0 0 0\n"
"MouseMoveEvent 206 33 0 0 0 0 0\n"
"RenderEvent 206 33 0 0 0 0 0\n"
"MouseMoveEvent 207 33 0 0 0 0 0\n"
"RenderEvent 207 33 0 0 0 0 0\n"
"MouseMoveEvent 211 33 0 0 0 0 0\n"
"RenderEvent 211 33 0 0 0 0 0\n"
"MouseMoveEvent 212 33 0 0 0 0 0\n"
"RenderEvent 212 33 0 0 0 0 0\n"
"MouseMoveEvent 214 33 0 0 0 0 0\n"
"RenderEvent 214 33 0 0 0 0 0\n"
"MouseMoveEvent 215 33 0 0 0 0 0\n"
"RenderEvent 215 33 0 0 0 0 0\n"
"MouseMoveEvent 216 33 0 0 0 0 0\n"
"RenderEvent 216 33 0 0 0 0 0\n"
"MouseMoveEvent 218 32 0 0 0 0 0\n"
"RenderEvent 218 32 0 0 0 0 0\n"
"MouseMoveEvent 220 32 0 0 0 0 0\n"
"RenderEvent 220 32 0 0 0 0 0\n"
"MouseMoveEvent 224 31 0 0 0 0 0\n"
"RenderEvent 224 31 0 0 0 0 0\n"
"MouseMoveEvent 228 31 0 0 0 0 0\n"
"RenderEvent 228 31 0 0 0 0 0\n"
"MouseMoveEvent 229 31 0 0 0 0 0\n"
"RenderEvent 229 31 0 0 0 0 0\n"
"MouseMoveEvent 232 31 0 0 0 0 0\n"
"RenderEvent 232 31 0 0 0 0 0\n"
"MouseMoveEvent 233 31 0 0 0 0 0\n"
"RenderEvent 233 31 0 0 0 0 0\n"
"MouseMoveEvent 236 31 0 0 0 0 0\n"
"RenderEvent 236 31 0 0 0 0 0\n"
"MouseMoveEvent 237 31 0 0 0 0 0\n"
"RenderEvent 237 31 0 0 0 0 0\n"
"MouseMoveEvent 238 31 0 0 0 0 0\n"
"RenderEvent 238 31 0 0 0 0 0\n"
"MouseMoveEvent 239 31 0 0 0 0 0\n"
"RenderEvent 239 31 0 0 0 0 0\n"
"MouseMoveEvent 240 31 0 0 0 0 0\n"
"RenderEvent 240 31 0 0 0 0 0\n"
"MouseMoveEvent 241 31 0 0 0 0 0\n"
"RenderEvent 241 31 0 0 0 0 0\n"
"MouseMoveEvent 242 31 0 0 0 0 0\n"
"RenderEvent 242 31 0 0 0 0 0\n"
"MouseMoveEvent 243 31 0 0 0 0 0\n"
"RenderEvent 243 31 0 0 0 0 0\n"
"MouseMoveEvent 244 31 0 0 0 0 0\n"
"RenderEvent 244 31 0 0 0 0 0\n"
"MouseMoveEvent 245 31 0 0 0 0 0\n"
"RenderEvent 245 31 0 0 0 0 0\n"
"MouseMoveEvent 246 31 0 0 0 0 0\n"
"RenderEvent 246 31 0 0 0 0 0\n"
"MouseMoveEvent 248 30 0 0 0 0 0\n"
"RenderEvent 248 30 0 0 0 0 0\n"
"MouseMoveEvent 249 30 0 0 0 0 0\n"
"RenderEvent 249 30 0 0 0 0 0\n"
"MouseMoveEvent 253 30 0 0 0 0 0\n"
"RenderEvent 253 30 0 0 0 0 0\n"
"MouseMoveEvent 254 30 0 0 0 0 0\n"
"RenderEvent 254 30 0 0 0 0 0\n"
"MiddleButtonReleaseEvent 254 30 0 0 0 0 0\n"
"MouseMoveEvent 253 31 0 0 0 0 0\n"
"MouseMoveEvent 160 43 0 0 0 0 0\n"
"MiddleButtonPressEvent 160 43 0 0 0 0 0\n"
"MouseMoveEvent 160 44 0 0 0 0 0\n"
"RenderEvent 160 44 0 0 0 0 0\n"
"MouseMoveEvent 160 45 0 0 0 0 0\n"
"RenderEvent 160 45 0 0 0 0 0\n"
"MouseMoveEvent 160 46 0 0 0 0 0\n"
"RenderEvent 160 46 0 0 0 0 0\n"
"MouseMoveEvent 160 47 0 0 0 0 0\n"
"RenderEvent 160 47 0 0 0 0 0\n"
"MouseMoveEvent 160 48 0 0 0 0 0\n"
"RenderEvent 160 48 0 0 0 0 0\n"
"MouseMoveEvent 160 49 0 0 0 0 0\n"
"RenderEvent 160 49 0 0 0 0 0\n"
"MouseMoveEvent 160 50 0 0 0 0 0\n"
"RenderEvent 160 50 0 0 0 0 0\n"
"MouseMoveEvent 160 51 0 0 0 0 0\n"
"RenderEvent 160 51 0 0 0 0 0\n"
"MouseMoveEvent 160 52 0 0 0 0 0\n"
"RenderEvent 160 52 0 0 0 0 0\n"
"MouseMoveEvent 160 53 0 0 0 0 0\n"
"RenderEvent 160 53 0 0 0 0 0\n"
"MouseMoveEvent 160 54 0 0 0 0 0\n"
"RenderEvent 160 54 0 0 0 0 0\n"
"MouseMoveEvent 160 55 0 0 0 0 0\n"
"RenderEvent 160 55 0 0 0 0 0\n"
"MouseMoveEvent 160 56 0 0 0 0 0\n"
"RenderEvent 160 56 0 0 0 0 0\n"
"MouseMoveEvent 160 57 0 0 0 0 0\n"
"RenderEvent 160 57 0 0 0 0 0\n"
"MiddleButtonReleaseEvent 160 57 0 0 0 0 0\n"
"MouseMoveEvent 159 57 0 0 0 0 0\n"
"MouseMoveEvent 148 36 0 0 0 0 0\n"
"MouseMoveEvent 147 36 0 0 0 0 0\n"
"LeftButtonPressEvent 147 36 0 0 0 0 0\n"
"MouseMoveEvent 147 37 0 0 0 0 0\n"
"RenderEvent 147 37 0 0 0 0 0\n"
"MouseMoveEvent 148 37 0 0 0 0 0\n"
"RenderEvent 148 37 0 0 0 0 0\n"
"MouseMoveEvent 148 38 0 0 0 0 0\n"
"RenderEvent 148 38 0 0 0 0 0\n"
"MouseMoveEvent 149 38 0 0 0 0 0\n"
"RenderEvent 149 38 0 0 0 0 0\n"
"MouseMoveEvent 150 38 0 0 0 0 0\n"
"RenderEvent 150 38 0 0 0 0 0\n"
"MouseMoveEvent 150 39 0 0 0 0 0\n"
"RenderEvent 150 39 0 0 0 0 0\n"
"MouseMoveEvent 151 40 0 0 0 0 0\n"
"RenderEvent 151 40 0 0 0 0 0\n"
"MouseMoveEvent 151 41 0 0 0 0 0\n"
"RenderEvent 151 41 0 0 0 0 0\n"
"MouseMoveEvent 152 43 0 0 0 0 0\n"
"RenderEvent 152 43 0 0 0 0 0\n"
"MouseMoveEvent 153 46 0 0 0 0 0\n"
"RenderEvent 153 46 0 0 0 0 0\n"
"MouseMoveEvent 154 47 0 0 0 0 0\n"
"RenderEvent 154 47 0 0 0 0 0\n"
"MouseMoveEvent 156 50 0 0 0 0 0\n"
"RenderEvent 156 50 0 0 0 0 0\n"
"MouseMoveEvent 156 53 0 0 0 0 0\n"
"RenderEvent 156 53 0 0 0 0 0\n"
"MouseMoveEvent 158 56 0 0 0 0 0\n"
"RenderEvent 158 56 0 0 0 0 0\n"
"MouseMoveEvent 159 60 0 0 0 0 0\n"
"RenderEvent 159 60 0 0 0 0 0\n"
"MouseMoveEvent 160 61 0 0 0 0 0\n"
"RenderEvent 160 61 0 0 0 0 0\n"
"MouseMoveEvent 161 65 0 0 0 0 0\n"
"RenderEvent 161 65 0 0 0 0 0\n"
"MouseMoveEvent 163 69 0 0 0 0 0\n"
"RenderEvent 163 69 0 0 0 0 0\n"
"MouseMoveEvent 163 71 0 0 0 0 0\n"
"RenderEvent 163 71 0 0 0 0 0\n"
"MouseMoveEvent 164 74 0 0 0 0 0\n"
"RenderEvent 164 74 0 0 0 0 0\n"
"MouseMoveEvent 165 76 0 0 0 0 0\n"
"RenderEvent 165 76 0 0 0 0 0\n"
"MouseMoveEvent 165 79 0 0 0 0 0\n"
"RenderEvent 165 79 0 0 0 0 0\n"
"MouseMoveEvent 166 83 0 0 0 0 0\n"
"RenderEvent 166 83 0 0 0 0 0\n"
"MouseMoveEvent 166 85 0 0 0 0 0\n"
"RenderEvent 166 85 0 0 0 0 0\n"
"MouseMoveEvent 167 89 0 0 0 0 0\n"
"RenderEvent 167 89 0 0 0 0 0\n"
"MouseMoveEvent 167 93 0 0 0 0 0\n"
"RenderEvent 167 93 0 0 0 0 0\n"
"MouseMoveEvent 167 96 0 0 0 0 0\n"
"RenderEvent 167 96 0 0 0 0 0\n"
"MouseMoveEvent 168 100 0 0 0 0 0\n"
"RenderEvent 168 100 0 0 0 0 0\n"
"MouseMoveEvent 168 102 0 0 0 0 0\n"
"RenderEvent 168 102 0 0 0 0 0\n"
"MouseMoveEvent 169 106 0 0 0 0 0\n"
"RenderEvent 169 106 0 0 0 0 0\n"
"MouseMoveEvent 169 109 0 0 0 0 0\n"
"RenderEvent 169 109 0 0 0 0 0\n"
"MouseMoveEvent 169 111 0 0 0 0 0\n"
"RenderEvent 169 111 0 0 0 0 0\n"
"MouseMoveEvent 169 114 0 0 0 0 0\n"
"RenderEvent 169 114 0 0 0 0 0\n"
"MouseMoveEvent 169 116 0 0 0 0 0\n"
"RenderEvent 169 116 0 0 0 0 0\n"
"MouseMoveEvent 170 119 0 0 0 0 0\n"
"RenderEvent 170 119 0 0 0 0 0\n"
"MouseMoveEvent 172 128 0 0 0 0 0\n"
"RenderEvent 172 128 0 0 0 0 0\n"
"MouseMoveEvent 172 130 0 0 0 0 0\n"
"RenderEvent 172 130 0 0 0 0 0\n"
"MouseMoveEvent 174 139 0 0 0 0 0\n"
"RenderEvent 174 139 0 0 0 0 0\n"
"MouseMoveEvent 174 145 0 0 0 0 0\n"
"RenderEvent 174 145 0 0 0 0 0\n"
"MouseMoveEvent 174 148 0 0 0 0 0\n"
"RenderEvent 174 148 0 0 0 0 0\n"
"MouseMoveEvent 174 159 0 0 0 0 0\n"
"RenderEvent 174 159 0 0 0 0 0\n"
"MouseMoveEvent 174 162 0 0 0 0 0\n"
"RenderEvent 174 162 0 0 0 0 0\n"
"MouseMoveEvent 174 168 0 0 0 0 0\n"
"RenderEvent 174 168 0 0 0 0 0\n"
"MouseMoveEvent 174 174 0 0 0 0 0\n"
"RenderEvent 174 174 0 0 0 0 0\n"
"MouseMoveEvent 172 180 0 0 0 0 0\n"
"RenderEvent 172 180 0 0 0 0 0\n"
"MouseMoveEvent 170 188 0 0 0 0 0\n"
"RenderEvent 170 188 0 0 0 0 0\n"
"MouseMoveEvent 170 192 0 0 0 0 0\n"
"RenderEvent 170 192 0 0 0 0 0\n"
"MouseMoveEvent 168 198 0 0 0 0 0\n"
"RenderEvent 168 198 0 0 0 0 0\n"
"MouseMoveEvent 168 202 0 0 0 0 0\n"
"RenderEvent 168 202 0 0 0 0 0\n"
"MouseMoveEvent 167 207 0 0 0 0 0\n"
"RenderEvent 167 207 0 0 0 0 0\n"
"MouseMoveEvent 167 210 0 0 0 0 0\n"
"RenderEvent 167 210 0 0 0 0 0\n"
"MouseMoveEvent 167 216 0 0 0 0 0\n"
"RenderEvent 167 216 0 0 0 0 0\n"
"MouseMoveEvent 166 218 0 0 0 0 0\n"
"RenderEvent 166 218 0 0 0 0 0\n"
"MouseMoveEvent 166 229 0 0 0 0 0\n"
"RenderEvent 166 229 0 0 0 0 0\n"
"MouseMoveEvent 164 238 0 0 0 0 0\n"
"RenderEvent 164 238 0 0 0 0 0\n"
"MouseMoveEvent 164 241 0 0 0 0 0\n"
"RenderEvent 164 241 0 0 0 0 0\n"
"MouseMoveEvent 162 249 0 0 0 0 0\n"
"RenderEvent 162 249 0 0 0 0 0\n"
"MouseMoveEvent 160 257 0 0 0 0 0\n"
"RenderEvent 160 257 0 0 0 0 0\n"
"MouseMoveEvent 159 259 0 0 0 0 0\n"
"RenderEvent 159 259 0 0 0 0 0\n"
"MouseMoveEvent 159 262 0 0 0 0 0\n"
"RenderEvent 159 262 0 0 0 0 0\n"
"MouseMoveEvent 159 263 0 0 0 0 0\n"
"RenderEvent 159 263 0 0 0 0 0\n"
"MouseMoveEvent 159 264 0 0 0 0 0\n"
"RenderEvent 159 264 0 0 0 0 0\n"
"MouseMoveEvent 159 265 0 0 0 0 0\n"
"RenderEvent 159 265 0 0 0 0 0\n"
"MouseMoveEvent 160 265 0 0 0 0 0\n"
"RenderEvent 160 265 0 0 0 0 0\n"
"MouseMoveEvent 161 265 0 0 0 0 0\n"
"RenderEvent 161 265 0 0 0 0 0\n"
"MouseMoveEvent 161 266 0 0 0 0 0\n"
"RenderEvent 161 266 0 0 0 0 0\n"
"MouseMoveEvent 162 266 0 0 0 0 0\n"
"RenderEvent 162 266 0 0 0 0 0\n"
"MouseMoveEvent 162 267 0 0 0 0 0\n"
"RenderEvent 162 267 0 0 0 0 0\n"
"MouseMoveEvent 163 267 0 0 0 0 0\n"
"RenderEvent 163 267 0 0 0 0 0\n"
"LeftButtonReleaseEvent 163 267 0 0 0 0 0\n"
"MouseMoveEvent 164 267 0 0 0 0 0\n"
"MouseMoveEvent 154 160 0 0 0 0 0\n"
"RenderEvent 154 160 0 0 0 0 0\n"
"MouseMoveEvent 152 160 0 0 0 0 0\n"
"MouseMoveEvent 150 160 0 0 0 0 0\n"
"MouseMoveEvent 148 160 0 0 0 0 0\n"
"MouseMoveEvent 147 160 0 0 0 0 0\n"
"MouseMoveEvent 145 160 0 0 0 0 0\n"
"RenderEvent 145 160 0 0 0 0 0\n"
"MouseMoveEvent 145 159 0 0 0 0 0\n"
"MouseMoveEvent 146 159 0 0 0 0 0\n"
"RenderEvent 146 159 0 0 0 0 0\n"
"MouseMoveEvent 147 159 0 0 0 0 0\n"
"MouseMoveEvent 148 159 0 0 0 0 0\n"
"MouseMoveEvent 149 159 0 0 0 0 0\n"
"MouseMoveEvent 150 159 0 0 0 0 0\n"
"MouseMoveEvent 151 159 0 0 0 0 0\n"
"MouseMoveEvent 152 159 0 0 0 0 0\n"
"MiddleButtonPressEvent 152 159 0 0 0 0 0\n"
"RenderEvent 152 159 0 0 0 0 0\n"
"MouseMoveEvent 152 158 0 0 0 0 0\n"
"RenderEvent 152 158 0 0 0 0 0\n"
"MouseMoveEvent 151 157 0 0 0 0 0\n"
"RenderEvent 151 157 0 0 0 0 0\n"
"MouseMoveEvent 151 156 0 0 0 0 0\n"
"RenderEvent 151 156 0 0 0 0 0\n"
"MouseMoveEvent 150 156 0 0 0 0 0\n"
"RenderEvent 150 156 0 0 0 0 0\n"
"MouseMoveEvent 149 155 0 0 0 0 0\n"
"RenderEvent 149 155 0 0 0 0 0\n"
"MouseMoveEvent 148 154 0 0 0 0 0\n"
"RenderEvent 148 154 0 0 0 0 0\n"
"MouseMoveEvent 147 153 0 0 0 0 0\n"
"RenderEvent 147 153 0 0 0 0 0\n"
"MouseMoveEvent 146 152 0 0 0 0 0\n"
"RenderEvent 146 152 0 0 0 0 0\n"
"MouseMoveEvent 145 151 0 0 0 0 0\n"
"RenderEvent 145 151 0 0 0 0 0\n"
"MouseMoveEvent 145 150 0 0 0 0 0\n"
"RenderEvent 145 150 0 0 0 0 0\n"
"MouseMoveEvent 144 150 0 0 0 0 0\n"
"RenderEvent 144 150 0 0 0 0 0\n"
"MouseMoveEvent 144 149 0 0 0 0 0\n"
"RenderEvent 144 149 0 0 0 0 0\n"
"MouseMoveEvent 143 149 0 0 0 0 0\n"
"RenderEvent 143 149 0 0 0 0 0\n"
"MouseMoveEvent 143 148 0 0 0 0 0\n"
"RenderEvent 143 148 0 0 0 0 0\n"
"MouseMoveEvent 142 147 0 0 0 0 0\n"
"RenderEvent 142 147 0 0 0 0 0\n"
"MouseMoveEvent 141 147 0 0 0 0 0\n"
"RenderEvent 141 147 0 0 0 0 0\n"
"MouseMoveEvent 140 146 0 0 0 0 0\n"
"RenderEvent 140 146 0 0 0 0 0\n"
"MouseMoveEvent 139 145 0 0 0 0 0\n"
"RenderEvent 139 145 0 0 0 0 0\n"
"MouseMoveEvent 138 145 0 0 0 0 0\n"
"RenderEvent 138 145 0 0 0 0 0\n"
"MouseMoveEvent 138 144 0 0 0 0 0\n"
"RenderEvent 138 144 0 0 0 0 0\n"
"MouseMoveEvent 137 144 0 0 0 0 0\n"
"RenderEvent 137 144 0 0 0 0 0\n"
"MouseMoveEvent 136 144 0 0 0 0 0\n"
"RenderEvent 136 144 0 0 0 0 0\n"
"MouseMoveEvent 135 143 0 0 0 0 0\n"
"RenderEvent 135 143 0 0 0 0 0\n"
"MouseMoveEvent 134 143 0 0 0 0 0\n"
"RenderEvent 134 143 0 0 0 0 0\n"
"MouseMoveEvent 132 143 0 0 0 0 0\n"
"RenderEvent 132 143 0 0 0 0 0\n"
"MouseMoveEvent 131 142 0 0 0 0 0\n"
"RenderEvent 131 142 0 0 0 0 0\n"
"MouseMoveEvent 130 142 0 0 0 0 0\n"
"RenderEvent 130 142 0 0 0 0 0\n"
"MouseMoveEvent 129 141 0 0 0 0 0\n"
"RenderEvent 129 141 0 0 0 0 0\n"
"MouseMoveEvent 127 141 0 0 0 0 0\n"
"RenderEvent 127 141 0 0 0 0 0\n"
"MouseMoveEvent 124 140 0 0 0 0 0\n"
"RenderEvent 124 140 0 0 0 0 0\n"
"MouseMoveEvent 122 140 0 0 0 0 0\n"
"RenderEvent 122 140 0 0 0 0 0\n"
"MouseMoveEvent 121 140 0 0 0 0 0\n"
"RenderEvent 121 140 0 0 0 0 0\n"
"MouseMoveEvent 118 140 0 0 0 0 0\n"
"RenderEvent 118 140 0 0 0 0 0\n"
"MouseMoveEvent 116 140 0 0 0 0 0\n"
"RenderEvent 116 140 0 0 0 0 0\n"
"MouseMoveEvent 114 140 0 0 0 0 0\n"
"RenderEvent 114 140 0 0 0 0 0\n"
"MouseMoveEvent 109 139 0 0 0 0 0\n"
"RenderEvent 109 139 0 0 0 0 0\n"
"MouseMoveEvent 105 139 0 0 0 0 0\n"
"RenderEvent 105 139 0 0 0 0 0\n"
"MouseMoveEvent 101 139 0 0 0 0 0\n"
"RenderEvent 101 139 0 0 0 0 0\n"
"MouseMoveEvent 96 140 0 0 0 0 0\n"
"RenderEvent 96 140 0 0 0 0 0\n"
"MouseMoveEvent 95 140 0 0 0 0 0\n"
"RenderEvent 95 140 0 0 0 0 0\n"
"MouseMoveEvent 94 140 0 0 0 0 0\n"
"RenderEvent 94 140 0 0 0 0 0\n"
"MouseMoveEvent 94 141 0 0 0 0 0\n"
"RenderEvent 94 141 0 0 0 0 0\n"
"MouseMoveEvent 92 141 0 0 0 0 0\n"
"RenderEvent 92 141 0 0 0 0 0\n"
"MouseMoveEvent 91 142 0 0 0 0 0\n"
"RenderEvent 91 142 0 0 0 0 0\n"
"MouseMoveEvent 89 143 0 0 0 0 0\n"
"RenderEvent 89 143 0 0 0 0 0\n"
"MouseMoveEvent 87 144 0 0 0 0 0\n"
"RenderEvent 87 144 0 0 0 0 0\n"
"MouseMoveEvent 86 145 0 0 0 0 0\n"
"RenderEvent 86 145 0 0 0 0 0\n"
"MouseMoveEvent 84 145 0 0 0 0 0\n"
"RenderEvent 84 145 0 0 0 0 0\n"
"MouseMoveEvent 82 146 0 0 0 0 0\n"
"RenderEvent 82 146 0 0 0 0 0\n"
"MouseMoveEvent 78 147 0 0 0 0 0\n"
"RenderEvent 78 147 0 0 0 0 0\n"
"MouseMoveEvent 76 149 0 0 0 0 0\n"
"RenderEvent 76 149 0 0 0 0 0\n"
"MouseMoveEvent 74 149 0 0 0 0 0\n"
"RenderEvent 74 149 0 0 0 0 0\n"
"MouseMoveEvent 73 149 0 0 0 0 0\n"
"RenderEvent 73 149 0 0 0 0 0\n"
"MouseMoveEvent 84 154 0 0 0 0 0\n"
"RenderEvent 84 154 0 0 0 0 0\n"
"MouseMoveEvent 83 155 0 0 0 0 0\n"
"RenderEvent 83 155 0 0 0 0 0\n"
"MouseMoveEvent 80 155 0 0 0 0 0\n"
"RenderEvent 80 155 0 0 0 0 0\n"
"MouseMoveEvent 75 156 0 0 0 0 0\n"
"RenderEvent 75 156 0 0 0 0 0\n"
"MouseMoveEvent 73 156 0 0 0 0 0\n"
"RenderEvent 73 156 0 0 0 0 0\n"
"MouseMoveEvent 72 156 0 0 0 0 0\n"
"RenderEvent 72 156 0 0 0 0 0\n"
"MouseMoveEvent 70 156 0 0 0 0 0\n"
"RenderEvent 70 156 0 0 0 0 0\n"
"MouseMoveEvent 68 157 0 0 0 0 0\n"
"RenderEvent 68 157 0 0 0 0 0\n"
"MouseMoveEvent 66 157 0 0 0 0 0\n"
"RenderEvent 66 157 0 0 0 0 0\n"
"MouseMoveEvent 65 157 0 0 0 0 0\n"
"RenderEvent 65 157 0 0 0 0 0\n"
"MouseMoveEvent 64 157 0 0 0 0 0\n"
"RenderEvent 64 157 0 0 0 0 0\n"
"MouseMoveEvent 64 156 0 0 0 0 0\n"
"RenderEvent 64 156 0 0 0 0 0\n"
"MouseMoveEvent 63 156 0 0 0 0 0\n"
"RenderEvent 63 156 0 0 0 0 0\n"
"MouseMoveEvent 62 156 0 0 0 0 0\n"
"RenderEvent 62 156 0 0 0 0 0\n"
"MouseMoveEvent 61 155 0 0 0 0 0\n"
"RenderEvent 61 155 0 0 0 0 0\n"
"MouseMoveEvent 59 154 0 0 0 0 0\n"
"RenderEvent 59 154 0 0 0 0 0\n"
"MouseMoveEvent 59 153 0 0 0 0 0\n"
"RenderEvent 59 153 0 0 0 0 0\n"
"MouseMoveEvent 58 152 0 0 0 0 0\n"
"RenderEvent 58 152 0 0 0 0 0\n"
"MouseMoveEvent 57 151 0 0 0 0 0\n"
"RenderEvent 57 151 0 0 0 0 0\n"
"MouseMoveEvent 56 149 0 0 0 0 0\n"
"RenderEvent 56 149 0 0 0 0 0\n"
"MouseMoveEvent 55 148 0 0 0 0 0\n"
"RenderEvent 55 148 0 0 0 0 0\n"
"MouseMoveEvent 54 146 0 0 0 0 0\n"
"RenderEvent 54 146 0 0 0 0 0\n"
"MouseMoveEvent 53 145 0 0 0 0 0\n"
"RenderEvent 53 145 0 0 0 0 0\n"
"MouseMoveEvent 53 144 0 0 0 0 0\n"
"RenderEvent 53 144 0 0 0 0 0\n"
"MouseMoveEvent 53 143 0 0 0 0 0\n"
"RenderEvent 53 143 0 0 0 0 0\n"
"MouseMoveEvent 53 142 0 0 0 0 0\n"
"RenderEvent 53 142 0 0 0 0 0\n"
"MouseMoveEvent 53 141 0 0 0 0 0\n"
"RenderEvent 53 141 0 0 0 0 0\n"
"MouseMoveEvent 53 140 0 0 0 0 0\n"
"RenderEvent 53 140 0 0 0 0 0\n"
"MouseMoveEvent 53 139 0 0 0 0 0\n"
"RenderEvent 53 139 0 0 0 0 0\n"
"MouseMoveEvent 53 138 0 0 0 0 0\n"
"RenderEvent 53 138 0 0 0 0 0\n"
"MouseMoveEvent 52 136 0 0 0 0 0\n"
"RenderEvent 52 136 0 0 0 0 0\n"
"MouseMoveEvent 51 135 0 0 0 0 0\n"
"RenderEvent 51 135 0 0 0 0 0\n"
"MouseMoveEvent 51 134 0 0 0 0 0\n"
"RenderEvent 51 134 0 0 0 0 0\n"
"MouseMoveEvent 51 133 0 0 0 0 0\n"
"RenderEvent 51 133 0 0 0 0 0\n"
"MouseMoveEvent 51 131 0 0 0 0 0\n"
"RenderEvent 51 131 0 0 0 0 0\n"
"MouseMoveEvent 51 130 0 0 0 0 0\n"
"RenderEvent 51 130 0 0 0 0 0\n"
"MouseMoveEvent 52 130 0 0 0 0 0\n"
"RenderEvent 52 130 0 0 0 0 0\n"
"MouseMoveEvent 53 131 0 0 0 0 0\n"
"RenderEvent 53 131 0 0 0 0 0\n"
"MouseMoveEvent 53 132 0 0 0 0 0\n"
"RenderEvent 53 132 0 0 0 0 0\n"
"MouseMoveEvent 53 133 0 0 0 0 0\n"
"RenderEvent 53 133 0 0 0 0 0\n"
"MouseMoveEvent 54 136 0 0 0 0 0\n"
"RenderEvent 54 136 0 0 0 0 0\n"
"MouseMoveEvent 55 138 0 0 0 0 0\n"
"RenderEvent 55 138 0 0 0 0 0\n"
"MouseMoveEvent 55 143 0 0 0 0 0\n"
"RenderEvent 55 143 0 0 0 0 0\n"
"MouseMoveEvent 55 147 0 0 0 0 0\n"
"RenderEvent 55 147 0 0 0 0 0\n"
"MouseMoveEvent 55 154 0 0 0 0 0\n"
"RenderEvent 55 154 0 0 0 0 0\n"
"MouseMoveEvent 55 157 0 0 0 0 0\n"
"RenderEvent 55 157 0 0 0 0 0\n"
"MouseMoveEvent 55 161 0 0 0 0 0\n"
"RenderEvent 55 161 0 0 0 0 0\n"
"MouseMoveEvent 56 166 0 0 0 0 0\n"
"RenderEvent 56 166 0 0 0 0 0\n"
"MouseMoveEvent 56 170 0 0 0 0 0\n"
"RenderEvent 56 170 0 0 0 0 0\n"
"MouseMoveEvent 56 174 0 0 0 0 0\n"
"RenderEvent 56 174 0 0 0 0 0\n"
"MouseMoveEvent 57 180 0 0 0 0 0\n"
"RenderEvent 57 180 0 0 0 0 0\n"
"MouseMoveEvent 62 190 0 0 0 0 0\n"
"RenderEvent 62 190 0 0 0 0 0\n"
"MouseMoveEvent 64 198 0 0 0 0 0\n"
"RenderEvent 64 198 0 0 0 0 0\n"
"MouseMoveEvent 69 206 0 0 0 0 0\n"
"RenderEvent 69 206 0 0 0 0 0\n"
"MouseMoveEvent 71 210 0 0 0 0 0\n"
"RenderEvent 71 210 0 0 0 0 0\n"
"MouseMoveEvent 67 204 0 0 0 0 0\n"
"RenderEvent 67 204 0 0 0 0 0\n"
"MouseMoveEvent 64 209 0 0 0 0 0\n"
"RenderEvent 64 209 0 0 0 0 0\n"
"MouseMoveEvent 73 216 0 0 0 0 0\n"
"RenderEvent 73 216 0 0 0 0 0\n"
"MouseMoveEvent 87 224 0 0 0 0 0\n"
"RenderEvent 87 224 0 0 0 0 0\n"
"MouseMoveEvent 113 233 0 0 0 0 0\n"
"RenderEvent 113 233 0 0 0 0 0\n"
"MouseMoveEvent 157 239 0 0 0 0 0\n"
"RenderEvent 157 239 0 0 0 0 0\n"
"MouseMoveEvent 186 241 0 0 0 0 0\n"
"RenderEvent 186 241 0 0 0 0 0\n"
"MouseMoveEvent 215 240 0 0 0 0 0\n"
"RenderEvent 215 240 0 0 0 0 0\n"
"MouseMoveEvent 237 236 0 0 0 0 0\n"
"RenderEvent 237 236 0 0 0 0 0\n"
"MouseMoveEvent 251 235 0 0 0 0 0\n"
"RenderEvent 251 235 0 0 0 0 0\n"
"MouseMoveEvent 257 235 0 0 0 0 0\n"
"RenderEvent 257 235 0 0 0 0 0\n"
"MouseMoveEvent 258 235 0 0 0 0 0\n"
"RenderEvent 258 235 0 0 0 0 0\n"
"MouseMoveEvent 259 234 0 0 0 0 0\n"
"RenderEvent 259 234 0 0 0 0 0\n"
"MouseMoveEvent 260 231 0 0 0 0 0\n"
"RenderEvent 260 231 0 0 0 0 0\n"
"MouseMoveEvent 261 228 0 0 0 0 0\n"
"RenderEvent 261 228 0 0 0 0 0\n"
"MouseMoveEvent 261 221 0 0 0 0 0\n"
"RenderEvent 261 221 0 0 0 0 0\n"
"MouseMoveEvent 262 216 0 0 0 0 0\n"
"RenderEvent 262 216 0 0 0 0 0\n"
"MouseMoveEvent 263 209 0 0 0 0 0\n"
"RenderEvent 263 209 0 0 0 0 0\n"
"MouseMoveEvent 264 204 0 0 0 0 0\n"
"RenderEvent 264 204 0 0 0 0 0\n"
"MouseMoveEvent 265 200 0 0 0 0 0\n"
"RenderEvent 265 200 0 0 0 0 0\n"
"MouseMoveEvent 265 198 0 0 0 0 0\n"
"RenderEvent 265 198 0 0 0 0 0\n"
"MouseMoveEvent 265 196 0 0 0 0 0\n"
"RenderEvent 265 196 0 0 0 0 0\n"
"MouseMoveEvent 265 195 0 0 0 0 0\n"
"RenderEvent 265 195 0 0 0 0 0\n"
"MouseMoveEvent 264 193 0 0 0 0 0\n"
"RenderEvent 264 193 0 0 0 0 0\n"
"MouseMoveEvent 263 191 0 0 0 0 0\n"
"RenderEvent 263 191 0 0 0 0 0\n"
"MouseMoveEvent 262 187 0 0 0 0 0\n"
"RenderEvent 262 187 0 0 0 0 0\n"
"MouseMoveEvent 260 184 0 0 0 0 0\n"
"RenderEvent 260 184 0 0 0 0 0\n"
"MouseMoveEvent 259 181 0 0 0 0 0\n"
"RenderEvent 259 181 0 0 0 0 0\n"
"MouseMoveEvent 257 178 0 0 0 0 0\n"
"RenderEvent 257 178 0 0 0 0 0\n"
"MouseMoveEvent 256 177 0 0 0 0 0\n"
"RenderEvent 256 177 0 0 0 0 0\n"
"MouseMoveEvent 255 174 0 0 0 0 0\n"
"RenderEvent 255 174 0 0 0 0 0\n"
"MouseMoveEvent 248 169 0 0 0 0 0\n"
"RenderEvent 248 169 0 0 0 0 0\n"
"MouseMoveEvent 246 166 0 0 0 0 0\n"
"RenderEvent 246 166 0 0 0 0 0\n"
"MouseMoveEvent 240 161 0 0 0 0 0\n"
"RenderEvent 240 161 0 0 0 0 0\n"
"MouseMoveEvent 235 156 0 0 0 0 0\n"
"RenderEvent 235 156 0 0 0 0 0\n"
"MouseMoveEvent 230 153 0 0 0 0 0\n"
"RenderEvent 230 153 0 0 0 0 0\n"
"MouseMoveEvent 227 151 0 0 0 0 0\n"
"RenderEvent 227 151 0 0 0 0 0\n"
"MouseMoveEvent 224 149 0 0 0 0 0\n"
"RenderEvent 224 149 0 0 0 0 0\n"
"MouseMoveEvent 219 146 0 0 0 0 0\n"
"RenderEvent 219 146 0 0 0 0 0\n"
"MouseMoveEvent 215 144 0 0 0 0 0\n"
"RenderEvent 215 144 0 0 0 0 0\n"
"MouseMoveEvent 212 143 0 0 0 0 0\n"
"RenderEvent 212 143 0 0 0 0 0\n"
"MouseMoveEvent 209 142 0 0 0 0 0\n"
"RenderEvent 209 142 0 0 0 0 0\n"
"MouseMoveEvent 207 142 0 0 0 0 0\n"
"RenderEvent 207 142 0 0 0 0 0\n"
"MouseMoveEvent 205 142 0 0 0 0 0\n"
"RenderEvent 205 142 0 0 0 0 0\n"
"MouseMoveEvent 203 142 0 0 0 0 0\n"
"RenderEvent 203 142 0 0 0 0 0\n"
"MouseMoveEvent 199 142 0 0 0 0 0\n"
"RenderEvent 199 142 0 0 0 0 0\n"
"MouseMoveEvent 197 142 0 0 0 0 0\n"
"RenderEvent 197 142 0 0 0 0 0\n"
"MouseMoveEvent 195 142 0 0 0 0 0\n"
"RenderEvent 195 142 0 0 0 0 0\n"
"MouseMoveEvent 193 142 0 0 0 0 0\n"
"RenderEvent 193 142 0 0 0 0 0\n"
"MouseMoveEvent 192 142 0 0 0 0 0\n"
"RenderEvent 192 142 0 0 0 0 0\n"
"MouseMoveEvent 190 142 0 0 0 0 0\n"
"RenderEvent 190 142 0 0 0 0 0\n"
"MouseMoveEvent 188 142 0 0 0 0 0\n"
"RenderEvent 188 142 0 0 0 0 0\n"
"MouseMoveEvent 184 142 0 0 0 0 0\n"
"RenderEvent 184 142 0 0 0 0 0\n"
"MouseMoveEvent 182 143 0 0 0 0 0\n"
"RenderEvent 182 143 0 0 0 0 0\n"
"MouseMoveEvent 180 143 0 0 0 0 0\n"
"RenderEvent 180 143 0 0 0 0 0\n"
"MouseMoveEvent 180 144 0 0 0 0 0\n"
"RenderEvent 180 144 0 0 0 0 0\n"
"MouseMoveEvent 179 144 0 0 0 0 0\n"
"RenderEvent 179 144 0 0 0 0 0\n"
"MouseMoveEvent 178 144 0 0 0 0 0\n"
"RenderEvent 178 144 0 0 0 0 0\n"
"MouseMoveEvent 177 145 0 0 0 0 0\n"
"RenderEvent 177 145 0 0 0 0 0\n"
"MouseMoveEvent 176 145 0 0 0 0 0\n"
"RenderEvent 176 145 0 0 0 0 0\n"
"MouseMoveEvent 174 145 0 0 0 0 0\n"
"RenderEvent 174 145 0 0 0 0 0\n"
"MouseMoveEvent 172 145 0 0 0 0 0\n"
"RenderEvent 172 145 0 0 0 0 0\n"
"MouseMoveEvent 169 146 0 0 0 0 0\n"
"RenderEvent 169 146 0 0 0 0 0\n"
"MouseMoveEvent 167 146 0 0 0 0 0\n"
"RenderEvent 167 146 0 0 0 0 0\n"
"MouseMoveEvent 166 146 0 0 0 0 0\n"
"RenderEvent 166 146 0 0 0 0 0\n"
"MouseMoveEvent 165 147 0 0 0 0 0\n"
"RenderEvent 165 147 0 0 0 0 0\n"
"MouseMoveEvent 164 148 0 0 0 0 0\n"
"RenderEvent 164 148 0 0 0 0 0\n"
"MouseMoveEvent 163 149 0 0 0 0 0\n"
"RenderEvent 163 149 0 0 0 0 0\n"
"MouseMoveEvent 162 150 0 0 0 0 0\n"
"RenderEvent 162 150 0 0 0 0 0\n"
"MouseMoveEvent 161 150 0 0 0 0 0\n"
"RenderEvent 161 150 0 0 0 0 0\n"
"MouseMoveEvent 160 150 0 0 0 0 0\n"
"RenderEvent 160 150 0 0 0 0 0\n"
"MouseMoveEvent 160 151 0 0 0 0 0\n"
"RenderEvent 160 151 0 0 0 0 0\n"
"MouseMoveEvent 158 151 0 0 0 0 0\n"
"RenderEvent 158 151 0 0 0 0 0\n"
"MouseMoveEvent 157 152 0 0 0 0 0\n"
"RenderEvent 157 152 0 0 0 0 0\n"
"MouseMoveEvent 156 152 0 0 0 0 0\n"
"RenderEvent 156 152 0 0 0 0 0\n"
"MouseMoveEvent 153 152 0 0 0 0 0\n"
"RenderEvent 153 152 0 0 0 0 0\n"
"MouseMoveEvent 149 153 0 0 0 0 0\n"
"RenderEvent 149 153 0 0 0 0 0\n"
"MouseMoveEvent 139 153 0 0 0 0 0\n"
"RenderEvent 139 153 0 0 0 0 0\n"
"MouseMoveEvent 135 153 0 0 0 0 0\n"
"RenderEvent 135 153 0 0 0 0 0\n"
"MouseMoveEvent 134 153 0 0 0 0 0\n"
"RenderEvent 134 153 0 0 0 0 0\n"
"MouseMoveEvent 133 153 0 0 0 0 0\n"
"RenderEvent 133 153 0 0 0 0 0\n"
"MouseMoveEvent 132 154 0 0 0 0 0\n"
"RenderEvent 132 154 0 0 0 0 0\n"
"MouseMoveEvent 129 156 0 0 0 0 0\n"
"RenderEvent 129 156 0 0 0 0 0\n"
"MouseMoveEvent 126 157 0 0 0 0 0\n"
"RenderEvent 126 157 0 0 0 0 0\n"
"MouseMoveEvent 121 159 0 0 0 0 0\n"
"RenderEvent 121 159 0 0 0 0 0\n"
"MouseMoveEvent 117 161 0 0 0 0 0\n"
"RenderEvent 117 161 0 0 0 0 0\n"
"MouseMoveEvent 113 162 0 0 0 0 0\n"
"RenderEvent 113 162 0 0 0 0 0\n"
"MouseMoveEvent 108 164 0 0 0 0 0\n"
"RenderEvent 108 164 0 0 0 0 0\n"
"MouseMoveEvent 105 166 0 0 0 0 0\n"
"RenderEvent 105 166 0 0 0 0 0\n"
"MouseMoveEvent 102 167 0 0 0 0 0\n"
"RenderEvent 102 167 0 0 0 0 0\n"
"MouseMoveEvent 100 168 0 0 0 0 0\n"
"RenderEvent 100 168 0 0 0 0 0\n"
"MouseMoveEvent 96 170 0 0 0 0 0\n"
"RenderEvent 96 170 0 0 0 0 0\n"
"MouseMoveEvent 95 170 0 0 0 0 0\n"
"RenderEvent 95 170 0 0 0 0 0\n"
"MouseMoveEvent 94 170 0 0 0 0 0\n"
"RenderEvent 94 170 0 0 0 0 0\n"
"MouseMoveEvent 91 172 0 0 0 0 0\n"
"RenderEvent 91 172 0 0 0 0 0\n"
"MouseMoveEvent 89 173 0 0 0 0 0\n"
"RenderEvent 89 173 0 0 0 0 0\n"
"MouseMoveEvent 86 174 0 0 0 0 0\n"
"RenderEvent 86 174 0 0 0 0 0\n"
"MouseMoveEvent 74 180 0 0 0 0 0\n"
"RenderEvent 74 180 0 0 0 0 0\n"
"MouseMoveEvent 71 182 0 0 0 0 0\n"
"RenderEvent 71 182 0 0 0 0 0\n"
"MouseMoveEvent 66 186 0 0 0 0 0\n"
"RenderEvent 66 186 0 0 0 0 0\n"
"MouseMoveEvent 65 186 0 0 0 0 0\n"
"RenderEvent 65 186 0 0 0 0 0\n"
"MouseMoveEvent 64 187 0 0 0 0 0\n"
"RenderEvent 64 187 0 0 0 0 0\n"
"MouseMoveEvent 64 189 0 0 0 0 0\n"
"RenderEvent 64 189 0 0 0 0 0\n"
"MouseMoveEvent 62 191 0 0 0 0 0\n"
"RenderEvent 62 191 0 0 0 0 0\n"
"MouseMoveEvent 58 193 0 0 0 0 0\n"
"RenderEvent 58 193 0 0 0 0 0\n"
"MouseMoveEvent 57 194 0 0 0 0 0\n"
"RenderEvent 57 194 0 0 0 0 0\n"
"MouseMoveEvent 56 195 0 0 0 0 0\n"
"RenderEvent 56 195 0 0 0 0 0\n"
"MouseMoveEvent 55 195 0 0 0 0 0\n"
"RenderEvent 55 195 0 0 0 0 0\n"
"MouseMoveEvent 54 195 0 0 0 0 0\n"
"RenderEvent 54 195 0 0 0 0 0\n"
"MouseMoveEvent 54 194 0 0 0 0 0\n"
"RenderEvent 54 194 0 0 0 0 0\n"
"MouseMoveEvent 54 192 0 0 0 0 0\n"
"RenderEvent 54 192 0 0 0 0 0\n"
"MouseMoveEvent 53 191 0 0 0 0 0\n"
"RenderEvent 53 191 0 0 0 0 0\n"
"MouseMoveEvent 52 189 0 0 0 0 0\n"
"RenderEvent 52 189 0 0 0 0 0\n"
"MouseMoveEvent 52 188 0 0 0 0 0\n"
"RenderEvent 52 188 0 0 0 0 0\n"
"MouseMoveEvent 52 187 0 0 0 0 0\n"
"RenderEvent 52 187 0 0 0 0 0\n"
"MouseMoveEvent 52 186 0 0 0 0 0\n"
"RenderEvent 52 186 0 0 0 0 0\n"
"MouseMoveEvent 52 185 0 0 0 0 0\n"
"RenderEvent 52 185 0 0 0 0 0\n"
"MouseMoveEvent 52 184 0 0 0 0 0\n"
"RenderEvent 52 184 0 0 0 0 0\n"
"MouseMoveEvent 52 183 0 0 0 0 0\n"
"RenderEvent 52 183 0 0 0 0 0\n"
"MouseMoveEvent 52 182 0 0 0 0 0\n"
"RenderEvent 52 182 0 0 0 0 0\n"
"MouseMoveEvent 52 181 0 0 0 0 0\n"
"RenderEvent 52 181 0 0 0 0 0\n"
"MouseMoveEvent 53 180 0 0 0 0 0\n"
"RenderEvent 53 180 0 0 0 0 0\n"
"MouseMoveEvent 54 178 0 0 0 0 0\n"
"RenderEvent 54 178 0 0 0 0 0\n"
"MouseMoveEvent 54 177 0 0 0 0 0\n"
"RenderEvent 54 177 0 0 0 0 0\n"
"MouseMoveEvent 54 176 0 0 0 0 0\n"
"RenderEvent 54 176 0 0 0 0 0\n"
"MouseMoveEvent 55 176 0 0 0 0 0\n"
"RenderEvent 55 176 0 0 0 0 0\n"
"MouseMoveEvent 55 175 0 0 0 0 0\n"
"RenderEvent 55 175 0 0 0 0 0\n"
"MouseMoveEvent 56 174 0 0 0 0 0\n"
"RenderEvent 56 174 0 0 0 0 0\n"
"MouseMoveEvent 56 173 0 0 0 0 0\n"
"RenderEvent 56 173 0 0 0 0 0\n"
"MouseMoveEvent 57 173 0 0 0 0 0\n"
"RenderEvent 57 173 0 0 0 0 0\n"
"MouseMoveEvent 57 172 0 0 0 0 0\n"
"RenderEvent 57 172 0 0 0 0 0\n"
"MouseMoveEvent 57 171 0 0 0 0 0\n"
"RenderEvent 57 171 0 0 0 0 0\n"
"MouseMoveEvent 57 169 0 0 0 0 0\n"
"RenderEvent 57 169 0 0 0 0 0\n"
"MouseMoveEvent 58 167 0 0 0 0 0\n"
"RenderEvent 58 167 0 0 0 0 0\n"
"MouseMoveEvent 58 166 0 0 0 0 0\n"
"RenderEvent 58 166 0 0 0 0 0\n"
"MouseMoveEvent 59 165 0 0 0 0 0\n"
"RenderEvent 59 165 0 0 0 0 0\n"
"MouseMoveEvent 58 165 0 0 0 0 0\n"
"RenderEvent 58 165 0 0 0 0 0\n"
"MouseMoveEvent 59 165 0 0 0 0 0\n"
"RenderEvent 59 165 0 0 0 0 0\n"
"MouseMoveEvent 56 165 0 0 0 0 0\n"
"RenderEvent 56 165 0 0 0 0 0\n"
"MouseMoveEvent 56 164 0 0 0 0 0\n"
"RenderEvent 56 164 0 0 0 0 0\n"
"MouseMoveEvent 57 164 0 0 0 0 0\n"
"RenderEvent 57 164 0 0 0 0 0\n"
"MouseMoveEvent 58 163 0 0 0 0 0\n"
"RenderEvent 58 163 0 0 0 0 0\n"
"MouseMoveEvent 59 163 0 0 0 0 0\n"
"RenderEvent 59 163 0 0 0 0 0\n"
"MouseMoveEvent 60 163 0 0 0 0 0\n"
"RenderEvent 60 163 0 0 0 0 0\n"
"MouseMoveEvent 62 162 0 0 0 0 0\n"
"RenderEvent 62 162 0 0 0 0 0\n"
"MouseMoveEvent 63 162 0 0 0 0 0\n"
"RenderEvent 63 162 0 0 0 0 0\n"
"MouseMoveEvent 64 162 0 0 0 0 0\n"
"RenderEvent 64 162 0 0 0 0 0\n"
"MouseMoveEvent 65 161 0 0 0 0 0\n"
"RenderEvent 65 161 0 0 0 0 0\n"
"MouseMoveEvent 66 161 0 0 0 0 0\n"
"RenderEvent 66 161 0 0 0 0 0\n"
"MouseMoveEvent 65 161 0 0 0 0 0\n"
"RenderEvent 65 161 0 0 0 0 0\n"
"MouseMoveEvent 65 160 0 0 0 0 0\n"
"RenderEvent 65 160 0 0 0 0 0\n"
"MouseMoveEvent 64 160 0 0 0 0 0\n"
"RenderEvent 64 160 0 0 0 0 0\n"
"MouseMoveEvent 63 160 0 0 0 0 0\n"
"RenderEvent 63 160 0 0 0 0 0\n"
"MiddleButtonReleaseEvent 63 160 0 0 0 0 0\n"
"RenderEvent 63 160 0 0 0 0 0\n"
"MouseMoveEvent 63 159 0 0 0 0 0\n"
"RenderEvent 63 159 0 0 0 0 0\n"
"MouseMoveEvent 63 158 0 0 0 0 0\n"
"MouseMoveEvent 63 157 0 0 0 0 0\n"
"MouseMoveEvent 63 155 0 0 0 0 0\n"
"MouseMoveEvent 63 152 0 0 0 0 0\n"
"MouseMoveEvent 65 146 0 0 0 0 0\n"
"MouseMoveEvent 65 138 0 0 0 0 0\n"
"RenderEvent 65 138 0 0 0 0 0\n"
  ;

int TestCaptionWidget( int argc, char *argv[] )
{
  // Create the RenderWindow, Renderer and both Actors
  //
  vtkSmartPointer<vtkRenderer> ren1 =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);

  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);

  // Create a test pipeline
  //
  vtkSmartPointer<vtkSphereSource> ss = vtkSmartPointer<vtkSphereSource>::New();
  ss->SetCenter(100,250,500);
  ss->Update();

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(ss->GetOutput());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // Create the widget
  vtkSmartPointer<vtkCaptionRepresentation> rep =
    vtkSmartPointer<vtkCaptionRepresentation>::New();
  rep->SetAnchorPosition(ss->GetOutput()->GetPoint(0));
  double startingPosition[3];
  rep->GetAnchorPosition(startingPosition);
  rep->GetCaptionActor2D()->
    SetCaption("This is a test caption\nAnd it has two lines");
  rep->GetCaptionActor2D()->
    GetTextActor()->GetTextProperty()->SetJustificationToCentered();
  rep->GetCaptionActor2D()->
    GetTextActor()->GetTextProperty()->SetVerticalJustificationToCentered();

  vtkSmartPointer<vtkCaptionWidget> widget =
    vtkSmartPointer<vtkCaptionWidget>::New();
  widget->SetInteractor(iren);
  widget->SetRepresentation(rep);

  // Print the widget and its representation
  rep->Print(std::cout);
  widget->Print(std::cout);

  // Add the actors to the renderer, set the background and size
  //
  ren1->AddActor(actor);
  ren1->SetBackground(0.1, 0.2, 0.4);
  renWin->SetSize(300, 300);

  // record events
  vtkSmartPointer<vtkInteractorEventRecorder> recorder =
    vtkSmartPointer<vtkInteractorEventRecorder>::New();
  recorder->SetInteractor(iren);

#ifdef RECORD
  recorder->SetFileName("record.log");
  recorder->On();
  recorder->Record();
#else
  recorder->ReadFromInputStringOn();
  recorder->SetInputString(eventLog);
#endif

  // render the image
  //
  iren->Initialize();
  renWin->Render();
  widget->On();
  renWin->Render();

#ifndef RECORD
  recorder->Play();
  recorder->Off();
#endif

  cout << "Setting new caption\n";
  rep->GetCaptionActor2D()->SetCaption("Okay the caption has now changed and the border should resize");
  rep->Modified();
  renWin->Render();

  iren->Start();

  double endingPosition[3];
  rep->GetAnchorPosition(endingPosition);
  std::cout << "Starting position of anchor: "
       << startingPosition[0] << ", "
       << startingPosition[1] << ", "
       << startingPosition[2] << std::endl;
  std::cout << "Ending position of anchor: "
       << endingPosition[0] << ", "
       << endingPosition[1] << ", "
       << endingPosition[2] << std::endl;

  return EXIT_SUCCESS;
}
