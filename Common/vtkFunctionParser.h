/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFunctionParser.h
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
// .NAME vtkFunctionParser - Parse and evaluate a mathematical expression
// .SECTION Description
// vtkFunctionParser is a class that takes in a mathematical expression as
// a char string, parses it, and evaluates it at the specified values of
// the variables in the input string.

#ifndef __vtkFunctionParser_h
#define __vtkFunctionParser_h

#include "vtkObject.h"

#define VTK_PARSER_IMMEDIATE 1
#define VTK_PARSER_UNARY_MINUS 2

// supported math functions
#define VTK_PARSER_ADD 3
#define VTK_PARSER_SUBTRACT 4
#define VTK_PARSER_MULTIPLY 5
#define VTK_PARSER_DIVIDE 6
#define VTK_PARSER_POWER 7
#define VTK_PARSER_ABSOLUTE_VALUE 8
#define VTK_PARSER_EXPONENT 9
#define VTK_PARSER_CEILING 10
#define VTK_PARSER_FLOOR 11
#define VTK_PARSER_LOGARITHM 12
#define VTK_PARSER_SQUARE_ROOT 13
#define VTK_PARSER_SINE 14
#define VTK_PARSER_COSINE 15
#define VTK_PARSER_TANGENT 16
#define VTK_PARSER_ARCSINE 17
#define VTK_PARSER_ARCCOSINE 18
#define VTK_PARSER_ARCTANGENT 19
#define VTK_PARSER_HYPERBOLIC_SINE 20
#define VTK_PARSER_HYPERBOLIC_COSINE 21
#define VTK_PARSER_HYPERBOLIC_TANGENT 22

// functions involving vectors
#define VTK_PARSER_VECTOR_UNARY_MINUS 23
#define VTK_PARSER_DOT_PRODUCT 24
#define VTK_PARSER_VECTOR_ADD 25
#define VTK_PARSER_VECTOR_SUBTRACT 26
#define VTK_PARSER_SCALAR_MULTIPLE 27
#define VTK_PARSER_MAGNITUDE 28
#define VTK_PARSER_NORMALIZE 29

// constants involving vectors
#define VTK_PARSER_IHAT 30
#define VTK_PARSER_JHAT 31
#define VTK_PARSER_KHAT 32

// codes for scalar variables come before those for vectors
#define VTK_PARSER_BEGIN_VARIABLES 33

// the value that is retuned as a result if there is an error
#define VTK_PARSER_ERROR_RESULT VTK_LARGE_FLOAT

class VTK_COMMON_EXPORT vtkFunctionParser : public vtkObject
{
public:
  static vtkFunctionParser *New();
  vtkTypeRevisionMacro(vtkFunctionParser, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Decription:
  // Set/Get input string to evaluate. 
  void SetFunction(const char *function);
  vtkGetStringMacro(Function);

  // Description:
  // Check whether the result is a scalar result.  If it isn't, then
  // either the result is a vector or an error has occurred.
  int IsScalarResult();

  // Description:
  // Check whether the result is a vector result.  If it isn't, then
  // either the result is scalar or an error has occurred.
  int IsVectorResult();

  // Description:
  // Get a scalar result from evaluating the input function.
  double GetScalarResult();

  // Description:
  // Get a vector result from evaluating the input function.
  double* GetVectorResult();
  void GetVectorResult(double result[3]) {
    double *r = this->GetVectorResult();
    result[0] = r[0]; result[1] = r[1]; result[2] = r[2]; };

  // Description:
  // Set the value of a scalar variable.  If a variable with this name
  // exists, then its value will be set to the new value.  If there is not
  // already a variable with this name, variableName will be added to the
  // list of variables, and its value will be set to the new value.
  void SetScalarVariableValue(const char* variableName, double value);
  void SetScalarVariableValue(int i, double value);

  // Description:
  // Get the value of a scalar variable.
  double GetScalarVariableValue(const char* variableName);
  double GetScalarVariableValue(int i);

  // Description:
  // Set the value of a vector variable.  If a variable with this name
  // exists, then its value will be set to the new value.  If there is not
  // already a variable with this name, variableName will be added to the
  // list of variables, and its value will be set to the new value.
  void SetVectorVariableValue(const char* variableName, double xValue,
                              double yValue, double zValue);
  void SetVectorVariableValue(const char* variableName, 
                              const double values[3]) {
    this->SetVectorVariableValue(variableName,values[0],values[1],values[2]);};
  void SetVectorVariableValue(int i, double xValue, double yValue,
                              double zValue);
  void SetVectorVariableValue(int i, const double values[3]) {
    this->SetVectorVariableValue(i,values[0],values[1],values[2]);};
  
  // Description:
  // Get the value of a vector variable.
  double* GetVectorVariableValue(const char* variableName);
  void GetVectorVariableValue(const char* variableName, double value[3]) {
    double *r = this->GetVectorVariableValue(variableName);
    value[0] = r[0]; value[1] = r[1]; value[2] = r[2]; };
  double* GetVectorVariableValue(int i);
  void GetVectorVariableValue(int i, double value[3]) {
    double *r = this->GetVectorVariableValue(i);
    value[0] = r[0]; value[1] = r[1]; value[2] = r[2]; };
  
  // Description:
  // Get the number of scalar variables.
  vtkGetMacro(NumberOfScalarVariables,int);

  // Description:
  // Get the number of vector variables.
  vtkGetMacro(NumberOfVectorVariables,int);

  // Description:
  // Get the ith scalar variable name.
  char* GetScalarVariableName(int i);
  
  // Description:
  // Get the ith vector variable name.
  char* GetVectorVariableName(int i);

  // Description:
  // Remove all the current variables.
  void RemoveAllVariables();
  
protected:
  vtkFunctionParser();
  ~vtkFunctionParser();
  
  int Parse();
  void Evaluate();

  int CheckSyntax();
  void RemoveSpaces();
  char* RemoveSpacesFrom(const char* variableName);
  
  int BuildInternalFunctionStructure();
  void BuildInternalSubstringStructure(int beginIndex, int endIndex);
  void AddInternalByte(unsigned char newByte);
  
  int IsSubstringCompletelyEnclosed(int beginIndex, int endIndex);
  int FindEndOfMathFunction(int beginIndex);
  int FindEndOfMathConstant(int beginIndex);
  
  int IsVariableName(int currentIndex);
  int IsElementaryOperator(int op);
  
  int GetMathFunctionNumber(int currentIndex);
  int GetMathFunctionStringLength(int mathFunctionNumber);
  int GetMathConstantNumber(int currentIndex);
  int GetMathConstantStringLength(int mathConstantNumber);
  int GetElementaryOperatorNumber(char op);
  int GetOperandNumber(int currentIndex);
  int GetVariableNameLength(int variableNumber);
  
  int DisambiguateOperators();
  
  char* Function;
  int FunctionLength;
  int NumberOfScalarVariables;
  int NumberOfVectorVariables;
  char** ScalarVariableNames;
  char** VectorVariableNames;
  double* ScalarVariableValues;
  double** VectorVariableValues;
  unsigned char *ByteCode;
  int ByteCodeSize;
  double *Immediates;
  int ImmediatesSize;
  double *Stack;
  int StackSize;
  int StackPointer;

  vtkTimeStamp FunctionMTime;
  vtkTimeStamp ParseMTime;
  vtkTimeStamp VariableMTime;
  vtkTimeStamp EvaluateMTime;
private:
  vtkFunctionParser(const vtkFunctionParser&);  // Not implemented.
  void operator=(const vtkFunctionParser&);  // Not implemented.
};

#endif
