/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkOrderStatistics.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
  -------------------------------------------------------------------------*/

#include "vtkToolkits.h"

#include "vtkOrderStatistics.h"
#include "vtkStatisticsAlgorithmPrivate.h"

#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"

#include <vtksys/stl/vector>
#include <vtksys/stl/map>
#include <vtksys/stl/set>
#include <vtksys/ios/sstream> 

#define VTK_STATISTICS_NUMBER_OF_VARIABLES 1

vtkCxxRevisionMacro(vtkOrderStatistics, "1.57");
vtkStandardNewMacro(vtkOrderStatistics);

// ----------------------------------------------------------------------
vtkOrderStatistics::vtkOrderStatistics()
{
  this->QuantileDefinition = vtkOrderStatistics::InverseCDFAveragedSteps;
  this->NumberOfIntervals = 4; // By default, calculate 5-points statistics

  this->AssessNames->SetNumberOfValues( 1 );
  this->AssessNames->SetValue( 0, "Quantile" );
}

// ----------------------------------------------------------------------
vtkOrderStatistics::~vtkOrderStatistics()
{
}

// ----------------------------------------------------------------------
void vtkOrderStatistics::PrintSelf( ostream &os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "NumberOfIntervals: " << this->NumberOfIntervals << endl;
  os << indent << "QuantileDefinition: " << this->QuantileDefinition << endl;
}

// ----------------------------------------------------------------------
int vtkOrderStatistics::FillInputPortInformation( int port, vtkInformation* info )
{
  int res; 
  if ( port == INPUT_MODEL )
    {
    info->Set( vtkAlgorithm::INPUT_IS_OPTIONAL(), 1 );
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet" );

    res = 1;
    }
  else
    {
    res = this->Superclass::FillInputPortInformation( port, info );
    }

  return res;
}

// ----------------------------------------------------------------------
int vtkOrderStatistics::FillOutputPortInformation( int port, vtkInformation* info )
{
  int res = this->Superclass::FillOutputPortInformation( port, info );
  if ( port == OUTPUT_MODEL )
    {
    info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    }
  else
    {
    info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkTable" );
    }
  
  return res;
}

// ----------------------------------------------------------------------
void vtkOrderStatistics::SetQuantileDefinition( int qd )
{
  switch ( qd )
    {
    case vtkOrderStatistics::InverseCDF:
      break;
    case vtkOrderStatistics::InverseCDFAveragedSteps:
      break;
    default:
      vtkWarningMacro( "Incorrect type of quantile definition: "
                       <<qd
                       <<". Ignoring it." );
      return;
    }
  
  this->QuantileDefinition =  static_cast<vtkOrderStatistics::QuantileDefinitionType>( qd );
  this->Modified();

  return;
}

// ----------------------------------------------------------------------
bool vtkOrderStatistics::SetParameter( const char* parameter,
                                       int vtkNotUsed( index ),
                                       vtkVariant value )
{
  if ( ! strcmp( parameter, "NumberOfIntervals" ) )
    {
    this->SetNumberOfIntervals( value.ToInt() );

    return true;
    }

  if ( ! strcmp( parameter, "QuantileDefinition" ) )
    {
    this->SetQuantileDefinition( value.ToInt() );

    return true;
    }

  return false;
}

// ----------------------------------------------------------------------
void vtkOrderStatistics::Learn( vtkTable* inData,
                                vtkTable* vtkNotUsed( inParameters ),
                                vtkDataObject* outMetaDO )
{
  vtkMultiBlockDataSet* outMeta = vtkMultiBlockDataSet::SafeDownCast( outMetaDO );
  if ( ! outMeta )
    {
    return;
    }

  // The primary statistics table
  vtkTable* primaryTab = vtkTable::New();

  vtkStringArray* stringCol = vtkStringArray::New();
  stringCol->SetName( "Variable" );
  primaryTab->AddColumn( stringCol );
  stringCol->Delete();

  vtkIdTypeArray* idTypeCol = vtkIdTypeArray::New();
  idTypeCol->SetName( "Cardinality" );
  primaryTab->AddColumn( idTypeCol );
  idTypeCol->Delete();

  if ( ! inData )
    {
    return;
    }

  vtkIdType n = inData->GetNumberOfRows();
  if ( n <= 0 )
    {
    return;
    }

  if ( ! this->Internals->Requests.size() )
    {
    return;
    }

  if ( inData->GetNumberOfColumns() <= 0 )
    {
    return;
    }
  
  vtkVariantArray* variantCol;
  double dq = 1. / static_cast<double>( this->NumberOfIntervals );
  for ( int i = 0; i <= this->NumberOfIntervals; ++ i )
    {
    variantCol = vtkVariantArray::New();
    div_t q = div( i << 2, this->NumberOfIntervals );
    
    if ( q.rem )
      {
      variantCol->SetName( vtkStdString( vtkVariant( i * dq ).ToString() + "-quantile" ).c_str() );
      }
    else
      {
      switch ( q.quot )
        {
        case 0:
          variantCol->SetName( "Minimum" );
          break;
        case 1:
          variantCol->SetName( "First Quartile" );
          break;
        case 2:
          variantCol->SetName( "Median" );
          break;
        case 3:
          variantCol->SetName( "Third Quartile" );
          break;
        case 4:
          variantCol->SetName( "Maximum" );
          break;
        default:
          variantCol->SetName( vtkStdString( vtkVariant( i * dq ).ToString() + "-quantile" ).c_str() );
          break;
        }
      }
    primaryTab->AddColumn( variantCol );
    variantCol->Delete();
    }

  // Loop over requests
  for ( vtksys_stl::set<vtksys_stl::set<vtkStdString> >::iterator rit = this->Internals->Requests.begin(); 
        rit != this->Internals->Requests.end(); ++ rit )
    {
    // Each request contains only one column of interest (if there are others, they are ignored)
    vtksys_stl::set<vtkStdString>::const_iterator it = rit->begin();
    vtkStdString col = *it;
    if ( ! inData->GetColumnByName( col ) )
      {
      vtkWarningMacro( "InData table does not have a column "
                       << col.c_str()
                       << ". Ignoring it." );
      continue;
      }

    vtkVariantArray* row = vtkVariantArray::New();

    // A row contains: variable name, cardinality, and NumberOfIntervals + 1 quantile values
    row->SetNumberOfValues( this->NumberOfIntervals + 3 );
    
    int i = 0;
    row->SetValue( i ++, col );
    row->SetValue( i ++, n );
    
    vtksys_stl::vector<double> quantileThresholds;
    double dh = n / static_cast<double>( this->NumberOfIntervals );
    for ( int j = 0; j < this->NumberOfIntervals; ++ j )
      {
      quantileThresholds.push_back( j * dh );
      }

    // Try to downcast column to eith data or string arrays for efficient data access
    vtkAbstractArray* arr = inData->GetColumnByName( col );

    if ( arr->IsA("vtkDataArray") ) 
      {
      vtkDataArray* darr = vtkDataArray::SafeDownCast( arr );

      vtksys_stl::map<double,vtkIdType> distr;
      for ( vtkIdType r = 0; r < n; ++ r )
        {
        ++ distr[darr->GetTuple1( r )];
        }

      vtkIdType sum = 0;
      vtksys_stl::vector<double>::iterator qit = quantileThresholds.begin();
      for ( vtksys_stl::map<double,vtkIdType>::iterator mit = distr.begin();
            mit != distr.end(); ++ mit  )
        {
        for ( sum += mit->second; qit != quantileThresholds.end() && sum >= *qit; ++ qit )
          {
          // Mid-point interpolation is available for numeric types only
          if ( sum == *qit
               && this->QuantileDefinition == vtkOrderStatistics::InverseCDFAveragedSteps )
            {
            vtksys_stl::map<double,vtkIdType>::iterator nit = mit;
            row->SetValue( i ++, ( (++ nit)->first + mit->first ) * .5 );
            }
          else
            {
            row->SetValue( i ++, mit->first );
            }
          }
        }
    
      row->SetValue( i, distr.rbegin()->first );
      primaryTab->InsertNextRow( row );
      }
    else if ( arr->IsA("vtkStringArray") ) 
      {
      vtkStringArray* sarr = vtkStringArray::SafeDownCast( arr ); 
      vtksys_stl::map<vtkStdString,vtkIdType> distr;
      for ( vtkIdType r = 0; r < n; ++ r )
        {
        ++ distr[sarr->GetValue( r )];
        }

      vtkIdType sum = 0;
      vtksys_stl::vector<double>::iterator qit = quantileThresholds.begin();
      for ( vtksys_stl::map<vtkStdString,vtkIdType>::iterator mit = distr.begin();
            mit != distr.end(); ++ mit  )
        {
        for ( sum += mit->second; qit != quantileThresholds.end() && sum >= *qit; ++ qit )
          {
          row->SetValue( i ++, mit->first );
          }
        }
      
      row->SetValue( i, distr.rbegin()->first );
      primaryTab->InsertNextRow( row );
      }
    else // column is of type vtkVariantArray, which is not supported by this filter
      {
      vtkWarningMacro("Type vtkVariantArray of column "
                      <<col.c_str()
                      << " not supported. Ignoring it." );
      }

    row->Delete();
    }

  // Finally set first block of output meta port to primary statistics table
  outMeta->SetNumberOfBlocks( 1 );
  outMeta->GetMetaData( static_cast<unsigned>( 0 ) )->Set( vtkCompositeDataSet::NAME(), "Primary Statistics" );
  outMeta->SetBlock( 0, primaryTab );

  // Clean up
  primaryTab->Delete();

  return;
}

// ----------------------------------------------------------------------
void vtkOrderStatistics::Derive( vtkDataObject* vtkNotUsed( inMeta ) )
{
}

// ----------------------------------------------------------------------
class TableColumnBucketingFunctor : public vtkStatisticsAlgorithm::AssessFunctor
{
public:
  vtkAbstractArray* Data;
  vtkVariantArray* Quantiles;

  TableColumnBucketingFunctor( vtkAbstractArray* vals,
                               vtkVariantArray* quantiles )
  {
    this->Data = vals;
    this->Quantiles = quantiles;
  }
  virtual ~TableColumnBucketingFunctor() { }
  virtual void operator() ( vtkVariantArray* result,
                            vtkIdType id )
  {
    vtkVariant x = this->Data->GetVariantValue( id );

    if ( x < this->Quantiles->GetValue( 2 ) ) // Value #0 is the variable name and #1 is the cardinality
      {
      // x is smaller than lower bound
      result->SetNumberOfValues( 1 );
      result->SetValue( 0, 0 );

      return;
      }

    vtkIdType n = this->Quantiles->GetNumberOfValues() + 2;
    vtkIdType q = 3;
    while ( q < n && x > this->Quantiles->GetValue( q ) )
      {
      ++ q;
      }

    result->SetNumberOfValues( 1 );
    result->SetValue( 0, q - 2 ); 
  }
};

// ----------------------------------------------------------------------
void vtkOrderStatistics::Assess( vtkTable* inData,
                                 vtkDataObject* inMetaDO,
                                 vtkTable* outData )
{
  if ( ! inData || inData->GetNumberOfColumns() <= 0 )
    {
    return;
    }

  vtkIdType nRowData = inData->GetNumberOfRows();
  if ( nRowData <= 0 )
    {
    return;
    }

  vtkMultiBlockDataSet* inMeta = vtkMultiBlockDataSet::SafeDownCast( inMetaDO );
  if ( ! inMeta || inMeta->GetNumberOfBlocks() < 2 )
    {
    return;
    }

  // Loop over requests
  for ( vtksys_stl::set<vtksys_stl::set<vtkStdString> >::const_iterator rit = this->Internals->Requests.begin(); 
        rit != this->Internals->Requests.end(); ++ rit )
    {
    // Each request contains only one column of interest (if there are others, they are ignored)
    vtksys_stl::set<vtkStdString>::const_iterator it = rit->begin();
    vtkStdString varName = *it;
    if ( ! inData->GetColumnByName( varName ) )
      {
      vtkWarningMacro( "InData table does not have a column "
                       << varName.c_str()
                       << ". Ignoring it." );
      continue;
      }

    vtkStringArray* varNames = vtkStringArray::New();
    varNames->SetNumberOfValues( VTK_STATISTICS_NUMBER_OF_VARIABLES );
    varNames->SetValue( 0, varName );

    // Store names to be able to use SetValueByName, and create the outData columns
    int nv = this->AssessNames->GetNumberOfValues();
    vtkStdString* names = new vtkStdString[nv];
    for ( int v = 0; v < nv; ++ v )
      {
      vtksys_ios::ostringstream assessColName;
      assessColName << this->AssessNames->GetValue( v )
                    << "("
                    << varName
                    << ")";

      names[v] = assessColName.str().c_str(); 

      vtkDoubleArray* assessValues = vtkDoubleArray::New(); 
      assessValues->SetName( names[v] ); 
      assessValues->SetNumberOfTuples( nRowData  ); 
      outData->AddColumn( assessValues ); 
      assessValues->Delete(); 
      }

    // Select assess functor
    AssessFunctor* dfunc;
    this->SelectAssessFunctor( outData,
                               inMeta,
                               varNames,
                               dfunc );

    if ( ! dfunc )
      {
      // Functor selection did not work. Do nothing.
      vtkWarningMacro( "AssessFunctors could not be allocated for column "
                       << varName.c_str()
                       << ". Ignoring it." );
      }
    else
      {
      // Assess each entry of the column
      vtkVariantArray* assessResult = vtkVariantArray::New();
      for ( vtkIdType r = 0; r < nRowData; ++ r )
        {
        (*dfunc)( assessResult, r );
        for ( int v = 0; v < nv; ++ v )
          {
          outData->SetValueByName( r, names[v], assessResult->GetValue( v ) );
          }
        }

      assessResult->Delete();
      }

    delete dfunc;
    delete [] names;
    varNames->Delete(); // Do not delete earlier! Otherwise, dfunc will be wrecked
    }
}

// ----------------------------------------------------------------------
void vtkOrderStatistics::SelectAssessFunctor( vtkTable* outData,
                                              vtkDataObject* inMetaDO,
                                              vtkStringArray* rowNames,
                                              AssessFunctor*& dfunc )
{
  vtkTable* inMeta = vtkTable::SafeDownCast( inMetaDO ); 
  if ( ! inMeta ) 
    { 
    return; 
    }

  vtkStdString varName = rowNames->GetValue( 0 );

  // Downcast meta columns to string arrays for efficient data access
  vtkStringArray* vars = vtkStringArray::SafeDownCast( inMeta->GetColumnByName( "Variable" ) );
  if ( ! vars )
    {
    dfunc = 0;
    return;
    }

  // Loop over parameters table until the requested variable is found
  vtkIdType nRowP = inMeta->GetNumberOfRows();
  for ( int r = 0; r < nRowP; ++ r )
    {
    if ( vars->GetValue( r ) == varName )
      {
      // Grab the data for the requested variable
      vtkAbstractArray* vals = outData->GetColumnByName( varName );
      if ( ! vals )
        {
        dfunc = 0;
        return;
        }

      dfunc = new TableColumnBucketingFunctor( vals, inMeta->GetRow( r ) );
      return;
      }
    }

  // The variable of interest was not found in the parameter table
  dfunc = 0;
}
