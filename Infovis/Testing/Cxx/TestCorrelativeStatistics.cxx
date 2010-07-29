/*
 * Copyright 2008 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
// .SECTION Thanks
// Thanks to Philippe Pebay from Sandia National Laboratories 
// for implementing this test.

#include "vtkDataObjectCollection.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkStringArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkTable.h"
#include "vtkTimerLog.h"
#include "vtkCorrelativeStatistics.h"

//=============================================================================
int TestCorrelativeStatistics( int, char *[] )
{
  int testStatus = 0;

  double mingledData[] = 
    {
      46,
      45,
      47,
      49,
      46,
      47,
      46,
      46,
      47,
      46,
      47,
      49,
      49,
      49,
      47,
      45,
      50,
      50,
      46,
      46,
      51,
      50,
      48,
      48,
      52,
      54,
      48,
      47,
      52,
      52,
      49,
      49,
      53,
      54,
      50,
      50,
      53,
      54,
      50,
      52,
      53,
      53,
      50,
      51,
      54,
      54,
      49,
      49,
      52,
      52,
      50,
      51,
      52,
      52,
      49,
      47,
      48,
      48,
      48,
      50,
      46,
      48,
      47,
      47,
    };
  int nVals1 = 32;

  vtkDoubleArray* dataset1Arr = vtkDoubleArray::New();
  dataset1Arr->SetNumberOfComponents( 1 );
  dataset1Arr->SetName( "Metric 0" );

  vtkDoubleArray* dataset2Arr = vtkDoubleArray::New();
  dataset2Arr->SetNumberOfComponents( 1 );
  dataset2Arr->SetName( "Metric 1" );

  vtkDoubleArray* dataset3Arr = vtkDoubleArray::New();
  dataset3Arr->SetNumberOfComponents( 1 );
  dataset3Arr->SetName( "Metric 2" );

  for ( int i = 0; i < nVals1; ++ i )
    {
    int ti = i << 1;
    dataset1Arr->InsertNextValue( mingledData[ti] );
    dataset2Arr->InsertNextValue( mingledData[ti + 1] );
    dataset3Arr->InsertNextValue( -1. );
    }

  vtkTable* datasetTable1 = vtkTable::New();
  datasetTable1->AddColumn( dataset1Arr );
  dataset1Arr->Delete();
  datasetTable1->AddColumn( dataset2Arr );
  dataset2Arr->Delete();
  datasetTable1->AddColumn( dataset3Arr );
  dataset3Arr->Delete();

  // Pairs of interest
  int nMetricPairs = 2;
  vtkStdString columnPairs[] = 
    { 
      "Metric 0", "Metric 1", // First pair
      "Metric 2", "Metric 1"  // Second pair
    };

  // Reference values
  // Means and variances for metrics 0, and 1, respectively
  double meansX1[] = { 49.21875, 49.5 };
  double varsX1[] = { 5.9828629, 7.548397 };

  // Means and variances for metrics 1 and 2, respectively
  double meansY1[] = { 49.5, -1. };
  double varsY1[] = { 7.548397, 0. };

  // Covariance matrix of (metric 0, metric 1) pair
  double covariance1[] = { 5.98286, 7.54839, 6.14516 }; 

  // Pearson r for each of the three pairs
  double correlations1[] = { 0.914433, 0. }; 

  // Threshold for outlier detection
  double threshold = 4.;

  // Set correlative statistics algorithm and its input data port
  vtkCorrelativeStatistics* cs1 = vtkCorrelativeStatistics::New();

  // First verify that absence of input does not cause trouble
  cout << "## Verifying that absence of input does not cause trouble... ";
  cs1->Update();
  cout << "done.\n";

  // Prepare first test with data
  cs1->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable1 );
  datasetTable1->Delete();

  // Select Column Pairs of Interest ( Learn Mode ) 
  // 1.1: a valid pair
  cs1->AddColumnPair( "Metric 0", "Metric 1" ); 
  // 1.2: the same valid pair, just reversed -- should thus be ignored
  cs1->AddColumnPair( "Metric 1", "Metric 0" );
  // 2: another valid pair
  cs1->AddColumnPair( "Metric 2", "Metric 1" ); 
  // 3: an invalid pair
  cs1->AddColumnPair( "Metric 1", "Metric 3" ); 

  // Test Learn, Derive, Test, and Assess options
  cs1->SetLearnOption( true );
  cs1->SetDeriveOption( true );
  cs1->SetAssessOption( true );
  cs1->SetTestOption( true );
  cs1->Update();

  // Get output data and meta tables
  vtkTable* outputData1 = cs1->GetOutput( vtkStatisticsAlgorithm::OUTPUT_DATA );
  vtkMultiBlockDataSet* outputMetaDS1 = vtkMultiBlockDataSet::SafeDownCast( cs1->GetOutputDataObject( vtkStatisticsAlgorithm::OUTPUT_MODEL ) );
  vtkTable* outputPrimary1 = vtkTable::SafeDownCast( outputMetaDS1->GetBlock( 0 ) );
  vtkTable* outputDerived1 = vtkTable::SafeDownCast( outputMetaDS1->GetBlock( 1 ) );
  vtkTable* outputTest1 = cs1->GetOutput( vtkStatisticsAlgorithm::OUTPUT_TEST );

  cout << "## Calculated the following primary statistics for first data set:\n";
  for ( vtkIdType r = 0; r < outputPrimary1->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputPrimary1->GetNumberOfColumns(); ++ i )
      {
      cout << outputPrimary1->GetColumnName( i )
           << "="
           << outputPrimary1->GetValue( r, i ).ToString()
           << "  ";
      }

    // Verify some of the calculated primary statistics
    if ( fabs ( outputPrimary1->GetValueByName( r, "Mean X" ).ToDouble() - meansX1[r] ) > 1.e-6 )
      {
      vtkGenericWarningMacro("Incorrect mean for X");
      testStatus = 1;
      }

    if ( fabs ( outputPrimary1->GetValueByName( r, "Mean Y" ).ToDouble() - meansY1[r] ) > 1.e-6 )
      {
      vtkGenericWarningMacro("Incorrect mean for Y");
      testStatus = 1;
      }
    cout << "\n";
    }

  cout << "\n## Calculated the following derived statistics for first data set:\n";
  for ( vtkIdType r = 0; r < outputDerived1->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputDerived1->GetNumberOfColumns(); ++ i )
      {
      cout << outputDerived1->GetColumnName( i )
           << "="
           << outputDerived1->GetValue( r, i ).ToString()
           << "  ";
      }

    // Verify some of the calculated derived statistics
    if ( fabs ( outputDerived1->GetValueByName( r, "Variance X" ).ToDouble() - varsX1[r] ) > 1.e-5 )
      {
      vtkGenericWarningMacro("Incorrect variance for X");
      testStatus = 1;
      }

    if ( fabs ( outputDerived1->GetValueByName( r, "Variance Y" ).ToDouble() - varsY1[r] ) > 1.e-5 )
      {
      vtkGenericWarningMacro("Incorrect variance for Y");
      testStatus = 1;
      }

    if ( fabs ( outputDerived1->GetValueByName( r, "Pearson r" ).ToDouble() - correlations1[r] ) > 1.e-6 )
      {
      vtkGenericWarningMacro("Incorrect correlation coefficient");
      testStatus = 1;
      }
    cout << "\n";
    }

  // Check some results of the Test option
  cout << "\n## Calculated the following Jarque-Bera-Srivastava statistics:\n";
  for ( vtkIdType r = 0; r < outputTest1->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputTest1->GetNumberOfColumns(); ++ i )
      {
      cout << outputTest1->GetColumnName( i )
           << "="
           << outputTest1->GetValue( r, i ).ToString()
           << "  ";
      }

    cout << "\n";
    }

  // Select Column Pairs of Interest ( Assess Mode ) 
  cs1->ResetRequests(); // Clear existing pairs
  cs1->AddColumnPair( columnPairs[0], columnPairs[1] ); // A valid pair

  cout << "\n## Searching for outliers with respect to this bivariate Gaussian distribution:\n"
       << "   (X, Y) = ("
       << columnPairs[0]
       << ", "
       << columnPairs[1]
       << "), mean=("
       << meansX1[0]
       << ", "
       << meansY1[0]
       << "), covariance=["
       << covariance1[0]
       << ", "
       << covariance1[2]
       << " ; "
       << covariance1[2]
       << ", "
       << covariance1[1]
       << "], Squared Mahalanobis > "
       << threshold
       << "\n";

  int nOutliers = 0;
  int tableIdx[] = { 0, 1, 3 };
  cout << "   Found the following outliers:\n";
  for ( int i = 0; i < 3; ++ i )
    {
    cout << "   "
         << outputData1->GetColumnName( tableIdx[i] );
    }
  cout << "\n";

  for ( vtkIdType r = 0; r < outputData1->GetNumberOfRows(); ++ r )
    {
    if ( outputData1->GetValue( r, tableIdx[2] ).ToDouble() > threshold )
      {
      ++ nOutliers;

      for ( int i = 0; i < 3; ++ i )
        {
        cout << "     "
             << outputData1->GetValue( r,  tableIdx[i] ).ToDouble()
             << "    ";
        }
      cout << "\n";
      }
    }

  if ( nOutliers != 3 )
    {
    vtkGenericWarningMacro("Expected 3 outliers, found " << nOutliers << ".");
    testStatus = 1;
    }

  // Test with a slight variation of initial data set (to test model aggregation)
  int nVals2 = 32;

  vtkDoubleArray* dataset4Arr = vtkDoubleArray::New();
  dataset4Arr->SetNumberOfComponents( 1 );
  dataset4Arr->SetName( "Metric 0" );

  vtkDoubleArray* dataset5Arr = vtkDoubleArray::New();
  dataset5Arr->SetNumberOfComponents( 1 );
  dataset5Arr->SetName( "Metric 1" );

  vtkDoubleArray* dataset6Arr = vtkDoubleArray::New();
  dataset6Arr->SetNumberOfComponents( 1 );
  dataset6Arr->SetName( "Metric 2" );

  for ( int i = 0; i < nVals2; ++ i )
    {
    int ti = i << 1;
    dataset4Arr->InsertNextValue( mingledData[ti] + 1. );
    dataset5Arr->InsertNextValue( mingledData[ti + 1] );
    dataset6Arr->InsertNextValue( 1. );
    }

  vtkTable* datasetTable2 = vtkTable::New();
  datasetTable2->AddColumn( dataset4Arr );
  dataset4Arr->Delete();
  datasetTable2->AddColumn( dataset5Arr );
  dataset5Arr->Delete();
  datasetTable2->AddColumn( dataset6Arr );
  dataset6Arr->Delete();

  // Set correlative statistics algorithm and its input data port
  vtkCorrelativeStatistics* cs2 = vtkCorrelativeStatistics::New();
  cs2->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable2 );
  datasetTable2->Delete();

  // Select all column pairs as pairs of interest
  for ( int i = 0; i< nMetricPairs; ++ i )
    {  // Add all valid pairs
    cs2->AddColumnPair( columnPairs[2 * i], columnPairs[ 2 * i + 1] );
    }

  // Update with Learn option only
  cs2->SetLearnOption( true );
  cs2->SetDeriveOption( false );
  cs2->SetAssessOption( false );
  cs2->Update();

  // Get output meta tables
  vtkMultiBlockDataSet* outputMetaDS2 = vtkMultiBlockDataSet::SafeDownCast( cs2->GetOutputDataObject( vtkStatisticsAlgorithm::OUTPUT_MODEL ) );
  vtkTable* outputPrimary2 = vtkTable::SafeDownCast( outputMetaDS2->GetBlock( 0 ) );

  cout << "\n## Calculated the following primary statistics for second data set:\n";
  for ( vtkIdType r = 0; r < outputPrimary2->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputPrimary2->GetNumberOfColumns(); ++ i )
      {
      cout << outputPrimary2->GetColumnName( i )
           << "="
           << outputPrimary2->GetValue( r, i ).ToString()
           << "  ";
      }
    cout << "\n";
    }

  // Now build a data object collection of the two obtained models
  vtkDataObjectCollection* doc = vtkDataObjectCollection::New();
  doc->AddItem( outputMetaDS1 );
  doc->AddItem( outputMetaDS2 );

  // And calculate the aggregated minimal statistics of the two models
  vtkCorrelativeStatistics* cs0 = vtkCorrelativeStatistics::New();
  vtkMultiBlockDataSet* aggregated = vtkMultiBlockDataSet::New();
  cs0->Aggregate( doc, aggregated );

  // Finally, calculate the derived statistics of the aggregated model
  cs0->SetInput( vtkStatisticsAlgorithm::INPUT_MODEL, aggregated );
  cs0->SetLearnOption( false );
  cs0->SetDeriveOption( true ); 
  cs0->SetAssessOption( false );
  cs0->Update();

  // Reference values
  // Means and variances for metrics 0 and 1, respectively
  double meansX0[] = { 49.71875 , 49.5 };
  double varsX0[] = { 6.1418651 , 7.548397 * 62. / 63. };

  // Means and variances for metrics 1 and 2, respectively
  double meansY0[] = { 49.5, 0. };
  double varsY0[] = { 7.548397 * 62. / 63., 64. / 63. };

  // Pearson r for each of the three pairs
  double correlations0[] = { 0.895327, 0. };

  // Get output meta tables
  vtkMultiBlockDataSet* outputMetaDS0 = vtkMultiBlockDataSet::SafeDownCast( cs0->GetOutputDataObject( vtkStatisticsAlgorithm::OUTPUT_MODEL ) );
  vtkTable* outputPrimary0 = vtkTable::SafeDownCast( outputMetaDS0->GetBlock( 0 ) );
  vtkTable* outputDerived0 = vtkTable::SafeDownCast( outputMetaDS0->GetBlock( 1 ) );

  cout << "\n## Calculated the following primary statistics for aggregated (first + second) data set:\n";
  for ( vtkIdType r = 0; r < outputPrimary0->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputPrimary0->GetNumberOfColumns(); ++ i )
      {
      cout << outputPrimary0->GetColumnName( i )
           << "="
           << outputPrimary0->GetValue( r, i ).ToString()
           << "  ";
      }

    // Verify some of the calculated primary statistics
    if ( fabs ( outputPrimary0->GetValueByName( r, "Mean X" ).ToDouble() - meansX0[r] ) > 1.e-6 )
      {
      vtkGenericWarningMacro("Incorrect mean for X");
      testStatus = 1;
      }

    if ( fabs ( outputPrimary0->GetValueByName( r, "Mean Y" ).ToDouble() - meansY0[r] ) > 1.e-6 )
      {
      vtkGenericWarningMacro("Incorrect mean for Y");
      testStatus = 1;
      }
    cout << "\n";
    }

  cout << "\n## Calculated the following derived statistics for aggregated (first + second) data set:\n";
  for ( vtkIdType r = 0; r < outputDerived0->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputDerived0->GetNumberOfColumns(); ++ i )
      {
      cout << outputDerived0->GetColumnName( i )
           << "="
           << outputDerived0->GetValue( r, i ).ToString()
           << "  ";
      }

    // Verify some of the calculated derived statistics
    if ( fabs ( outputDerived0->GetValueByName( r, "Variance X" ).ToDouble() - varsX0[r] ) > 1.e-5 )
      {
      vtkGenericWarningMacro("Incorrect variance for X");
      testStatus = 1;
      }

    if ( fabs ( outputDerived0->GetValueByName( r, "Variance Y" ).ToDouble() - varsY0[r] ) > 1.e-5 )
      {
      vtkGenericWarningMacro("Incorrect variance for Y");
      testStatus = 1;
      }

    if ( fabs ( outputDerived0->GetValueByName( r, "Pearson r" ).ToDouble() - correlations0[r] ) > 1.e-6 )
      {
      vtkGenericWarningMacro("Incorrect correlation coefficient");
      testStatus = 1;
      }
    cout << "\n";
    }

  // Clean up
  cs0->Delete();
  cs1->Delete();
  cs2->Delete();
  doc->Delete();
  aggregated->Delete();

  // ************** Pseudo-random sample to exercise Jarque-Bera-Srivastava test *********
  int nVals = 10000;

  vtkDoubleArray* datasetBinormalX = vtkDoubleArray::New();
  datasetBinormalX->SetNumberOfComponents( 1 );
  datasetBinormalX->SetName( "Standard Binormal X" );

  vtkDoubleArray* datasetBinormalY = vtkDoubleArray::New();
  datasetBinormalY->SetNumberOfComponents( 1 );
  datasetBinormalY->SetName( "Standard Binormal Y" );

  vtkDoubleArray* datasetUniform = vtkDoubleArray::New();
  datasetUniform->SetNumberOfComponents( 1 );
  datasetUniform->SetName( "Standard Uniform" );

  vtkDoubleArray* datasetLaplace = vtkDoubleArray::New();
  datasetLaplace->SetNumberOfComponents( 1 );
  datasetLaplace->SetName( "Standard Laplace" );

  // Seed random number generator
  vtkMath::RandomSeed( static_cast<int>( vtkTimerLog::GetUniversalTime() ) );

  // Pre-set Pearson correlation coefficient
  double rho = .8;
  double ror = sqrt( 1. - rho * rho );
  double x, y;
  for ( int i = 0; i < nVals; ++ i )
    {
    x = vtkMath::Gaussian();
    y = rho * x + ror * vtkMath::Gaussian();
    datasetBinormalX->InsertNextValue( x );
    datasetBinormalY->InsertNextValue( y );
    datasetUniform->InsertNextValue( vtkMath::Random() );
    double u = vtkMath::Random() - .5;
    datasetLaplace->InsertNextValue( ( u < 0. ? 1. : -1. ) * log ( 1. - 2. * fabs( u ) ) );
    }

  vtkTable* testTable = vtkTable::New();
  testTable->AddColumn( datasetBinormalX );
  datasetBinormalX->Delete();
  testTable->AddColumn( datasetBinormalY );
  datasetBinormalY->Delete();
  testTable->AddColumn( datasetUniform );
  datasetUniform->Delete();
  testTable->AddColumn( datasetLaplace );
  datasetLaplace->Delete();

  // Set descriptive statistics algorithm and its input data port
  vtkCorrelativeStatistics* cs4 = vtkCorrelativeStatistics::New();
  cs4->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, testTable );
  testTable->Delete();

  // Select Column Pairs of Interest ( Learn Mode )
  cs4->AddColumnPair( "Standard Binormal X", "Standard Binormal Y" );
  cs4->AddColumnPair( "Standard Binormal X", "Standard Uniform" );
  cs4->AddColumnPair( "Standard Laplace", "Standard Binormal Y" );
  cs4->AddColumnPair( "Standard Uniform", "Standard Laplace" );

  // Test Learn, Derive, and Test options only
  cs4->SetLearnOption( true );
  cs4->SetDeriveOption( true );
  cs4->SetTestOption( true );
  cs4->SetAssessOption( false );
  cs4->Update();

  // Get output data and meta tables
  vtkMultiBlockDataSet* outputMetaCS4 = vtkMultiBlockDataSet::SafeDownCast( cs4->GetOutputDataObject( vtkStatisticsAlgorithm::OUTPUT_MODEL ) );
  vtkTable* outputPrimary4 = vtkTable::SafeDownCast( outputMetaCS4->GetBlock( 0 ) );
  vtkTable* outputDerived4 = vtkTable::SafeDownCast( outputMetaCS4->GetBlock( 1 ) );
  vtkTable* outputTest4 = cs4->GetOutput( vtkStatisticsAlgorithm::OUTPUT_TEST );

  cout << "\n## Calculated the following primary statistics for pseudo-random variables (n="
       << nVals
       << "):\n";
  for ( vtkIdType r = 0; r < outputPrimary4->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputPrimary4->GetNumberOfColumns(); ++ i )
      {
      cout << outputPrimary4->GetColumnName( i )
           << "="
           << outputPrimary4->GetValue( r, i ).ToString()
           << "  ";
      }

    cout << "\n";
    }

  cout << "\n## Calculated the following derived statistics for pseudo-random variables (n="
       << nVals
       << "):\n";
  for ( vtkIdType r = 0; r < outputDerived4->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int i = 0; i < outputDerived4->GetNumberOfColumns(); ++ i )
      {
      cout << outputDerived4->GetColumnName( i )
           << "="
           << outputDerived4->GetValue( r, i ).ToString()
           << "  ";
      }

    cout << "\n";
    }

  // Check some results of the Test option
  cout << "\n## Calculated the following Jarque-Bera-Srivastava statistics for pseudo-random variables (n="
       << nVals
       << "):\n";

#ifdef VTK_USE_GNU_R
  int nNonGaussian = 3;
  int nRejected = 0;
  double alpha = .01;
#endif // VTK_USE_GNU_R

  // Loop over Test table
  for ( vtkIdType r = 0; r < outputTest4->GetNumberOfRows(); ++ r )
    {
    cout << "   ";
    for ( int c = 0; c < outputTest4->GetNumberOfColumns(); ++ c )
      {
      cout << outputTest4->GetColumnName( c )
           << "="
           << outputTest4->GetValue( r, c ).ToString()
           << "  ";
      }

#ifdef VTK_USE_GNU_R
    // Check if null hypothesis is rejected at specified significance level
    double p = outputTest4->GetValueByName( r, "P" ).ToDouble();
    // Must verify that p value is valid (it is set to -1 if R has failed)
    if ( p > -1 && p < alpha )
      {
      cout << "Null hypothesis (normality) rejected at "
           << alpha
           << " significance level";

      ++ nRejected;
      }
#endif // VTK_USE_GNU_R

    cout << "\n";
    }

#ifdef VTK_USE_GNU_R
  if ( nRejected < nNonGaussian )
    {
    vtkGenericWarningMacro("Rejected only "
                           << nRejected
                           << " null hypotheses of normality whereas "
                           << nNonGaussian
                           << " variables are not Gaussian");
    testStatus = 1;
    }
#endif // VTK_USE_GNU_R

  // Clean up
  cs4->Delete();

  return testStatus;
}
