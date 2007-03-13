/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTemporalStreamTracer.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTemporalStreamTracer - A Parallel Particle tracer for unsteady vector fields
// .SECTION Description
// vtkTemporalStreamTracer is a filter that integrates a vector field to generate
// 
//
// .SECTION See Also
// vtkRibbonFilter vtkRuledSurfaceFilter vtkInitialValueProblemSolver 
// vtkRungeKutta2 vtkRungeKutta4 vtkRungeKutta45 vtkStreamTracer

#ifndef __vtkTemporalStreamTracer_h
#define __vtkTemporalStreamTracer_h

#include "vtkToolkits.h" // For VTK_USE_MPI 
#include "vtkSmartPointer.h"
#include "vtkStreamTracer.h"

//BTX
#include <vtkstd/vector> // I'll remove these soon.
#include <vtkstd/list>
//ETX
//#define JB_H5PART_PARTICLE_OUTPUT 1

#ifdef VTK_USE_MPI
  class vtkMultiProcessController;
#endif

class vtkCompositeDataSet;
class vtkDataArray;
class vtkDoubleArray;
class vtkGenericCell;
class vtkIntArray;
class vtkTemporalInterpolatedVelocityField;
class vtkPoints;
class vtkCellArray;
class vtkDoubleArray;
class vtkFloatArray;
class vtkIntArray;
#ifdef JB_H5PART_PARTICLE_OUTPUT
  class vtkH5PartWriter;
#endif

//BTX
class vtkTemporalStreamTracerInternals : public vtkObjectBase {
  public:
    static vtkTemporalStreamTracerInternals *New() { return new vtkTemporalStreamTracerInternals(); }
    //
};
//ETX

//BTX
    typedef struct { double x[4]; } Position;
    typedef struct {
      // These are used during iteration
      int           Counter;
      int           Index;
      bool          Wrap;
      Position      CurrentPosition;
      int           CachedDataSet[2];
      vtkIdType     CachedCellId[2];
      // These are computed scalars we might display
      int           SourceID;
      int           InjectedPointId;      
      float         UniqueParticleId;
      float         vorticity;
      float         rotation;
      float         angularVel;
    } ParticleInformation;

    struct ParticleLifetime {
      ParticleInformation      Information;
      vtkstd::vector<Position> Coordinates;
    };

    typedef vtkstd::vector<ParticleInformation>  ParticleList;
    typedef vtkstd::list<ParticleLifetime>       ParticleDataList;
    typedef ParticleDataList::iterator           ParticleIterator;
//ETX

class VTK_PARALLEL_EXPORT vtkTemporalStreamTracer : public vtkStreamTracer
{
public:

    vtkTypeRevisionMacro(vtkTemporalStreamTracer,vtkStreamTracer);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Construct object to start from position (0,0,0), integrate forward,
    // terminal speed 1.0E-12, vorticity computation on, integration
    // step length 0.5 (unit cell length), maximum number of steps 2000,
    // using 2nd order Runge Kutta and maximum propagation 1.0 (unit length).
    static vtkTemporalStreamTracer *New();

    vtkSetMacro(TimeStep,unsigned int);
    vtkGetMacro(TimeStep,unsigned int);

    // Description:
    // If the Particles being traced animate for a long time, the
    // trails or traces will become long and stringy. Setting
    // the MaxTraceTimeLength will limit how much of the trace
    // is displayed. Tracks longer then the Max will disappear
    // and the trace will apppear like a snake of fixed length
    // which progresses as the particle moves
    vtkSetMacro(MaxTrackLength,unsigned int);
    vtkGetMacro(MaxTrackLength,unsigned int);

    // Description:
    // If storing multiple steps along the track, use TrackMarkSpaceRatio
    // to determine which steps to save. Use a value >1. 
    // TrackMarkSpaceRatio=1 and MaxTrackLength=2 wil always
    // store the last two points. Using TrackMarkSpaceRatio=2 
    // and MaxTrackLength=4 will store the 1,3,5,7 last points
    // et cetera.
    // NB. When crossing domain boundaries from one processor to another
    // the stored track lists may be corrupted(=reset). This is a ToDo!
    vtkSetMacro(TrackMarkSpaceRatio,unsigned int);
    vtkGetMacro(TrackMarkSpaceRatio,unsigned int);

    // Description:
    // If the data source does not have the correct time values 
    // present on each time step - setting this value to non unity can
    // be used to adjust the time step size from 1s pre step to
    // 1xTimeStepResolution
    vtkSetMacro(TimeStepResolution,double);
    vtkGetMacro(TimeStepResolution,double);

    // Description:
    // When animating particles, it is nice to inject new ones every Nth step
    // to produce a continuous flow. Setting ForceReinjectionEveryNSteps to a 
    // non zero value will cause the particle source to reinject particles
    // every Nth step even if it is otherwise unchanged.
    // Note that if the particle source is also animated, this flag will be
    // redundant as the particles will be reinjected whenever the source changes 
    // anyway
    vtkSetMacro(ForceReinjectionEveryNSteps,int);
    vtkGetMacro(ForceReinjectionEveryNSteps,int);

    // Description:
    // Internal Memory allocation for particles will grow by
    // this amount when new memory is needed. If large numbers of particles
    // are likely to be generated during N timesteps, then This number should
    // be set quite large. The default is 1000, but values of 10000 or 100000 
    // may be a good idea for very large simulation/animations
//    vtkSetMacro(MemoryAllocationSize,int);
//    vtkGetMacro(MemoryAllocationSize,int);

    // Description:
    // Specify an alternative Geometry object as the source of particles
    // This method exists so that in the ParaView GUI we can either
    // select a widget as source, or by using this input, a dataset generated
    // during the session.
    // Old style. Do not use.
    void SetSource2(vtkDataSet *source);
    vtkDataSet *GetSource2();

    // Description:
    // Specify an alternative Geometry object as the source of particles
    // This method exists so that in the ParaView GUI we can either
    // select a widget as source, or by using this input, a dataset generated
    // during the session.
    // New Style. Do use.
    void SetSource2Connection(vtkAlgorithmOutput* algOutput);

    // Description:
    vtkSetMacro(EnableSource1,int);
    vtkGetMacro(EnableSource1,int);
    vtkBooleanMacro(EnableSource1,int);

    // Description:
    vtkSetMacro(EnableSource2,int);
    vtkGetMacro(EnableSource2,int);
    vtkBooleanMacro(EnableSource2,int);

//BTX
  #ifdef VTK_USE_MPI
//ETX
    // Description:
    // Set/Get the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);
//BTX
  #endif
//ETX

  protected:

     vtkTemporalStreamTracer();
    ~vtkTemporalStreamTracer();

    //
    // Make sure the pipeline knows what type we expect as input
    //
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    //
    // The usual suspects
    //
    virtual int ProcessRequest(vtkInformation* request,
                               vtkInformationVector** inputVector,
                               vtkInformationVector* outputVector);

    //
    // Store any information we need in the output and fetch what we can
    // from the input
    //
    virtual int RequestInformation(vtkInformation* request,
                                  vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector);

    //
    // Compute input time steps given the output step
    //
    virtual int RequestUpdateExtent(vtkInformation* request,
                                    vtkInformationVector** inputVector,
                                    vtkInformationVector* outputVector);

    //
    // Generate output
    //
    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);


    int InitializeInterpolator(double times[2]);
    int SetupInputs(vtkInformation* inInfo, vtkInformation* outInfo);

//
//BTX
//
    // Description : Tests points to see if they are inside this region
    // Pass in either a source object or an input list, one
    // parameter should be valid, the other NULL
    void InjectSeeds(vtkDataSet *source, int sourceID, int injectionID, 
      ParticleList *inputlist,
      ParticleList &candidates, ParticleList *outofdomain);

    void UpdateSeeds(ParticleList &candidates);

    void TransmitReceiveParticles(
      ParticleList &outofdomain, ParticleList &received, bool removeself);

    void IntegrateParticle(
      ParticleIterator &it, 
      double currenttime, double terminationtime,
      vtkInitialValueProblemSolver* integrator);
//
//ETX
//

    void GenerateOutputLines(vtkPolyData *output);
    bool DoParticleSendTasks(ParticleLifetime &info, double point1[4], double velocity[3], double delT);
    bool DoParticleSendTasks(ParticleLifetime &info, double point1[4], double delT);
    bool ComputeDomainExitLocation(
      double pos[4], double p2[4], double intersection[4],
      vtkGenericCell *cell);
    void AddParticleToMPISendList(ParticleLifetime &info);

    int UpdatePiece;
    int UpdateNumPieces;

    // Support pipeline time 
    unsigned int        TimeStep;
    unsigned int        ActualTimeStep;
    unsigned int        NumberOfInputTimeSteps;
    unsigned int        MaxTrackLength;
    unsigned int        TrackMarkSpaceRatio;
    //
    int                 EnableSource1;      
    int                 EnableSource2;
    //
    int                 AllFixedGeometry;
    int                 NoFixedGeometry;

    //BTX
    vtkstd::vector<double>  InputTimeValues;
    vtkstd::vector<double>  OutputTimeValues;
    //ETX

    // internal data variables
    int                 MaxCellSize;
    double              EarliestTime;
    double              CurrentTimeSteps[2];
    double              TimeStepResolution;
    int                 ForceReinjectionEveryNSteps;
    bool                ReinjectionFlag;  
    int                 ReinjectionCounter;
    int                 MemoryAllocationSize;
    //
    vtkTimeStamp  ParticleInjectionTime;
    vtkTimeStamp  SeedInjectionTime;

//BTX
    unsigned int      NumberOfParticles;
    ParticleDataList  ParticleHistories;
#ifdef JB_H5PART_PARTICLE_OUTPUT
    vtkH5PartWriter *HDF5ParticleWriter;
#endif
//ETX
/*
    vtkstd::vector<Position>            MysteryCoordinates;
*/
//BTX
    //
    vtkstd::vector<double> weights;
    //
    // These are the final Points/Cells that are generated from the above lists
    //
    vtkSmartPointer<vtkPoints>        OutputCoordinates;
    vtkSmartPointer<vtkCellArray>     ParticlePolyLines;
    //
    // Scalar arrays that are generated as each particle is updated
    //
    vtkSmartPointer<vtkDoubleArray>   time;
    vtkSmartPointer<vtkIntArray>      retVals;
    vtkSmartPointer<vtkDoubleArray>   cellVectors;
    vtkSmartPointer<vtkDoubleArray>   vorticity;
    vtkSmartPointer<vtkDoubleArray>   rotation;
    vtkSmartPointer<vtkDoubleArray>   angularVel;
    //
    ParticleList                      MPISendList;
    //
    vtkSmartPointer<vtkTemporalInterpolatedVelocityField>  Interpolator;
    vtkCompositeDataSet                                   *InputDataT[2];
    vtkGenericCell                                        *GenericCell;

    // info about each dataset we will use repeatedly
    typedef struct {
      double b[6];
    } bounds;
    vtkstd::vector<bounds> CachedBounds[2];
    vtkstd::vector<bool>   GeometryFixed[2];

    bool InsideBounds(double point[]);

    //
//ETX

//BTX
  #ifdef VTK_USE_MPI
//ETX
    vtkMultiProcessController* Controller;
//BTX
  #endif
//ETX
  static vtkIdType UniqueIdCounter;

private:
  // Description:
  // Hide this because we require a new interpolator type
  void SetInterpolatorPrototype(vtkInterpolatedVelocityField*) {};

private:
  vtkTemporalStreamTracer(const vtkTemporalStreamTracer&);  // Not implemented.
  void operator=(const vtkTemporalStreamTracer&);  // Not implemented.
};


#endif


