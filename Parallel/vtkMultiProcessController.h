/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMultiProcessController.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMultiProcessController - Multiprocessing communication superclass
// .SECTION Description
// vtkMultiProcessController is used to control multiple processes
// in a distributed computing environment. It has
// methods for executing single/multiple method(s) on multiple processors,
// triggering registered callbacks (Remote Methods) (AddRMI(), TriggerRMI())
// and communication. Please note that the communication is done using
// the communicator which is accessible to the user. Therefore it is
// possible to get the communicator with GetCommunicator() and use
// it to send and receive data. This is the encoured communication method.
// The internal (RMI) communications are done using a second internal
// communicator (called RMICommunicator).
//
// .SECTION see also
// vtkMPIController
// vtkCommunicator vtkMPICommunicator

#ifndef __vtkMultiProcessController_h
#define __vtkMultiProcessController_h

#include "vtkObject.h"

#include "vtkCommunicator.h" // Needed for direct access to communicator

class vtkDataSet;
class vtkImageData;
class vtkCollection;
class vtkOutputWindow;
class vtkDataObject;
class vtkMultiProcessController;

//BTX
// The type of function that gets called when new processes are initiated.
typedef void (*vtkProcessFunctionType)(vtkMultiProcessController *controller, 
                                       void *userData);

// The type of function that gets called when an RMI is triggered.
typedef void (*vtkRMIFunctionType)(void *localArg, 
                                   void *remoteArg, int remoteArgLength, 
                                   int remoteProcessId);
//ETX


class VTK_PARALLEL_EXPORT vtkMultiProcessController : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkMultiProcessController,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // This method is for setting up the processes.
  // If a subclass needs to initialize process communication (i.e. MPI)
  // it would over ride this method.
  virtual void Initialize(int* vtkNotUsed(argc), char*** vtkNotUsed(argv))=0;

  // Description:
  // This method is for setting up the processes.
  // If a subclass needs to initialize process communication (i.e. MPI)
  // it would over ride this method.  Provided for initialization outside vtk.
  virtual void Initialize(int* vtkNotUsed(argc), char*** vtkNotUsed(argv),
                          int initializedExternally)=0;

  // Description:
  // This method is for cleaning up.
  // If a subclass needs to clean up process communication (i.e. MPI)
  // it would over ride this method.
  virtual void Finalize()=0;

  // Description:
  // This method is for cleaning up.
  // If a subclass needs to clean up process communication (i.e. MPI)
  // it would over ride this method.  Provided for finalization outside vtk.
  virtual void Finalize(int finalizedExternally)=0;

  // Description:
  // Set the number of processes you will be using.  This defaults
  // to the maximum number available.  If you set this to a value
  // higher than the default, you will get an error.
  void SetNumberOfProcesses(int num);
  int GetNumberOfProcesses();

  //BTX
  // Description:
  // Set the SingleMethod to f() and the UserData of the
  // for the method to be executed by all of the processes
  // when SingleMethodExecute is called.  All the processes will
  // start by calling this function.
  void SetSingleMethod(vtkProcessFunctionType, void *data);
  //ETX

  // Description:
  // Execute the SingleMethod (as define by SetSingleMethod) using
  // this->NumberOfProcesses processes.  This will only return when
  // all the processes finish executing their methods.
  virtual void SingleMethodExecute() = 0;
  
  //BTX
  // Description:
  // Set the MultipleMethod to f() and the UserData of the
  // for the method to be executed by the process index
  // when MultipleMethodExecute is called.  This is for having each 
  // process start with a different function and data argument.
  void SetMultipleMethod(int index, vtkProcessFunctionType, void *data); 
  //ETX

  // Description:
  // Execute the MultipleMethods (as define by calling SetMultipleMethod
  // for each of the required this->NumberOfProcesses methods) using
  // this->NumberOfProcesses processes.
  virtual void MultipleMethodExecute() = 0;

  // Description:
  // Tells you which process [0, NumProcess) you are in.
  int GetLocalProcessId();

  // Description:
  // This convenience method returns the controller associated with the 
  // local process.  It returns NULL until the processes are spawned.
  // It is better if you hang on to the controller passed as an argument to the
  // SingleMethod or MultipleMethod functions.
  static vtkMultiProcessController *GetGlobalController();

  // Description:
  // This method can be used to tell the controller to create
  // a special output window in which all messages are preceded
  // by the process id.
  virtual void CreateOutputWindow() = 0;
  
  //------------------ RMIs --------------------
  //BTX
  // Description:
  // Register remote method invocation in the receiving process
  // which makes the call.  It must have a unique tag as an RMI id.
  // The vtkRMIFunctionType has several arguments: localArg (same as passed in),
  // remoteArg, remoteArgLength (memory passed by process triggering the RMI),
  // remoteProcessId.
  // Since only one callback can be registered per tag, this method will remove
  // any previously registered callback for the given tag.
  // Returns a unique Id for the RMI registration which can be used to
  // unregister the callback. RemoveRMI() should be preferred over
  // RemoveFirstRMI() since it avoid accidental removal of callbacks.
  unsigned long AddRMI(vtkRMIFunctionType, void *localArg, int tag);
  
  // Description:
  // Remove the first RMI matching the tag.
  int RemoveFirstRMI(int tag);

  // Description:
  // Remove the  RMI matching the id. The id is the same id returned by
  // AddRMI().
  int RemoveRMI(unsigned long id);

  // Description:
  // Take an RMI away.
  void RemoveRMI(vtkRMIFunctionType f, void *arg, int tag)
    {f = f; arg = arg; tag = tag; vtkErrorMacro("RemoveRMI Not Implemented Yet");};
  //ETX
  
  // Description:
  // A method to trigger a method invocation in another process.
  void TriggerRMI(int remoteProcessId, void *arg, int argLength, int tag);

  // Description:
  // A conveniance method.  Called on process 0 to break "ProcessRMIs" loop
  // on all other processes.
  void TriggerBreakRMIs();

  // Description:
  // Convenience method when the arg is a string. 
  void TriggerRMI(int remoteProcessId, const char *arg, int tag) 
    { this->TriggerRMI(remoteProcessId, (void*)arg, 
                       static_cast<int>(strlen(arg))+1, tag); }

  // Description:
  // Convenience method when there is no argument.
  void TriggerRMI(int remoteProcessId, int tag)
    { this->TriggerRMI(remoteProcessId, NULL, 0, tag); }

  // Description:
  // Calling this method gives control to the controller to start
  // processing RMIs. Possible return values are:
  // RMI_NO_ERROR,
  // RMI_TAG_ERROR : rmi tag could not be received,
  // RMI_ARG_ERROR : rmi arg could not be received.
  // If reportErrors is false, no vtkErrorMacro is called.
  // ProcessRMIs() calls ProcessRMIs(int) with reportErrors = 0.
  // If dont_loop is 1, this call just process one RMI message
  // and exits.
  int ProcessRMIs(int reportErrors, int dont_loop = 0);
  int ProcessRMIs();
 
  // Description:
  // Setting this flag to 1 will cause the ProcessRMIs loop to return.
  // This also causes vtkUpStreamPorts to return from
  // their WaitForUpdate loops.
  vtkSetMacro(BreakFlag, int);
  vtkGetMacro(BreakFlag, int);

  // Description:
  // Returns the communicator associated with this controller.
  // A default communicator is created in constructor.
  vtkGetObjectMacro(Communicator, vtkCommunicator);

  // Description:
  // Accessor to some default tags.
  static int GetBreakRMITag() { return BREAK_RMI_TAG; }
  static int GetRMITag() { return RMI_TAG; }
  static int GetRMIArgTag() { return RMI_ARG_TAG; }  
  
//BTX

  enum Errors 
  {
    RMI_NO_ERROR,
    RMI_TAG_ERROR,
    RMI_ARG_ERROR
  };

  enum Consts 
  {
    MAX_PROCESSES  = 8192,
    ANY_SOURCE     = -1,
    INVALID_SOURCE = -2
  };

  enum Tags 
  {
    RMI_TAG        = 1,
    RMI_ARG_TAG    = 2,
    BREAK_RMI_TAG  = 3,
    XML_WRITER_DATA_INFO = 4
  };

//ETX

  // Description:
  // This method can be used to synchronize processes.
  virtual void Barrier() = 0;

  static void SetGlobalController(vtkMultiProcessController *controller);

  //------------------ Communication --------------------
  
  // Description:
  // This method sends data to another process.  Tag eliminates ambiguity
  // when multiple sends or receives exist in the same process.
  int Send(const int* data, vtkIdType length, int remoteProcessId, int tag);
  int Send(const unsigned long* data, vtkIdType length, int remoteProcessId, 
           int tag);
  int Send(const char* data, vtkIdType length, int remoteProcessId, int tag);
  int Send(const unsigned char* data, vtkIdType length, int remoteProcessId, int tag);
  int Send(const float* data, vtkIdType length, int remoteProcessId, int tag);
  int Send(const double* data, vtkIdType length, int remoteProcessId, int tag);
#ifdef VTK_USE_64BIT_IDS
  int Send(const vtkIdType* data, vtkIdType length, int remoteProcessId, int tag);
#endif
  int Send(vtkDataObject *data, int remoteId, int tag);
  int Send(vtkDataArray *data, int remoteId, int tag);

  // Description:
  // This method receives data from a corresponding send. It blocks
  // until the receive is finished.  It calls methods in "data"
  // to communicate the sending data.
  int Receive(int* data, vtkIdType length, int remoteProcessId, int tag);
  int Receive(unsigned long* data, vtkIdType length, int remoteProcessId, 
              int tag);
  int Receive(char* data, vtkIdType length, int remoteProcessId, int tag);
  int Receive(unsigned char* data, vtkIdType length, int remoteProcessId, int tag);
  int Receive(float* data, vtkIdType length, int remoteProcessId, int tag);
  int Receive(double* data, vtkIdType length, int remoteProcessId, int tag);
#ifdef VTK_USE_64BIT_IDS
  int Receive(vtkIdType* data, vtkIdType length, int remoteProcessId, int tag);
#endif
  int Receive(vtkDataObject* data, int remoteId, int tag);
  int Receive(vtkDataArray* data, int remoteId, int tag);

  vtkDataObject *ReceiveDataObject(int remoteId, int tag);

  //---------------------- Collective Operations ----------------------

  // Description:
  // Broadcast sends the array in the process with id \c srcProcessId to all of
  // the other processes.  All processes must call these method with the same
  // arguments in order for it to complete.
  int Broadcast(int *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
  int Broadcast(unsigned long *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
  int Broadcast(unsigned char *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
  int Broadcast(char *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
  int Broadcast(float *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
  int Broadcast(double *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int Broadcast(vtkIdType *data, vtkIdType length, int srcProcessId) {
    return this->Communicator->Broadcast(data, length, srcProcessId);
  }
#endif
  int Broadcast(vtkDataObject *data, int srcProcessId) {
    return this->Communicator->Broadcast(data, srcProcessId);
  }
  int Broadcast(vtkDataArray *data, int srcProcessId) {
    return this->Communicator->Broadcast(data, srcProcessId);
  }

  // Description:
  // Gather collects arrays in the process with id \c destProcessId.  Each
  // process (including the destination) sends the contents of its send buffer
  // to the destination process.  The destination process receives the
  // messages and stores them in rank order.  The \c length argument
  // (which must be the same on all processes) is the length of the
  // sendBuffers.  The \c recvBuffer (on te destination process) must be of
  // length length*numProcesses.  Gather is the inverse operation of Scatter.
  int Gather(const int *sendBuffer, int *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
  int Gather(const unsigned long *sendBuffer, unsigned long *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
  int Gather(const unsigned char *sendBuffer, unsigned char *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
  int Gather(const char *sendBuffer, char *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
  int Gather(const float *sendBuffer, float *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
  int Gather(const double *sendBuffer, double *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int Gather(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
             vtkIdType length, int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, length,
                                      destProcessId);
  }
#endif
  int Gather(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer,
             int destProcessId) {
    return this->Communicator->Gather(sendBuffer, recvBuffer, destProcessId);
  }

  // Description:
  // GatherV is the vector variant of Gather.  It extends the functionality of
  // Gather by allowing a varying count of data from each process.
  // GatherV collects arrays in the process with id \c destProcessId.  Each
  // process (including the destination) sends the contents of its send buffer
  // to the destination process.  The destination process receives the
  // messages and stores them in rank order.  The \c sendLength argument
  // defines how much the local process sends to \c destProcessId and
  // \c recvLengths is an array containing the amount \c destProcessId
  // receives from each process, in rank order.
  int GatherV(const int* sendBuffer, int* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
  int GatherV(const unsigned long* sendBuffer, unsigned long* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
  int GatherV(const unsigned char* sendBuffer, unsigned char* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
  int GatherV(const char* sendBuffer, char* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
  int GatherV(const float* sendBuffer, float* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
  int GatherV(const double* sendBuffer, double* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int GatherV(const vtkIdType* sendBuffer, vtkIdType* recvBuffer, 
              vtkIdType sendLength, vtkIdType* recvLengths, vtkIdType* offsets,
              int destProcessId) {
    return this->Communicator->GatherV(sendBuffer, recvBuffer,
                                       sendLength, recvLengths,
                                       offsets, destProcessId);
  }
#endif

  // Description:
  // Scatter takes an array in the process with id \c srcProcessId and
  // distributes it.  Each process (including the source) receives a portion of
  // the send buffer.  Process 0 receives the first \c length values, process 1
  // receives the second \c length values, and so on.  Scatter is the inverse
  // operation of Gather.
  int Scatter(const int *sendBuffer, int *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
  int Scatter(const unsigned long *sendBuffer, unsigned long *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
  int Scatter(const unsigned char *sendBuffer, unsigned char *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
  int Scatter(const char *sendBuffer, char *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
  int Scatter(const float *sendBuffer, float *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
  int Scatter(const double *sendBuffer, double *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int Scatter(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
              vtkIdType length, int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, length,
                                       srcProcessId);
  }
#endif
  int Scatter(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer,
              int srcProcessId) {
    return this->Communicator->Scatter(sendBuffer, recvBuffer, srcProcessId);
  }

  // Description:
  // ScatterV is the vector variant of Scatter.  It extends the functionality of
  // Scatter by allowing a varying count of data to each process.
  // ScatterV takes an array in the process with id \c srcProcessId and
  // distributes it.  Each process (including the source) receives a portion of
  // the send buffer defined by the \c sendLengths and \c offsets arrays.
  int ScatterV(const int *sendBuffer, int *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
  int ScatterV(const unsigned long *sendBuffer, unsigned long *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
  int ScatterV(const unsigned char *sendBuffer, unsigned char *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
  int ScatterV(const char *sendBuffer, char *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
  int ScatterV(const float *sendBuffer, float *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
  int ScatterV(const double *sendBuffer, double *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int ScatterV(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
               vtkIdType *sendLengths, vtkIdType *offsets,
               vtkIdType recvLength, int srcProcessId) {
    return this->Communicator->ScatterV(sendBuffer, recvBuffer,
                                        sendLengths, offsets, recvLength,
                                        srcProcessId);
  }
#endif

  // Description:
  // Same as gather except that the result ends up on all processes.
  int AllGather(const int *sendBuffer, int *recvBuffer, vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
  int AllGather(const unsigned long *sendBuffer,
                unsigned long *recvBuffer, vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
  int AllGather(const unsigned char *sendBuffer,
                unsigned char *recvBuffer, vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
  int AllGather(const char *sendBuffer, char *recvBuffer, vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
  int AllGather(const float *sendBuffer, float *recvBuffer, vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
  int AllGather(const double *sendBuffer,
                double *recvBuffer, vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
#ifdef VTK_USE_64BIT_IDS
  int AllGather(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
                vtkIdType length) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer, length);
  }
#endif
  int AllGather(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer) {
    return this->Communicator->AllGather(sendBuffer, recvBuffer);
  }

  // Description:
  // Same as GatherV except that the result is placed in all processes.
  int AllGatherV(const int* sendBuffer, int* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
  int AllGatherV(const unsigned long* sendBuffer, unsigned long* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
  int AllGatherV(const unsigned char* sendBuffer, unsigned char* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
  int AllGatherV(const char* sendBuffer, char* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
  int AllGatherV(const float* sendBuffer, float* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
  int AllGatherV(const double* sendBuffer, double* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
#ifdef VTK_USE_64BIT_IDS
  int AllGatherV(const vtkIdType* sendBuffer, vtkIdType* recvBuffer, 
                 vtkIdType sendLength, vtkIdType* recvLengths,
                 vtkIdType* offsets) {
    return this->Communicator->AllGatherV(sendBuffer, recvBuffer,
                                          sendLength, recvLengths,
                                          offsets);
  }
#endif

  // Description:
  // Reduce an array to the given destination process.  This version of Reduce
  // takes an identifier defined in the
  // vtkCommunicator::StandardOperations enum to define the operation.
  int Reduce(const int *sendBuffer, int *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const unsigned long *sendBuffer, unsigned long *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const unsigned char *sendBuffer, unsigned char *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const char *sendBuffer, char *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const float *sendBuffer, float *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const double *sendBuffer, double *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int Reduce(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
             vtkIdType length, int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
#endif
  int Reduce(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer,
             int operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer,
                                      operation, destProcessId);
  }

//BTX
  // Description:
  // Reduce an array to the given destination process.  This version of Reduce
  // takes a custom operation as a subclass of vtkCommunicator::Operation.
  int Reduce(const int *sendBuffer, int *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const unsigned long *sendBuffer, unsigned long *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const unsigned char *sendBuffer, unsigned char *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const char *sendBuffer, char *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const float *sendBuffer, float *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
  int Reduce(const double *sendBuffer, double *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
#ifdef VTK_USE_64BIT_IDS
  int Reduce(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
             vtkIdType length, vtkCommunicator::Operation *operation,
             int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer, length,
                                      operation, destProcessId);
  }
#endif
  int Reduce(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer,
             vtkCommunicator::Operation *operation, int destProcessId) {
    return this->Communicator->Reduce(sendBuffer, recvBuffer,
                                      operation, destProcessId);
  }
//ETX

  // Description:
  // Same as Reduce except that the result is placed in all of the processes.
  int AllReduce(const int *sendBuffer, int *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const unsigned long *sendBuffer, unsigned long *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const unsigned char *sendBuffer, unsigned char *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const char *sendBuffer, char *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const float *sendBuffer, float *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const double *sendBuffer, double *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
#ifdef VTK_USE_64BIT_IDS
  int AllReduce(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
                vtkIdType length, int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
#endif
  int AllReduce(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer,
                int operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, operation);
  }
//BTX
  int AllReduce(const int *sendBuffer, int *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const unsigned long *sendBuffer, unsigned long *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const unsigned char *sendBuffer, unsigned char *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const char *sendBuffer, char *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const float *sendBuffer, float *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
  int AllReduce(const double *sendBuffer, double *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
#ifdef VTK_USE_64BIT_IDS
  int AllReduce(const vtkIdType *sendBuffer, vtkIdType *recvBuffer,
                vtkIdType length, vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, length,
                                         operation);
  }
#endif
  int AllReduce(vtkDataArray *sendBuffer, vtkDataArray *recvBuffer,
                vtkCommunicator::Operation *operation) {
    return this->Communicator->AllReduce(sendBuffer, recvBuffer, operation);
  }
//ETX

// Internally implemented RMI to break the process loop.

protected:
  vtkMultiProcessController();
  ~vtkMultiProcessController();
  
  vtkProcessFunctionType      SingleMethod;
  void                       *SingleData;
  vtkProcessFunctionType      MultipleMethod[MAX_PROCESSES];
  void                       *MultipleData[MAX_PROCESSES];  
  
  vtkCollection *RMIs;
  
  // This is a flag that can be used by the ports to break
  // their update loop. (same as ProcessRMIs)
  int BreakFlag;

  void ProcessRMI(int remoteProcessId, void *arg, int argLength, int rmiTag);

  // This method implements "GetGlobalController".  
  // It needs to be virtual and static.
  virtual vtkMultiProcessController *GetLocalController();

  
  // This flag can force deep copies during send.
  int ForceDeepCopy;

  vtkOutputWindow* OutputWindow;

  // Note that since the communicators can be created differently
  // depending on the type of controller, the subclasses are
  // responsible of deleting them.
  vtkCommunicator* Communicator;

  // Communicator which is a copy of the current user
  // level communicator except the context; i.e. even if the tags 
  // are the same, the RMI messages will not interfere with user 
  // level messages. 
  // Note that since the communicators can be created differently
  // depending on the type of controller, the subclasses are
  // responsible of deleting them.
  vtkCommunicator* RMICommunicator;

private:
  vtkMultiProcessController(const vtkMultiProcessController&);  // Not implemented.
  void operator=(const vtkMultiProcessController&);  // Not implemented.

  unsigned long RMICount;
};


inline int vtkMultiProcessController::Send(vtkDataObject *data, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(vtkDataArray *data, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(const int* data, vtkIdType length, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(const unsigned long* data, 
                                           vtkIdType length,
                                           int remoteProcessId,
                                           int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(const char* data, vtkIdType length, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(const unsigned char* data,
                                           vtkIdType length, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(const float* data, vtkIdType length, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Send(const double* data, vtkIdType length,
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

#ifdef VTK_USE_64BIT_IDS
inline int vtkMultiProcessController::Send(const vtkIdType* data,
                                           vtkIdType length, 
                                           int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Send(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}
#endif

inline int vtkMultiProcessController::Receive(vtkDataObject* data, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline vtkDataObject* vtkMultiProcessController::ReceiveDataObject(
  int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->ReceiveDataObject(remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(vtkDataArray* data, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(int* data, vtkIdType length, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(unsigned long* data, 
                                              vtkIdType length,
                                              int remoteProcessId, 
                                              int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(char* data, vtkIdType length, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(unsigned char* data,
                                              vtkIdType length, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(float* data, vtkIdType length, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

inline int vtkMultiProcessController::Receive(double* data, vtkIdType length, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}

#ifdef VTK_USE_64BIT_IDS
inline int vtkMultiProcessController::Receive(vtkIdType* data,
                                              vtkIdType length, 
                                              int remoteProcessId, int tag)
{
  if (this->Communicator)
    {
    return this->Communicator->Receive(data, length, remoteProcessId, tag);
    }
  else
    {
    return 0;
    }
}
#endif

#endif
