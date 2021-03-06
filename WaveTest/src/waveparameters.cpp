#include <waveparameters.h>
#include <mpicommunicator.h>

// DO NOT MODIFY - This file is automatically generated during compilation.

Result WaveParameters::broadcastParameters(){
  MPICommunicator *comm = MPICommunicator::getInstance();

  unsigned int bcInitialConditions = static_cast<unsigned int>(mInitialConditions);
  comm->broadcastParameter(&bcInitialConditions);
  mInitialConditions = static_cast<InitialConditions>(bcInitialConditions);

  comm->broadcastParameter(&mGridPointsX);

  comm->broadcastParameter(&mGridPointsY);

  comm->broadcastParameter(&mDomainMinX);

  comm->broadcastParameter(&mDomainMaxX);

  comm->broadcastParameter(&mDomainMinY);

  comm->broadcastParameter(&mDomainMaxY);

  comm->broadcastParameter(&mProjectName);

  comm->broadcastParameter(&mTimeStart);

  comm->broadcastParameter(&mTimeEnd);

  comm->broadcastParameter(&mGhostPoints);

  comm->broadcastParameter(&mMinCFL);

  comm->broadcastParameter(&mMaxCFL);

  comm->broadcastParameter(&mErrorTolerance);

  return SUCCESS;}
