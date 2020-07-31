#include <mpicommunicator.h>
#include <cstring>

MPICommunicator MPICommunicator::comm = MPICommunicator();

// Constructor {{{
MPICommunicator::MPICommunicator(){
  isInitialized = false;
  worldSize = -1;
  worldRank = -1;
  rootRank = -1;
  processorName = "";
}
// }}}

// Destructor {{{
MPICommunicator::~MPICommunicator(){

}
// }}}

// init {{{
Result MPICommunicator::init(){
  if(isInitialized){
    return SUCCESS;
  }
  if(MPI_Init(NULL, NULL) != MPI_SUCCESS){
    return FAILURE;
  }

  // Get the number of ranks and the current rank.
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  // Get the processor name.
  char buffer[MPI_MAX_PROCESSOR_NAME];
  int bufferLength;
  MPI_Get_processor_name(buffer, &bufferLength);
  processorName = std::string(buffer, bufferLength);
  rootRank = 0;

  isInitialized = true;

  return SUCCESS;
}
// }}}

// cleanup {{{
Result MPICommunicator::cleanup(){
  if(!isInitialized){
    return UNINITIALIZED;
  }
  if(MPI_Finalize() != MPI_SUCCESS){
    return FAILURE;
  }
  worldRank = -1;
  worldSize = -1;
  rootRank = -1;
  isInitialized = false;
  return SUCCESS;
}
// }}}

// broadcastParameter {{{
Result MPICommunicator::broadcastParameter(int *par){
  if(!isInitialized){
    return UNINITIALIZED;
  }
  int result = MPI_Bcast(par, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
  if(result){
    return FAILURE;
  }
  return SUCCESS;
}

Result MPICommunicator::broadcastParameter(unsigned int *par){
  if(!isInitialized){
    return UNINITIALIZED;
  }
  int result = MPI_Bcast(par, 1, MPI_UNSIGNED, rootRank, MPI_COMM_WORLD);
  if(result){
    return FAILURE;
  }
  return SUCCESS;
}

Result MPICommunicator::broadcastParameter(double *par){
  if(!isInitialized){
    return UNINITIALIZED;
  }
  int result = MPI_Bcast(par, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
  if(result){
    return FAILURE;
  }
  return SUCCESS;
}

Result MPICommunicator::broadcastParameter(std::string *par){
  if(!isInitialized){
    return UNINITIALIZED;
  }
  // Strings are a little different. We first need to convert it to a C-style string.
  char *str;
  unsigned int len;
  if(worldRank == rootRank){
    str = new char[par->length()+1];
    std::strcpy(str,par->c_str());
    len = par->length()+1;
  }
  int result = MPI_Bcast(&len, 1, MPI_UNSIGNED, rootRank, MPI_COMM_WORLD);
  if(!result){
    if(worldRank != rootRank){
      str = new char[len];
    }
    result = MPI_Bcast(str, len, MPI_CHAR, rootRank, MPI_COMM_WORLD);
  }
  Result r;
  if(result){
    r = FAILURE;
  }
  else{
    r = SUCCESS;
    if(worldRank != rootRank){
      (*par) = std::string(str);
      delete[] str;
    }
  }
  
  if(worldRank == rootRank){
    delete[] str;
  }
  return r;
} 
// }}}

// sendData {{{
/*
Result MPICommunicator::sendData(double *buffer, int size, int dest){
  MPI_Request request;
  MPI_Status status;
  return SUCCESS;
}
*/
// }}}

// allSum {{{
Result MPICommunicator::allSum(double sendData, double &recvData){
  if(!isInitialized){
    return UNINITIALIZED;
  }

  //MPI_Allreduce(&sendData, &recvData, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double *sum;
  if(worldRank == rootRank){
    sum = new double[worldSize];
  }

  int result = MPI_Gather(&sendData, 1, MPI_DOUBLE, sum, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);

  if(result){
    if(worldRank == rootRank){
      delete[] sum;
      return FAILURE;
    }
  }

  if(worldRank == rootRank){
    double prev;
    recvData = 0.0;
    double err = 0.0;
    for(unsigned int i = 0; i < worldSize; i++){
      prev = recvData;
      recvData += sum[i] + err;
      err = sum[i] - ((recvData - prev) - err);
    }
  }

  result = MPI_Bcast(&recvData, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
  Result r = (result == 0) ? SUCCESS : FAILURE;

  if(worldRank == rootRank){
    delete[] sum;
  }
  return r;
} 
// }}}

// findMax {{{
Result MPICommunicator::findMax(double in, double &out){
  if(!isInitialized){
    return UNINITIALIZED;
  }

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return SUCCESS;
}
// }}}

// findMin {{{
Result MPICommunicator::findMin(double in, double &out){
  if(!isInitialized){
    return UNINITIALIZED;
  }

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  
  return SUCCESS;
}
// }}}
