#include <mpicommunicator.h>
#include <testparser.h>
#include <testparameters.h>
#include <iostream>
#include <cstdio>

using namespace std;

void printTest(char *name, bool success,int rank){
  if(success){
    printf("\033[1;32m%s test passed on rank %d.\033[0m\n\n",name,rank);
  }
  else{
    printf("\033[1;31m%s test failed on rank %d.\033[0m\n\n",name,rank);
  }
}

int main(int argc, char *argv[]){
  MPICommunicator *comm = MPICommunicator::getInstance();
  Result result = comm->init();
  if(result != SUCCESS){
    cout << "There was an error initializing MPICommunicator.\n";
    return 0;
  }

  int rank = comm->getRank();
  int root = comm->getRootRank();
  
  // Load the parameter file.
  TestParameters params;
  if(rank == root){
    if(argc < 2){
      std::cout << "Usage: ./ParameterTest <parameter file>\n";
      return 0;
    }
    TestParser parser;
    parser.updateParameters(argv[1],&params);
  }

  // Broadcast all the parameters.
  params.broadcastParameters();

  // Validate the results.
  printTest("InitialConditions",params.getInitialConditions() == TestParameters::FLAT,rank);
  printTest("GridPointsX",params.getGridPointsX() == 201,rank);
  printTest("GridPointsY",params.getGridPointsY() == 201,rank);
  printTest("DomainMinX",params.getDomainMinX() == -0.5,rank);
  printTest("DomainMinY",params.getDomainMinY() == -0.5,rank);
  printTest("DomainMaxX",params.getDomainMaxX() == 0.5,rank);
  printTest("DomainMaxY",params.getDomainMaxY() == 0.5,rank);
  printTest("ProjectName",params.getProjectName().compare("Working Project") == 0,rank);

  result = comm->cleanup();
  if(result != SUCCESS){
    cout << "There was an error cleaning up MPICommunicator.\n";
  }
  return 0;
}
