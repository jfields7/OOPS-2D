#include <domain.h>
#include <mpicommunicator.h>
#include <meshparser.h>
#include <meshparameters.h>
#include <iostream>
#include <cstdio>

using namespace std;

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
  MeshParameters params;
  if(rank == root){
    if(argc < 2){
      std::cout << "Usage: ./ParameterTest <parameter file>\n";
      return 0;
    }
    MeshParser parser;
    parser.updateParameters(argv[1],&params);
  }

  // Broadcast all the parameters.
  params.broadcastParameters();

  // Let's set up a domain!
  Domain domain;
  pair2<double> bounds;
  bounds[0][0] = params.getDomainMinX();
  bounds[0][1] = params.getDomainMaxX();
  bounds[1][0] = params.getDomainMinY();
  bounds[1][1] = params.getDomainMaxY();
  domain.setBounds(bounds);

  unsigned int shp[2] = {params.getGridPointsX(), params.getGridPointsY()};

  result = domain.buildMesh(shp);

  if(result != SUCCESS){
    std::cout << "There was an error building the mesh.\n";
    std::cout << "  Error Code " << result << "\n";
  }

  const Grid *grid = domain.getGrid();
  const pair2<double> gbounds = grid->getBounds();
  printf("Grid acquired on rank %d", rank);
  printf("  Bounds: [%g,%g]x[%g,%g]\n",gbounds[0][0], gbounds[0][1], gbounds[1][0], gbounds[1][1]);

  comm->cleanup();
  return 0;
}
