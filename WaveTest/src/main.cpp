#include <mpicommunicator.h>
#include <waveparser.h>
#include <waveparameters.h>
#include <wave.h>
#include <rk4.h>
#include <iostream>
#include <cstdio>
#include <cmath>

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
  WaveParameters params;
  if(rank == root){
    if(argc < 2){
      std::cout << "Usage: ./WaveTest <parameter file>\n";
      return 0;
    }
    WaveParser parser;
    parser.updateParameters(argv[1], &params);
  }

  // Broadcast all the parameters.
  params.broadcastParameters();

  // Let's set up our PDE!
  Domain domain = Domain();
  pair2<double> bounds;
  bounds[0][0] = params.getDomainMinX();
  bounds[0][1] = params.getDomainMaxX();
  bounds[1][0] = params.getDomainMinY();
  bounds[1][1] = params.getDomainMaxY();
  domain.setBounds(bounds);
  domain.setGhostPoints(params.getGhostPoints());

  unsigned int shp[2] = {params.getGridPointsX(), params.getGridPointsY()};

  result = domain.buildMesh(shp);

  if(result != SUCCESS){
    std::cout << "There was an error building the mesh.\n";
    std::cout << "  Error Code " << result << "\n";
  }

  // Set up our ODE system.
  RK4 rk4 = RK4();
  Wave ode = Wave(&domain, &rk4);
  ode.setParameters(&params);
  ode.initData();
  ode.setVariableOutput("Evolution",0,true);
  ode.setVariableOutput("Evolution",1,true);

  double ti = params.getTimeStart();
  double tf = params.getTimeEnd();
  auto dx = domain.getGrid()->getSpacing();
  double dt = domain.getCFL()*fmin(dx[0], dx[1]);
  unsigned int M = (tf - ti)/dt;

  //ode.dumpField("Evolution","phi00000.csv", 0, 0);
  ode.outputVTK("phi00000", 0);
  for(unsigned int i = 0; i < M; i++){
    double t = (i + 1)*dt;
    ode.evolveStep(dt);

    char buffer[12];
    sprintf(buffer, "phi%05d",i+1);
    //ode.dumpField("Evolution",buffer, t, 0);
    ode.outputVTK(buffer, t);
  }

  result = comm->cleanup();
  if(result != SUCCESS){
    cout << "There was an error cleaning up MPICommunicator.\n";
  }
  return 0;
}
