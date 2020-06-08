#include <wave.h>
#include <iostream>
#include <waveparameters.h>
#include <cmath>

// Constructor
Wave::Wave(Domain* d, Solver* s) : ODE(2, d, s){
  params = nullptr;

  addField("Evolution", nEqs, true, true);
  setVariableName("Evolution",0,"U_PHI");
  setVariableName("Evolution",1,"U_PI");

  reallocateData();
}

// Destructor
Wave::~Wave(){

}

void Wave::rhs(std::shared_ptr<FieldMap>& fieldMap){
  const Grid& grid = fieldMap->getGrid();
  unsigned int nx = grid.getSize()[0];
  unsigned int ny = grid.getSize()[1];

  unsigned int nb = domain->getGhostPoints();

  double dx = grid.getSpacing()[0];
  double dy = grid.getSpacing()[1];

  double **dudt = fieldMap->getSolverField("Evolution")->getCurrentRHS();
  double **u = fieldMap->getSolverField("Evolution")->getIntermediateData();

  pair2<int> commPartners = domain->getCommPartners();
  unsigned int xstart = (domain->hasBoundary(LEFT)) ? nb + 1 : nb;
  unsigned int xend = (domain->hasBoundary(RIGHT)) ? nx - nb - 1 : nx - nb;
  unsigned int ystart = (domain->hasBoundary(DOWN)) ? nb + 1 : nb;
  unsigned int yend = (domain->hasBoundary(UP)) ? ny - nb - 1 : ny - nb;

  // Zero out the righthand side. This makes fixing the boundaries easier
  // because we don't have to treat the corners specially.
  for(unsigned int j = 0; j < ny; j++){
    for(unsigned int i = 0; i < nx; i++){
      unsigned int xy = grid.getIndex(i,j);
      dudt[U_PHI][xy] = 0.0;
      dudt[U_PI ][xy] = 0.0;
    }
  }

  // Loop over all the interior points.
  for(unsigned int j = ystart; j < yend; j++){
    for(unsigned int i = xstart; i < xend; i++){
      unsigned int xy = grid.getIndex(i,j);
      dudt[U_PHI][xy] = u[U_PI][xy];
      dudt[U_PI ][xy] = (u[U_PHI][grid.getIndex(i+1,j)] - 2.0*u[U_PHI][xy] +
                         u[U_PHI][grid.getIndex(i-1,j)])/(dx*dx) +
                        (u[U_PHI][grid.getIndex(i,j+1)] - 2.0*u[U_PHI][xy] + 
                         u[U_PHI][grid.getIndex(i,j-1)])/(dy*dy);
    }
  }

  // Apply outflow boundary conditions.
  // Left side.
  if(domain->hasBoundary(LEFT)){
    for(unsigned int j = 0; j < ny; j++){
      unsigned int xy = grid.getIndex(nb,j);
      unsigned int x1y = grid.getIndex(nb + 1, j);
      unsigned int x2y = grid.getIndex(nb + 2, j);
      dudt[U_PHI][xy] += (-3.0*u[U_PHI][xy] + 4.0*u[U_PHI][x1y] - u[U_PHI][x2y])/(2.0*dx);
      dudt[U_PI ][xy] += (-3.0*u[U_PI ][xy] + 4.0*u[U_PI ][x1y] - u[U_PI ][x2y])/(2.0*dx);
    }
  }
  // Right side.
  if(domain->hasBoundary(RIGHT)){
    for(unsigned int j = 0; j < ny; j++){
      unsigned int xy = grid.getIndex(nx - nb - 1,j);
      unsigned int x1y = grid.getIndex(nx - nb - 2, j);
      unsigned int x2y = grid.getIndex(nx - nb - 3, j);
      dudt[U_PHI][xy] += (-3.0*u[U_PHI][xy] + 4.0*u[U_PHI][x1y] - u[U_PHI][x2y])/(2.0*dx);
      dudt[U_PI ][xy] += (-3.0*u[U_PI ][xy] + 4.0*u[U_PI ][x1y] - u[U_PI ][x2y])/(2.0*dx);
    }
  }
  // Bottom side.
  if(domain->hasBoundary(DOWN)){
    for(unsigned int i = 0; i < nx; i++){
      unsigned int xy = grid.getIndex(i,nb);
      unsigned int x1y = grid.getIndex(i, nb + 1);
      unsigned int x2y = grid.getIndex(i, nb + 2);
      dudt[U_PHI][xy] += (-3.0*u[U_PHI][xy] + 4.0*u[U_PHI][x1y] - u[U_PHI][x2y])/(2.0*dy);
      dudt[U_PI ][xy] += (-3.0*u[U_PI ][xy] + 4.0*u[U_PI ][x1y] - u[U_PI ][x2y])/(2.0*dy);
    }
  }
  // Top side.
  if(domain->hasBoundary(UP)){
    for(unsigned int i = 0; i < nx; i++){
      unsigned int xy = grid.getIndex(i,ny - nb - 1);
      unsigned int x1y = grid.getIndex(i, ny - nb - 2);
      unsigned int x2y = grid.getIndex(i, ny - nb - 3);
      dudt[U_PHI][xy] += (-3.0*u[U_PHI][xy] + 4.0*u[U_PHI][x1y] - u[U_PHI][x2y])/(2.0*dy);
      dudt[U_PI ][xy] += (-3.0*u[U_PI ][xy] + 4.0*u[U_PI ][x1y] - u[U_PI ][x2y])/(2.0*dy);
    }
  }
}

void Wave::initData(){
  pair2<double> bounds = domain->getBounds();
  double x0 = 0.5*(bounds[0][1] + bounds[0][0]);
  double y0 = 0.5*(bounds[1][1] + bounds[1][0]);
  double amp = 1.0;
  double sigma = 0.125;
  // Just assume a Gaussian for the time being.
  unsigned int nx = domain->getGrid()->getSize()[0];
  unsigned int ny = domain->getGrid()->getSize()[1];

  auto evol = (*fieldData)["Evolution"];
  auto points = domain->getGrid()->getPoints();
  const double *x = points[0];
  const double *y = points[1];
  double **u = evol->getData();
  for(unsigned int j = 0; j < ny; j++){
    for(unsigned int i = 0; i < nx; i++){
      unsigned int xy = domain->getGrid()->getIndex(i,j);
      u[U_PHI][xy] = amp*std::exp(-((x[i] - x0)*(x[i] - x0) + (y[j] - y0)*(y[j] - y0))/(sigma*sigma));
      u[U_PI ][xy] = 0.0;
    }
  }
}

void Wave::setParameters(WaveParameters *p){
  params = p;
}

WaveParameters* Wave::getParameters(){
  return params;
}
