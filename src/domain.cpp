#include <domain.h>
#include <types.h>
#include <mpicommunicator.h>
#include <cmath>

// Constructor {{{
Domain::Domain(){
  bounds[0][0] = 0.0;
  bounds[0][1] = 1.0;
  bounds[1][0] = 0.0;
  bounds[1][1] = 1.0;
  d_cfl = 0.5;
  nghosts = 3;
  unsigned int shp[2] = {11, 11};
  grid = new Grid(bounds, shp, nghosts);
  commPartners[0][0] = -1;
  commPartners[0][1] = -1;
  commPartners[1][0] = -1;
  commPartners[1][1] = -1;
  location[0] = 0;
  location[1] = 0;
  bflag = 0;
}
// }}}

// Destructor {{{
Domain::~Domain(){
  delete grid;
}
// }}}

// setBounds {{{
void Domain::setBounds(const pair2<double>& bounds){
  // Make sure that the bounds are ordered correctly.
  // Otherwise, flip them.
  if(bounds[0][0] < bounds[0][1]){
    this->bounds[0][0] = bounds[0][0];
    this->bounds[0][1] = bounds[0][1];
  }
  else{
    this->bounds[0][0] = bounds[0][1];
    this->bounds[0][1] = bounds[0][0];
  }
  if(bounds[1][0] < bounds[1][1]){
    this->bounds[1][0] = bounds[1][0];
    this->bounds[1][1] = bounds[1][1];
  }
  else{
    this->bounds[1][0] = bounds[1][1];
    this->bounds[1][1] = bounds[1][0];
  }
}
// }}}

// buildGrid {{{
Result Domain::buildGrid(unsigned int shp[2], pair2<double> &bounds){
  if(bounds[0][0] < bounds[0][1]){
    if(bounds[0][0] < this->bounds[0][0] || bounds[0][1] > this->bounds[0][1]){
      return OUT_OF_BOUNDS;
    }
  }
  if(bounds[1][0] < bounds[1][1]){
    if(bounds[1][0] < this->bounds[1][0] || bounds[1][1] > this->bounds[1][1]){
      return OUT_OF_BOUNDS;
    }
  }
  //grid = Grid(bounds, shp, nghosts);
  delete grid;
  grid = new Grid(bounds, shp, nghosts);

  return SUCCESS;
}
// }}}

// MPI Routines {{{
// buildMesh {{{
Result Domain::buildMesh(unsigned int shp[2]){
  MPICommunicator *comm = MPICommunicator::getInstance();
  if(!comm->getIsInitialized()){
    return UNINITIALIZED;
  }
  int numRanks = comm->getWorldSize();
  int rank = comm->getRank();
  //int root = comm->getRootRank();

  // Divide up the grid based on the number of processors.
  int dims[2] = {0};
  MPI_Dims_create(numRanks, 2, dims);
  
  // Calculate the location and assign this processor's
  // communication partners.
  location[0] = rank % dims[0];
  location[1] = (rank - location[0])/dims[0];
  assignCommunicationPartners(dims);

  // Divide up the mesh and return if it was successful.
  return divideGrids(shp, dims);
}
// }}}

// assignCommunicationPartners {{{
void Domain::assignCommunicationPartners(int dims[2]){
  if(location[0] == 0){
    commPartners[0][0] = -1;
    bflag = bflag | (1 << LEFT);
  }
  else{
    commPartners[0][0] = location[0]-1 + dims[0]*location[1];
  }
  if(location[0] == dims[0]-1){
    commPartners[0][1] = -1;
    bflag = bflag | (1 << RIGHT);
  }
  else{
    commPartners[0][1] = location[0]+1 + dims[0]*location[1];
  }
  if(location[1] == 0){
    commPartners[1][0] = -1;
    bflag = bflag | (1 << DOWN);
  }
  else{
    commPartners[1][0] = location[0] + dims[0]*(location[1] - 1);
  }
  if(location[1] == dims[1]-1){
    commPartners[1][1] = -1;
    bflag = bflag | (1 << UP);
  }
  else{
    commPartners[1][1] = location[0] + dims[0]*(location[1] + 1);
  }
}
// }}}

// divideGrids {{{
Result Domain::divideGrids(unsigned int shp[2], int dims[2]){
  // Check that the number of number of points in each direction
  // will result in equal dimensions.
  double dx = (bounds[0][1] - bounds[0][0])/(shp[0] - 1);
  double dy = (bounds[1][1] - bounds[1][0])/(shp[1] - 1);
  if(fabs(dx - dy) > 1e-15){
    std::cout << "Grid spacing is not equal!\n";
    std::cout << "  dx = " << dx << "\n";
    std::cout << "  dy = " << dy << "\n";
    return UNEQUAL_SPACING;
  }
  // Divide up the mesh.
  pair2<double> grid_bounds;
  // The bounds are best calculated from the pointwise dimensions.
  // Most likely, the number of points in each direction will not
  // be evenly divisible by the the number of processors in each
  // direction. Therefore, the last processor in each direction
  // will be allocated the extra points.
  unsigned int sz[2] = {0};
  unsigned int nx = shp[0]/dims[0];
  unsigned int ny = shp[1]/dims[1];
  grid_bounds[0][0] = bounds[0][0] + location[0]*nx*dx;
  pointBounds[0][0] = location[0]*nx;
  grid_bounds[1][0] = bounds[1][0] + location[1]*ny*dx;
  pointBounds[1][0] = location[1]*ny;
  if(commPartners[0][1] == -1){
    grid_bounds[0][1] = bounds[0][1];
    sz[0] = nx + (shp[0] % nx);
    pointBounds[0][1] = shp[0]-1;
  }
  else{
    grid_bounds[0][1] = bounds[0][0] + (location[0] + 1)*nx*dx;
    pointBounds[0][1] = (location[0] + 1)*nx;
    sz[0] = nx + 1;
  }
  if(commPartners[1][1] == -1){
    grid_bounds[1][1] = bounds[1][1];
    sz[1] = ny + (shp[1] % ny);
  }
  else{
    grid_bounds[1][1] = bounds[1][0] + (location[1] + 1)*ny*dx;
    pointBounds[1][1] = (location[1] + 1)*ny;
    sz[1] = ny + 1;
  }

  return buildGrid(sz, grid_bounds);
}
// }}}

// }}}
