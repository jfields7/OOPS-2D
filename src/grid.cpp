#include <grid.h>
#include <iostream>
#include <new>
#include <cmath>

// Constructor {{{
Grid::Grid(const pair2<double>& bounds, unsigned int sz[2], unsigned int nghosts){
  // Make sure the bounds are ordered correctly. If not, swap them.
  if(bounds[0][0] < bounds[0][1]){
    grid_bounds[0][0] = bounds[0][0];
    grid_bounds[0][1] = bounds[0][1];
  }
  else{
    grid_bounds[0][0] = bounds[0][1];
    grid_bounds[0][1] = bounds[0][0];
  }
  if(bounds[1][0] < bounds[1][1]){
    grid_bounds[1][0] = bounds[1][0];
    grid_bounds[1][1] = bounds[1][1];
  }
  else{
    grid_bounds[1][0] = bounds[1][1];
    grid_bounds[1][1] = bounds[1][0];
  }

  // Let's calculate the grid spacing. Assume it's a cell-edge grid.
  dx = (grid_bounds[0][1] - grid_bounds[0][0])/(double)(sz[0]-1);
  if(fabs((grid_bounds[1][1] - grid_bounds[1][0])/(double)(sz[1]-1) - dx) > 1e-15){
    // Grid spacing is not equal. Complain.
    std::cout << "Grid spacing is not equal! Adjusting y-max bound.\n";
    grid_bounds[1][1] = grid_bounds[1][0] + dx*(sz[1] - 1);
  }

  // Construct the grid. Start by adding ghost points.
  shp[0] = sz[0] + 2*nghosts;
  shp[1] = sz[1] + 2*nghosts;
  try{
    points[0] = new double[shp[0]];
    points[1] = new double[shp[1]];
    for(unsigned int i = 0; i < shp[0]; i++){
      points[0][i] = grid_bounds[0][0] + ((int)i - (int)nghosts)*dx;
    }
    for(unsigned int i = 0; i < shp[1]; i++){
      points[1][i] = grid_bounds[1][0] + ((int)i - (int)nghosts)*dx;
    }
  }
  catch(std::bad_alloc& ba){
    std::cout << "Failed to allocate memory for grid.\n";
    shp[0] = 0;
    shp[1] = 0;
    dx = 0;
    points[0] = NULL;
    points[1] = NULL;
  }
}
// }}}

// Copy constructor {{{
Grid::Grid(const Grid& other){
  grid_bounds[0][0] = other.getBounds()[0][0];
  grid_bounds[0][1] = other.getBounds()[0][1];
  grid_bounds[1][0] = other.getBounds()[1][0];
  grid_bounds[1][1] = other.getBounds()[1][1];
  dx = other.getSpacing();
  shp[0] = other.getSize()[0];
  shp[1] = other.getSize()[1];
  try{
    points[0] = new double[shp[0]];
    points[1] = new double[shp[1]];
    for(unsigned int i = 0; i < shp[0]; i++){
      points[0][i] = other.getPoints()[0][i];
    }
    for(unsigned int i = 0; i < shp[1]; i++){
      points[1][i] = other.getPoints()[1][i];
    }
  }
  catch(std::bad_alloc& ba){
    std::cout << "Failed to allocate memory for grid.\n";
    shp[0] = 0;
    shp[1] = 0;
    dx = 0;
    points[0] = NULL;
    points[1] = NULL;
  }
  //std::cout << "Memory allocated.\n";
}
// }}}

// Destructor {{{
Grid::~Grid(){
  // Clear all the memory from the points.
    delete[] points[0];
    delete[] points[1];
}
// }}}
