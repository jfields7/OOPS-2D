#include <odedata.h>
#include <new>
#include <iostream>

// ODEData {{{
ODEData::ODEData(unsigned int eqCount, const Grid& grid, unsigned int lines) : mGrid(grid){
  nEq = eqCount;
  this->lines = lines;

  unsigned int shp[2] = {0};
  shp[0] = grid.getSize()[0];
  shp[1] = grid.getSize()[1];
  unsigned int n = shp[0]*shp[1];
  unsigned int nmax = (shp[0] > shp[1]) ? shp[0] : shp[1];

  // Try to allocate memory for the array.
  if(lines > 0){
    line = new double**[lines];
  }
  try{
    data = new double*[eqCount];
    for(unsigned int i = 0; i < lines; i++){
      line[i] = new double*[eqCount];
      for(unsigned int m = 0; m < eqCount; m++){
        line[i][m] = new double[nmax];
      }
    }
    for(unsigned int m = 0; m < eqCount; m++){
      data[m] = new double[n];
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for ODE data.\n";
    nEq = 0;
    data = NULL;
    line = NULL;
  }
}
// }}}

// Copy Constructor {{{
ODEData::ODEData(const ODEData& other) : mGrid(other.getGrid()){
  std::cout << "ODEData copy constructor: This shouldn't be getting called, but it is.\n";
}
// }}}

// ~ODEData {{{
ODEData::~ODEData(){
  for(unsigned int i = 0; i < nEq; i++){
    delete[] data[i];
  }
  for(unsigned int i = 0; i < lines; i++){
    for(unsigned int m = 0; m < nEq; m++){
      delete[] line[i][m];
    }
    delete[] line[i];
  }
  delete[] data;
}
// }}}
