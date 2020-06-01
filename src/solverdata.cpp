#include <solverdata.h>
#include <new>
#include <iostream>

// Constructor {{{
SolverData::SolverData(unsigned int eqCount, unsigned int nStages, const Grid& grid) : ODEData(eqCount, grid){
  this->nStages = nStages;
  unsigned int shp[2] = {0};
  shp[0] = grid.getSize()[0];
  shp[1] = grid.getSize()[1];
  unsigned int n = shp[0]*shp[1];
  currStage = 0;
  // Try to allocate memory for the arrays.
  try{
    data_int = new double*[eqCount];
    for(unsigned int i = 0; i < eqCount; i++){
      data_int[i] = new double[n];
    }
    work = new double**[nStages];
    for(int i = 0; i < nStages; i++){
      work[i] = new double*[eqCount];
      for(int j = 0; j < eqCount; j++){
        work[i][j] = new double[n];
      }
    }
  }
  catch(std::bad_alloc& ba){
    std::cerr << "Failed to allocate memory for solver data.\n";
    nEq = 0;
    data_int = NULL;
    work = NULL;
    currStage = 0;
  }
}
// }}}

// Constructor {{{
SolverData::SolverData(const SolverData& other) : ODEData(other){
  std::cout << "SolverData copy constructor: This shouldn't be getting called, but it is.\n";
}
// }}}

// Destructor {{{
SolverData::~SolverData(){
  for(unsigned int i = 0; i < nStages; i++){
    for(unsigned int j = 0; j < nEq; j++){
      delete[] work[i][j];
    }
    delete[] work[i];
  }
  for(unsigned int i = 0; i < nEq; i++){
    delete[] data_int[i];
  }
  delete[] work;
  delete[] data_int;
}
// }}}
