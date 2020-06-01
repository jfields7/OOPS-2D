#include <ode.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <mpicommunicator.h>
#include <output.h>
#include <sstream>
#include <fstream>
#include <stdio.h>

ODE::ODE(const unsigned int n, Domain *d, Solver *s) : nEqs(n), domain(d), solver(s){
  max_dx = 0.0;
  time = 0.0;
}

ODE::~ODE(){
}

// reallocateData {{{
Result ODE::reallocateData(){
  //fieldData = std::unique_ptr<FieldMap>(new FieldMap(domain->getGrid(), fieldList));
  fieldData = std::make_shared<FieldMap>(domain->getGrid(), fieldList);

  return SUCCESS;
}
// }}}

// addField {{{
Result ODE::addField(std::string name, unsigned int eqs, bool isEvolved, bool isComm){
  if(fieldList.find(name) != fieldList.end()){
    return FIELD_EXISTS;
  }
  else{
    unsigned int stages = (isEvolved ? solver->getNStages() : 0);
    fieldList.insert({name, FieldInfo(name, eqs, stages, isComm)});
  }

  return SUCCESS;
}
// }}}

// removeField {{{
Result ODE::removeField(std::string name){
  if(fieldList.find(name) == fieldList.end()){
    return UNRECOGNIZED_FIELD;
  }
  else{
    fieldList.erase(name);
  }
  return SUCCESS;
}
// }}}

// setDomain {{{
Result ODE::setDomain(Domain *d){
  domain = d;

  return reallocateData();
}
// }}}

// setSolver {{{
Result ODE::setSolver(Solver *s){
  solver = s;

  return reallocateData();
}
// }}}

// evolveStep {{{
Result ODE::evolveStep(double dt){
  // Keep track of the current time as well as the original time.
  double old_time = time;

  solver->setStageTime(old_time, time, dt, 0);
  // Copy the current solution into the intermediate variables.
  for(auto fieldpair : fieldData->getSolverFields()){
    auto field = fieldpair.second;
    field->setCurrentStage(0);
    double **data0 = field->getData();
    double **dataint = field->getIntermediateData();
    unsigned int vars = field->getEqCount();
    unsigned int shp = field->getGrid().getSize()[0]*field->getGrid().getSize()[1];
    for(unsigned int m = 0; m < vars; m++){
      for(unsigned int i = 0; i < shp; i++){
        dataint[m][i] = data0[m][i];
      }
    }
  }
  // Calculate the first stage.
  solver->calcStage(this, fieldData, dt, 0);

  // Perform the grid exchange and apply boundary conditions.
  doAfterStage();
  performGridExchange();
  doAfterExchange();
  applyBoundaries();
  doAfterBoundaries();

  for(unsigned int i = 1; i < solver->getNStages(); i++){
    solver->setStageTime(old_time, time, dt, i);
    // Loop over every data set in the domain to update the stage.
    for(auto fieldpair : fieldData->getSolverFields()){
      fieldpair.second->setCurrentStage(i);
    }
    solver->calcStage(this, fieldData, dt, i);

    // Perform the grid exchange and apply boundary conditions.
    doAfterStage();
    performGridExchange();
    doAfterExchange();
    applyBoundaries();
    doAfterBoundaries();
  }
  
  /*auto evol = fieldData->getSolverField("Evolution");
  // Loop over every stage for the solver.
  for(unsigned int i = 0; i < solver->getNStages(); i++){
    solver->setStageTime(old_time, time, dt, i);
    solver->calcStage(this,evol->getData(), evol->getIntermediateData(), (evol->getWorkData())[i],
                      evol->getGrid(), dt, i);

    doAfterStage(true);
    performGridExchange(true);
    doAfterExchange(true);
    applyBoundaries(true);
    doAfterBoundaries(true);
  }*/

  time = old_time + dt;
  //solver->combineStages(evol->getWorkData(), evol->getData(), evol->getGrid(), dt, nEqs);
  solver->combineStages(fieldData, dt);
  doAfterStage();
  performGridExchange();
  doAfterExchange();
  applyBoundaries();
  doAfterBoundaries();

  // Swap all the intermediate data back into the solution array.
  for(auto fieldpair : fieldData->getSolverFields()){
    auto field = fieldpair.second;
    double **data0 = field->getData();
    double **dataint = field->getIntermediateData();
    unsigned int vars = field->getEqCount();
    unsigned int shp = field->getGrid().getSize()[0]*field->getGrid().getSize()[1];
    for(unsigned int m = 0; m < vars; m++){
      for(unsigned int i = 0; i < shp; i++){
        data0[m][i] = dataint[m][i];
      }
    }
  }

  return SUCCESS;
}
// }}}

// performGridExchange {{{
void ODE::performGridExchange(){
  // Collect some information before getting started.
  MPICommunicator *comm = MPICommunicator::getInstance();
  int rank = comm->getRank();
  int size = comm->getWorldSize();

  // If there's only one rank, there's no reason to exchange data.
  if(size == 1){
    return;
  }

  // Some information to make this easier to generalize later.
  const unsigned int dim = 2;
  const unsigned int neighbors = 1 << dim;

  // Count up the number of fields we need to transfer.
  unsigned int nfields = 0;
  for(auto field : fieldList){
    FieldInfo &info = field.second;
    if(info.isComm){
      nfields += info.nEqs;
    }
  }
  
  // First, let's create data buffers for our fields.
  const Grid* grid = domain->getGrid();
  pair2<int> commPartners = domain->getCommPartners();
  unsigned int shp[dim];
  for(unsigned int i = 0; i < dim; i++){
    shp[i] = grid->getSize()[i];
  }
  unsigned int nb = domain->getGhostPoints();
  double *sendBuffer[neighbors];
  double *recvBuffer[neighbors];
  MPI_Request request_recv[neighbors], request_send[neighbors];
  MPI_Status status[neighbors];
  // FIXME: This unholy construction here has five loops.
  // If there's not a better way to do this, we really need
  // to move part of it into a new function.
  for(unsigned int d = 0; d < dim; d++){
    // Figure out how much memory to allocate.
    unsigned int size = 1;
    for(unsigned int k = 0; k < dim; k++){
      if(k != d){
        size *= shp[k];
      }
    }
    for(unsigned int b = 0; b < 2; b++){
      // Get the buffer index.
      unsigned int pos = b + 2*d;
      // Allocate the buffers.
      sendBuffer[pos] = new double[nfields*nb*size];
      recvBuffer[pos] = new double[nfields*nb*size];
    }
    // Copy the physical boundaries into the send buffers.
    unsigned int curfield = 0;
    for(auto field : fieldList){
      FieldInfo &info = field.second;
      if(!info.isComm){
        continue;
      }
      double** data;
      if(info.nStages == 0){
        data = (*fieldData)[field.first]->getData();
      }
      else{
        data = fieldData->getSolverField(field.first)->getIntermediateData();
      }
      unsigned int eq = info.nEqs;
      for(unsigned int m = 0; m < eq; m++){
        unsigned int offset = (curfield + m)*nb*size;
        if(commPartners[d][0] != -1){
          for(unsigned int j = 0; j < size; j++){
            for(unsigned int i = nb+1; i < 2*nb+1; i++){
              unsigned int pp = grid->getIndex(i, j, d);
              unsigned int pb = (i - nb - 1) + nb*j;
              sendBuffer[2*d][offset + pb] = data[m][pp];
            }
          }
        }
        if(commPartners[d][1] != -1){
          for(unsigned int j = 0; j < size; j++){
            for(unsigned int i = shp[d] - 2*nb - 1; i < shp[d] - nb; i++){
              unsigned int pp = grid->getIndex(i, j, d);
              unsigned int pb = (i - shp[d] + 2*nb + 1) + nb*j;
              sendBuffer[1 + 2*d][offset + pb] = data[m][pp];
            }
          }
        }
      }
      curfield += info.nEqs;
    }
  }

  // Retrieve the physical points from a neighboring processor which overlap
  // this processor's ghost region. Send physical points from this processor.
  for(unsigned int d = 0; d < dim; d++){
    unsigned int size = 1;
    for(unsigned int k = 0; k < dim; k++){
      if(k != d){
        size *= shp[k];
      }
    }
    unsigned int length = nfields*nb*size;
    for(unsigned int b = 0; b < 2; b++){
      unsigned int pos = b + 2*d;
      if(commPartners[d][b] != -1){
        MPI_Irecv(recvBuffer[pos], length, MPI_DOUBLE, commPartners[d][b], 
        commPartners[d][b], MPI_COMM_WORLD, &request_recv[pos]);
        MPI_Isend(sendBuffer[pos], length, MPI_DOUBLE, commPartners[d][b], rank, MPI_COMM_WORLD, 
                  &request_send[pos]);
      }
    }
  }

  for(unsigned int d = 0; d < dim; d++){
    for(unsigned int b = 0; b < 2; b++){
      if(commPartners[d][b] != -1){
        unsigned int p = b + 2*d;
        MPI_Wait(&request_send[p], &status[p]);
        MPI_Wait(&request_recv[p], &status[p]);
      }
    }
  }

  // Move the received data back into the right dataset.
  for(unsigned int d = 0; d < dim; d++){
    // Figure out how much data to copy.
    unsigned int size = 1;
    for(unsigned int k = 0; k < dim; k++){
      if(k != d){
        size *= shp[k];
      }
    }
    unsigned int curfield = 0;
    for(auto field : fieldList){
      FieldInfo &info = field.second;
      if(!info.isComm){
        continue;
      }
      double **data;
      if(info.nStages == 0){
        data = (*fieldData)[field.first]->getData();
      }
      else{
        data = fieldData->getSolverField(field.first)->getIntermediateData();
      }
      unsigned int eq = info.nEqs;
      for(unsigned int m = 0; m < eq; m++){
        unsigned int offset = (curfield + m)*nb*size;

        if(commPartners[d][0] != -1){
          for(unsigned int j = 0; j < size; j++){
            for(unsigned int i = 0; i < nb; i++){
              unsigned int pp = grid->getIndex(i, j, d);
              unsigned int pb = i + nb*j;
              data[m][pp] = recvBuffer[2*d][offset + pb];
            }
          }
        }
        if(commPartners[d][1] != -1){
          for(unsigned int j = 0; j < size; j++){
            for(unsigned int i = shp[d]-nb; i < shp[d]; i++){
              unsigned int pp = grid->getIndex(i, j, d);
              unsigned int pb = (i - shp[d] + nb) + nb*j;
              data[m][pp] = recvBuffer[1 + 2*d][offset + pb];
            }
          }
        }
      }
      curfield += eq;
    }
  }

  // Delete memory we just allocated.
  for(unsigned int i = 0; i < neighbors; i++){
    delete[] sendBuffer[i];
    delete[] recvBuffer[i];
  }
}
// }}}

// dumpField {{{
void ODE::dumpField(std::string field, char* name, double t, unsigned int var){
  unsigned int nb = domain->getGhostPoints();
  const Grid *grid = domain->getGrid();
  double **data = (*fieldData)[field]->getData();
  auto points = grid->getPoints();

  unsigned int shp[2] = {0};
  shp[0] = grid->getSize()[0];
  shp[1] = grid->getSize()[1];

  std::stringstream ss;
  ss << name;
  MPICommunicator *comm = MPICommunicator::getInstance();
  int rank = comm->getRank();
  int root = comm->getRootRank();
  int size = comm->getWorldSize();
  if(rank != root){
    ss << "." << rank;
  }
  char newname[256];
  strcpy(newname,ss.str().c_str());

  output::outputCSV(newname, data[var], points[0], points[1], shp, t, nb);

  // Collate the output.
  if(rank == root && size > 1){
    MPI_Request *request = new MPI_Request[size-1];
    MPI_Status *status = new MPI_Status[size-1];
    int *done = new int[size-1];
    // Open the primary data file.
    std::ofstream out(name, std::ios::app);
    if(!out.is_open()){
      std::cout << "There was an error while saving data.";
      return;
    }
    // Request a status from all the processors.
    for(int i = 1; i < size; i++){
      MPI_Irecv(&done[i-1], 1, MPI_INT, i, i, MPI_COMM_WORLD, &request[i-1]);
    }
    // Append all the new ranks to it.
    for(int i = 1; i < size; i++){
      // Wait for the process to finish before collating its output.
      MPI_Wait(&request[i-1], &status[i-1]);
      std::stringstream file;
      file << name << "." << i;
      std::ifstream in(file.str());
      if(!in.is_open()){
        std::cout << "There was an error while saving data.";
        return;
      }
      out << in.rdbuf();
      // Close the stream. Delete the temporary file.
      in.close();
      remove(file.str().c_str());
    }
    out.close();
    delete[] done;
    delete[] status;
    delete[] request;
  }
  else if (size > 1){
    // Send a message saying data has been sent.
    MPI_Request request;
    MPI_Status status;
    int done = 1;
    MPI_Isend(&done, 1, MPI_INT, root, rank, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);
  }
}
// }}}
