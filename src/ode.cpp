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
Result ODE::addField(std::string name, unsigned int eqs, bool isEvolved, bool isComm, unsigned int lines){
  if(fieldList.find(name) != fieldList.end()){
    return FIELD_EXISTS;
  }
  else{
    unsigned int stages = (isEvolved ? solver->getNStages() : 0);
    fieldList.insert({name, FieldInfo(name, eqs, stages, lines, isComm)});
    fieldOutput.insert({name,std::vector<varPair>()});
    std::stringstream ss;
    for(unsigned int i = 0; i < eqs; i++){
      ss.str(std::string(""));
      ss << name << "::var" << i;
      fieldOutput[name].push_back(std::make_pair(false,ss.str()));
    }
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
    fieldOutput.erase(name);
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
  if(size == 1 && !domain->isPeriodic()){
    return;
  }
  else if(size == 1 && domain->isPeriodic()){
    // We need to make sure that the periodic theta boundary gets
    // applied, though.
    performPeriodicExchange();
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

// performPeriodicExchange {{{
void ODE::performPeriodicExchange(){
  const Grid* grid = domain->getGrid();
  auto shp = grid->getSize();
  double **data;
  unsigned int nb = domain->getGhostPoints();
  for(auto field : fieldList){
    FieldInfo &info = field.second;
    if(!info.isComm){
      continue;
    }
    if(info.nStages == 0){
      data = (*fieldData)[field.first]->getData();
    }
    else{
      data = fieldData->getSolverField(field.first)->getIntermediateData();
    }
    unsigned int eq = info.nEqs;
    for(unsigned int m = 0; m < eq; m++){
      for(unsigned int j = 0; j < nb; j++){
        for(unsigned int i = 0; i < shp[0]; i++){
          unsigned int ps = grid->getIndex(i, j);
          unsigned int pc = grid->getIndex(i, j + shp[1] - 2*nb - 1);
          data[m][ps] = data[m][pc];

          ps = grid->getIndex(i, shp[1] - nb + j);
          pc = grid->getIndex(i, nb + j + 1);
          data[m][ps] = data[m][pc];
        }
      }
    }
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

  output::outputCSV(newname, data[var], points[0], points[1], shp, t, nb,domain->getCoordinates() == Domain::POLAR);

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

// outputVTK {{{
void ODE::outputVTK(char* name, double t){
  MPICommunicator *comm = MPICommunicator::getInstance();
  FILE *f;
  int rank = comm->getRank();
  int size = comm->getWorldSize();
  int root = comm->getRootRank();
  // If we have more than one rank, we need multiple files.
  std::stringstream ss;
  ss << name;
  ss << ".vts";
  if(size > 1){
    ss << "." << rank;
  }
  char newname[256];
  strcpy(newname, ss.str().c_str());
  f = fopen(newname, "w");
  
  unsigned int nb = domain->getGhostPoints();
  const Grid *grid = domain->getGrid();
  auto points = grid->getPoints();

  unsigned int domsz[2] = {0};
  domsz[0] = domain->getSize()[0];
  domsz[1] = domain->getSize()[1];
  unsigned int shp[2] = {0};
  shp[0] = grid->getSize()[0];
  shp[1] = grid->getSize()[1];

  bool polar = domain->getCoordinates() == Domain::POLAR;

  output::writeVTKHeader(f, domsz, domain->getPointBounds());
  output::writeVTKTime(f, t);
  ss.str(std::string(""));
  char varnames[256];
  for(auto field : fieldList){
    const FieldInfo& info = field.second;
    for(unsigned int i = 0; i < info.nEqs; i++){
      const varPair& v = fieldOutput[field.first][i];
      if(v.first == true){
        ss << v.second.c_str() << " ";
      }
    }
  }
  strcpy(varnames, ss.str().c_str());
  fprintf(f,"      <PointData Scalars=\"%s\">\n",varnames);
  for(auto field : fieldList){
    const FieldInfo& info = field.second;
    double **data = (*fieldData)[field.first]->getData();
    for(unsigned int i = 0; i < info.nEqs; i++){
      const varPair& v = fieldOutput[field.first][i];
      if(v.first == true){
        output::writeVTKScalar(f, v.second.c_str(), data[i], shp, nb);
      }
    }
  }
  fprintf(f,"      </PointData>\n");
  output::writeVTKPoints(f, points[0], points[1], shp, nb, polar);
  output::writeVTKFooter(f);

  fclose(f);

  // If we have multiple ranks, we need a PVTS file as well.
  if(size > 1){
    if(rank == root){
      MPI_Request *request = new MPI_Request[size];
      MPI_Status *status = new MPI_Status[size];
      unsigned int data[4] = {0};
      pair2<unsigned int> bnds;
      ss.str(std::string());
      ss << name;
      ss << ".pvts";

      strcpy(newname, ss.str().c_str());
      f = fopen(newname, "w");

      output::writePVTKHeader(f, domsz);
      fprintf(f,"    <PPointData Scalars=\"%s\">\n",varnames);
      for(auto field : fieldList){
        const FieldInfo& info = field.second;
        for(unsigned int i = 0; i < info.nEqs; i++){
          const varPair& v = fieldOutput[field.first][i];
          if(v.first == true){
            output::writePVTKScalar(f, v.second.c_str());
          }
        }
      }
      fprintf(f,"    </PPointData>\n");
      output::writePVTKPoints(f);
      for(int i = 0; i < size; i++){
        ss.str(std::string());
        ss << name;
        ss << ".vts." << i;
        strcpy(newname, ss.str().c_str());
        if(i == rank){
          output::writePVTKPiece(f, newname, domain->getPointBounds());
        }
        else{
          // Get the point bounds from the other processors.
          MPI_Irecv(data, 4, MPI_UNSIGNED, i, i, MPI_COMM_WORLD, &request[i]);
          MPI_Wait(&request[i],&status[i]);
          bnds[0][0] = data[0];
          bnds[0][1] = data[1];
          bnds[1][0] = data[2];
          bnds[1][1] = data[3];
          output::writePVTKPiece(f, newname, bnds);
        }
      }
      output::writePVTKFooter(f);

      fclose(f);
      delete[] status;
      delete[] request;
    }
    else{
      unsigned int data[4] = {0};
      auto bnds = domain->getPointBounds();
      data[0] = bnds[0][0];
      data[1] = bnds[0][1];
      data[2] = bnds[1][0];
      data[3] = bnds[1][1];
      MPI_Request request;
      MPI_Status status;
      MPI_Isend(data, 4, MPI_UNSIGNED, root, rank, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
    }
  }
}
// }}}
