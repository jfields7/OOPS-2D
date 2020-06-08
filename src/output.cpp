#include <output.h>
#include <cstdio>
#include <cmath>
#include <base64.h>

void output::outputCSV(const char *name, const double *v, const double *x, const double *y, 
                       const unsigned int shp[2], double time, unsigned int ghost, bool polar){
  FILE *f;
  f = fopen(name, "w");
  if(!polar){
    for(unsigned int j = ghost; j < shp[1] - ghost; j++){
      for(unsigned int i = ghost; i < shp[0] - ghost; i++){
        unsigned int pp = i + shp[0]*j;
        fprintf(f,"%08g,   %08g,   %08g,   %08g\n", time, x[i], y[j], v[pp]);
      }
    }
  }
  else{
    for(unsigned int j = ghost; j < shp[1] - ghost; j++){
      for(unsigned int i = ghost; i < shp[0] - ghost; i++){
        unsigned int pp = i + shp[0]*j;
        fprintf(f,"%08g,   %08g,   %08g,   %08g\n", time, x[i]*cos(y[j]), x[i]*sin(y[j]), v[pp]);
      }
    }
  }

  fclose(f);
}

void output::writeVTKHeader(FILE *f, const unsigned int shp[2], const pair2<unsigned int>& bnds){
  unsigned int nx, ny, nxst, nxend, nyst, nyend;
  //nx = shp[0] - 2*ghost - 1;
  //ny = shp[1] - 2*ghost - 1;
  nx = shp[0] - 1;
  ny = shp[1] - 1;
  nxst = bnds[0][0];
  nxend = bnds[0][1];
  nyst = bnds[1][0];
  nyend = bnds[1][1];
  fprintf(f,"<?xml version=\"1.0\"?>\n");
  fprintf(f,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(f,"  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
    0, nx, 0, ny, 0, 0);
  fprintf(f,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",
    nxst, nxend, nyst, nyend, 0, 0);
}

void output::writeVTKTime(FILE *f, double time){
  fprintf(f,"      <FieldData>\n");
  fprintf(f,"        <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" format=\"binary\">");
  fprintf(f,"%s",Base64::encode(time).c_str());
  fprintf(f,"</DataArray>\n");
  //fprintf(f,"        <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" format=\"ascii\">\n");
  //fprintf(f,"          %g\n",time);
  //fprintf(f,"        </DataArray>\n");
  fprintf(f,"      </FieldData>\n");
}

void output::writeVTKScalar(FILE *f, const char *name, const double *v, const unsigned int shp[2],
                            unsigned int ghost){
  //fprintf(f,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"binary\">\n",name);
  fprintf(f,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",name);
  unsigned int pp;
  for(unsigned int j = ghost; j < shp[1] - ghost; j++){
    for(unsigned int i = ghost; i < shp[0] - ghost; i++){
      pp = i + shp[0]*j;
      fprintf(f,"%g ", v[pp]);
      //fprintf(f,"%s", Base64::encode(v[pp]).c_str());
    }
  }
  fprintf(f,"\n        </DataArray>\n");
  //fprintf(f,"</DataArray>\n");
}

void output::writeVTKPoints(FILE *f, const double *x, const double *y, const unsigned int shp[2],
                            unsigned int ghost, bool polar){
  double buffer[3] = {0.0};
  fprintf(f,"      <Points>\n");
  //fprintf(f,"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">");
  fprintf(f,"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  if(!polar){
    for(unsigned int j = ghost; j < shp[1] - ghost; j++){
      for(unsigned int i = ghost; i < shp[0] - ghost; i++){
        buffer[0] = x[i];
        buffer[1] = y[j];
        //fwrite(buffer, sizeof(double), sizeof(buffer), f);
        fprintf(f,"%g %g %g ",x[i], y[j], 0.0);
      }
    }
  }
  else{
    for(unsigned int j = ghost; j < shp[1] - ghost; j++){
      for(unsigned int i = ghost; i < shp[0] - ghost; i++){
        buffer[0] = x[i]*cos(y[j]);
        buffer[1] = x[i]*sin(y[j]);
        //fwrite(buffer, sizeof(double), sizeof(buffer), f);
        fprintf(f,"%g %g %g ",buffer[0], buffer[1], 0.0);
      }
    }
  }
  fprintf(f,"        </DataArray>\n");
  //fprintf(f,"<DataArray>\n");
  fprintf(f,"      </Points>\n");
}

void output::writeVTKFooter(FILE *f){
  fprintf(f,"    </Piece>\n");
  fprintf(f,"  </StructuredGrid>\n");
  fprintf(f,"</VTKFile>\n");
}

void output::writePVTKHeader(FILE *f, const unsigned int shp[2]){
  fprintf(f,"<?xml version=\"1.0\"?>\n");
  fprintf(f,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(f,"  <PStructuredGrid WholeExtent=\"%d %d %d %d %d %d\" GhostLevel=\"0\">\n",
    0, shp[0]-1, 0, shp[1]-1, 0, 0);
}

void output::writePVTKScalar(FILE *f, const char* name){
  fprintf(f,"      <PDataArray type=\"Float64\" Name=\"%s\"/>\n",name);
}

void output::writePVTKPoints(FILE *f){
  fprintf(f,"    <PPoints>\n");
  fprintf(f,"      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
  fprintf(f,"    </PPoints>\n");
}

void output::writePVTKPiece(FILE *f, const char* name, const pair2<unsigned int>& bnds){
  fprintf(f,"    <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>\n",
    bnds[0][0], bnds[0][1], bnds[1][0], bnds[1][1], 0, 0, name);
}

void output::writePVTKFooter(FILE *f){
  fprintf(f,"  </PStructuredGrid>\n");
  fprintf(f,"</VTKFile>\n");
}
