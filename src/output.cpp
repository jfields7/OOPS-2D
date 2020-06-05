#include <output.h>
#include <cstdio>
#include <cmath>

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

void output::writeVTKHeader(FILE *f, const double *x, const double *y, const unsigned int shp[2],
                            unsigned int ghost){
  fprintf(f,"<?xml version=\"1.0\"?>\n");
  fprintf(f,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(f,"  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
    0, shp[0] - 2*ghost - 1, 0, shp[1] - 2*ghost - 1, 0, 0);
  fprintf(f,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",
    0, shp[0] - 2*ghost - 1, 0, shp[1] - 2*ghost - 1, 0, 0);
}

void output::writeVTKTime(FILE *f, double time){
  fprintf(f,"      <FieldData>\n");
  //fprintf(f,"        <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" format=\"binary\">");
  //fwrite(&time,sizeof(double),1,f);
  //fprintf(f,"</DataArray>\n");
  fprintf(f,"        <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" format=\"ascii\">\n");
  fprintf(f,"          %g\n",time);
  fprintf(f,"        </DataArray>\n");
  fprintf(f,"      </FieldData>\n");
}

void output::writeVTKScalar(FILE *f, const char *name, const double *v, const unsigned int shp[2],
                            unsigned int ghost){
  //fprintf(f,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"binary\">",name);
  fprintf(f,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",name);
  unsigned int pp;
  for(unsigned int j = ghost; j < shp[1] - ghost; j++){
    for(unsigned int i = ghost; i < shp[0] - ghost; i++){
      pp = i + shp[0]*j;
      //fwrite(&v[pp], sizeof(double), 1, f);
      fprintf(f,"%g ", v[pp]);
    }
  }
  fprintf(f,"        </DataArray>\n");
  //fprintf(f,"<DataArray>\n");
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
        fwrite(buffer, sizeof(double), sizeof(buffer), f);
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
