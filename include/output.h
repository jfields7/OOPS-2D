#ifndef OUTPUT_H
#define OUTPUT_H
#include <cstdio>

namespace output{
  void outputCSV(const char *name, const double *v, const double *x, const double *y, const unsigned int shp[2], 
                 double time, unsigned int ghost, bool polar=false); 
  void writeVTKHeader(FILE *f, const double *x, const double *y, const unsigned int shp[2],
                      unsigned int ghost);
  void writeVTKTime(FILE *f, double time);
  void writeVTKScalar(FILE *f, const char* name, const double *v, const unsigned int shp[2],
                      const unsigned int ghost);
  void writeVTKPoints(FILE *f, const double *x, const double *y, const unsigned int shp[2],
                      unsigned int ghost, bool polar);
  void writeVTKFooter(FILE *f);
};

#endif
