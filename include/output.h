#ifndef OUTPUT_H
#define OUTPUT_H

namespace output{
  void outputCSV(const char *name, const double *v, const double *x, const double *y, const unsigned int shp[2], 
                 double time, unsigned int ghost, bool polar=false); 
};

#endif
