#include <output.h>
#include <cstdio>

void output::outputCSV(const char *name, const double *v, const double *x, const double *y, 
                       const unsigned int shp[2], double time, unsigned int ghost){
  FILE *f;
  f = fopen(name, "w");
  for(unsigned int j = ghost; j < shp[1] - ghost; j++){
    for(unsigned int i = ghost; i < shp[0] - ghost; i++){
      unsigned int pp = i + shp[0]*j;
      fprintf(f,"%08g,   %08g,   %08g,   %08g\n", time, x[i], y[j], v[pp]);
    }
  }

  fclose(f);
}
