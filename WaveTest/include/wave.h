#ifndef WAVE_H
#define WAVE_H

#include <ode.h>
#include <waveparameters.h>

class Wave : public ODE{
  private:
    // Variable labels
    static const unsigned int U_PHI = 0;
    static const unsigned int U_PI  = 1;

    WaveParameters *params;

  protected:
    virtual void rhs(const Grid& grid, double **u, double **dudt);

  public:
    Wave(Domain* d, Solver* s);
    virtual ~Wave();

    virtual void initData();

    void setParameters(WaveParameters *p);
    WaveParameters* getParameters();
};
#endif
