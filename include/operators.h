#ifndef OPERATORS_H
#define OPERATORS_H

namespace operators{
  /**
   * A 2nd-order centered Laplacian operator.
   */
  inline double laplacian_2(const double x[3], const double y[3], const double dx, const double dy){
    return (x[0] - 2.0*x[1] + x[2])/(dx*dx) + (y[0] -2.0*y[1] + y[2])/(dy*dy);
  }

  /**
   * A 2nd-order centered first-derivative operator.
   */
  inline double dx_2(const double x[3], const double dx){
    return (x[2] - x[0])/(2.0*dx);
  }

  /**
   * A 2nd-order forward-difference first-derivative operator.
   */
  inline double dx_2off(const double x[3], const double dx){
    return (-3.0*x[0] + 4.0*x[1] - x[2])/(2.0*dx);
  }
}

#endif
