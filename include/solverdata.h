#ifndef SOLVERDATA_H
#define SOLVERDATA_H

#include "types.h"
#include "grid.h"
#include "odedata.h"

class SolverData : public ODEData{
  private:
    // Private members {{{
    /**
     * The number of stages used by the solver.
     */
    unsigned int nStages;

    /**
     * The intermediate date for the solution.
     */
    double **data_int;

    /**
     * The storage for the results from the ODE solver.
     */
    double ***work;
    // }}}

    /**
     * A copy constructor for SolverData. We make this private so that
     * the object cannot be copied.
     */
    SolverData(const SolverData& other);

  public:
  /**
   * The constructor for SolverData.
   * @param eqCount - the number of equations in the system.
   * @param stages  - the number of stages used by the solver.
   * @param grid    - the Grid this data belongs to.
   */
  SolverData(unsigned int eqCount, unsigned int stages, const Grid& grid);

  /**
   * A destructor for SolverData that deletes all the memory that was allocated.
   */
  virtual ~SolverData();

  // Inline getter functions {{{
  /**
   * Get the array containing the intermediate data. This is stored as [vars][points].
   */
  inline double** getIntermediateData() const{
    return data_int;
  }

  /**
   * Get the array containing the work data. This is stored as [stage][vars][points].
   */
  inline double*** getWorkData() const{
    return work;
  }
  // }}}
};

#endif
