#ifndef GRID_H
#define GRID_H

#include "types.h"
#include <array>

/******************************************************************************
 *
 * Class name: Grid
 * Author: Jacob Fields
 * Date Modified: 5-21-2020
 * 
 * Description: This defines a 2d grid with fixed spacing.
 *
 * Usage: The grid has no default constructor. Construct a new grid with
 *        Grid(const double bounds[2][2], unsigned int sz[2], unsigned nghosts)
 *
 *****************************************************************************/

class Grid{
  private:
    // Grid spacing.
    double dx;
    // Number of grid points.
    unsigned int shp[2];
    // The boundaries of the grid.
    pair2<double> grid_bounds;
    // The grid points themselves.
    double *points[2];
  public:
    /**
     * A grid constructor.
     * @param bounds A pair of ordered pairs describing the left and rightmost
     *   bounds of the grid. This should not include ghost points.
     * @param sz A pair (nx, ny) describing the number of points on the grid.
     * @param nghosts The number of ghost points on the grid.
     */
    Grid(const pair2<double> &bounds, unsigned int sz[2], unsigned nghosts);

    /**
     * A copy constructor. Because of the memory allocated for our points,
     * we need to define a copy constructor that performs a deep copy.
     */
    Grid(const Grid& other);

    /**
     * A grid destructor. It frees any memory used.
     */
    ~Grid();

    /**
     * Get the data points on the grid, including ghost points.
     * @returns a const double** containing the x-coordinates and y-coordinates
     *   of the grid points.
     */
    inline const double* const* getPoints() const{
      return points;
    }

    /**
     * Get the number of grid points. This does include ghost points.
     * @returns A const int[2] containing the number of grid points.
     */
    inline const unsigned int* getSize() const{
      return shp;
    }

    /**
     * Get the bounds (without ghost points) of this grid.
     * @returns A pair of ordered pairs specifiying the left and right bounds of
     *   the x and y directions, respectively.
     */
    inline const pair2<double>& getBounds() const{
      return grid_bounds;
    }

    /**
     * Get the grid spacing for this grid. Because it is assumed to be a uniform grid,
     * This grid spacing applies to both the x and y directions.
     */
    inline double getSpacing() const{
      return dx;
    }
    
    /**
     * Calculate a Fortran-style array index based on this grid.
     */
    inline unsigned int getIndex(unsigned int i, unsigned int j) const{
      return i + shp[0]*j;
    }
    
    /**
     * Calculate a Fortran-style array index based on a particular direction.
     */
    inline unsigned int getIndex(unsigned int i, unsigned int j, unsigned int d) const{
      switch(d){
        case 0:
          return i + shp[0]*j;
          break;
        case 1:
          return j + shp[0]*i;
          break;
        default:
          return 0;
          break;
      }
    }

};

#endif
