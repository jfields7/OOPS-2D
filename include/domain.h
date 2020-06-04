#ifndef DOMAIN_H
#define DOMAIN_H

#include "types.h"
#include "grid.h"

/*****************************************************************************
 *
 * Class name: Domain
 * Author: Jacob Fields
 * Date Modified: 5-20-2020
 * 
 * Description: This defines the 2d domain used for the solver. It should
 *              hold onto the primary grid in use for this processor.
 *
 ****************************************************************************/

class Domain{
  public:
    enum GridCoordinates{
      CARTESIAN,
      POLAR
    };
  private:
    /**
     * The bounds of the problem.
     */
    pair2<double> bounds;

    /**
     * The CFL to use on the grid.
     */
    double d_cfl;

    /**
     * The number of ghost points permitted on the grids.
     */
    unsigned int nghosts;

    /**
     * The Grid associated with this particular MPI rank.
     */
    Grid *grid;

    /**
     * This processor's physical location on the domain.
     */
    int location[2];

    /**
     * This processor's boundaries in point indices.
     */
    pair2<int> pointBounds;

    /**
     * The communication partners for this processor. An entry of -1 indicates
     * the processor has no partner and actually contains the physical boundary.
     */
    pair2<int> commPartners;

    /**
     * The coordinate system for Grid objects on this Domain.
     */
    GridCoordinates coords;

    /**
     * Whether or not the theta variable is a periodic boundary.
     */
    bool periodic;


    /**
     * Boundary flags
     */
    unsigned int bflag;

    // MPI Routines {{{
    void assignCommunicationPartners(int dims[2]);

    Result divideGrids(unsigned int shp[2], int dims[2]);
    // }}}

  public:
    /**
     * The default constructor. It sets up the domain with a set of default
     *   parameters, i.e., a region from [0,0] to [1,1] with 3 ghost points
     *   and a CFL of 0.5.
     */
    Domain();

    /**
     * The destructor.
     */
    ~Domain();


    // Getters and Setters {{{
    /**
     * Set the CFL for this domain.
     * @param cfl The CFL to use.
     */
    inline void setCFL(double cfl){
      d_cfl = cfl;
    }

    /**
     * Get the CFL for this domain.
     * @returns A double indicating the CFL.
     */
    inline double getCFL() const{
      return d_cfl;
    }

    /**
     * Set the bounds for this domain.
     * WARNING: This will not perform bounds-checking on the grid.
     * @param bounds A pair of ordered pairs describing the bounds.
     */
    void setBounds(const pair2<double> &bounds);

    /**
     * Get the bounds for this domain.
     * @returns A pair of ordered pairs specifying the bounds.
     */
    inline const pair2<double>& getBounds() const{
      return bounds;
    }

    /**
     * Set the number of ghost points to use on this domain. This will automatically
     * reallocate the grid.
     * @param n The number of ghost points.
     */
    inline void setGhostPoints(unsigned int n){
      nghosts = n;
    }

    /**
     * Get the number of ghost points used by the domain.
     * @returns an unsigned integer indicating the number of ghost points.
     */
    inline unsigned int getGhostPoints() const{
      return nghosts;
    }

    /**
     * Get a reference to the Grid on this domain.
     * @returns a const reference to a Grid object.
     */
    inline const Grid* getGrid() const{
      return grid;
    }

    /**
     * Get the communication partners for this particular rank.
     * @returns a pair2<int> containing the communication partners.
     */
    inline const pair2<int>& getCommPartners() const{
      return commPartners;
    }

    /**
     * Check if the Domain on this processor has a physical boundary.
     */
    inline bool hasBoundary(Boundary b) const{
      return (bflag & (1 << b)) > 0;
    }

    /**
     * Check if this Domain's Grid is in polar or cartesian coordinates.
     */
    inline GridCoordinates getCoordinates() const{
      return coords;
    }

    /**
     * Check if this Domain has a periodic polar boundary.
     */
    inline bool isPeriodic() const{
      return periodic;
    }
    // }}}

    /**
     * Set up the grid on the domain with a different number of points
     * and the specified bounds.
     * @param shp A pair of integers specifying the x and y dimensions,
     *            respectively.
     * @param bounds A pair of ordered pairs describing the bounds for this
     *               Grid.
     * @returns UNEQUAL_SPACING if the resulting grid increments are not equal
     *                          in the x and y directions. Otherwise, SUCCESS.
     */
    Result buildGrid(unsigned int shp[2], pair2<double> &bounds);

    // MPI Routines {{{
    /**
     * Build a mesh with the specified size to cover the entire domain. Divide it
     * up across all the processors.
     * @param shp A pair of integers specifying the number of points in the
     *            x and y dimensions, respectively.
     * @returns SUCCESS if the mesh was build successfully, otherwise it returns
     *          an error code.
     */
    Result buildMesh(unsigned int shp[2], GridCoordinates coord=CARTESIAN, bool periodic=false);
    // }}}
};

#endif
