#ifndef MPICOMMUNICATOR_H
#define MPICOMMUNICATOR_H

#include <mpi.h>
#include <string.h>
#include "types.h"

/**********************************************************************************
 * 
 * Class name: MPICommunicator
 * Author: Jacob Fields
 * Date Modified: 5-20-2020
 *
 * Description: The MPICommunicator class handles the most important details of
 *              setting up MPI to work on each processor. It is constructed as a
 *              singleton so that multiple instances cannot exist.
 *
 *********************************************************************************/
class MPICommunicator{
  private:
    static MPICommunicator comm;

    // We make construction and destruction completely private. We can't have multi
    MPICommunicator();
    ~MPICommunicator();
    MPICommunicator(const MPICommunicator&){};

    bool isInitialized;

    int worldSize;
    int worldRank;
    int rootRank;
    std::string processorName;
  public:
    /**
     * Get a reference to the MPICommunicator singleton instance.
     * @return An MPICommunicator reference.
     */
    static inline MPICommunicator* getInstance(){
      return &comm;
    }

    /**
     * Initialize MPI
     */
    Result init();
    /**
     * Cleanup MPI
     */
    Result cleanup();

    /**
     * Broadcast the parameters to the other processors.
     * @param par A pointer to the parameter to broadcast.
     * @return UNIINITIALIZED if uninitialized, FAILURE if an
     *         error occurred, and SUCCESS if no issues occurred.
     */
    Result broadcastParameter(int *par);
    Result broadcastParameter(unsigned int *par);
    Result broadcastParameter(double *par);
    Result broadcastParameter(std::string *par);

    /**
     * Send data to another processor.
     */
    //Result sendData(double *buffer, int size, int dest);
    /**
     * Receive data from a processor.
     */
    //Result getData(double *buffer, int size, int src);

    /**
     * Perform a sum operation on all processors. Rather than using
     * a reduce operation, it actually gathers all the data on the
     * root so it can perform compensated summation and then 
     * broadcasts it back to all the other processors.
     */
    Result allSum(double sendData, double &recvData);

    /**
     * Find the maximum across all processors.
     */
    Result findMax(double in, double &out);

    /**
     * Find the minimum across all processors.
     */
    Result findMin(double in, double &out);

    /**
     * Get the current MPI rank.
     * @return the rank of this processor.
     */
    inline int getRank() const{
      return worldRank;
    }

    /**
     * Get the MPI world size.
     * @return the number of processors currently in use.
     */
    inline int getWorldSize() const{
      return worldSize;
    }
    
    /**
     * Get the processor name.
     * @return a const string containing the processor name.
     */
    inline const std::string& getProcessorName() const{
      return processorName;
    }

    /**
     * Find out if MPI has been initialized.
     * @return false if uninitialized, true if initialized.
     */
    inline bool getIsInitialized() const {
      return isInitialized;
    }

    /**
     * Get the rank of the root process.
     * @return -1 if uninitialized, 0 otherwise.
     */
    inline int getRootRank(){
      return rootRank;
    }
};

#endif
