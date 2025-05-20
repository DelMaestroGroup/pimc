/**
 * @file cmc.h
 * @author Adrian Del Maestro
 * @date 08.13.2012
 *
 * @brief ClassicalMonteCarlo class definition
 */

#ifndef CMC_H 
#define CMC_H

#include "common.h"
#include "communicator.h"
#include "potential.h"

// ========================================================================  
// ClassicalMonteCarlo Class
// ========================================================================  
/** 
 * Pre-equilibration via classical Monte Carlo.
 *
 * We perform a series of classical Monte Carlo updates with the hopes of
 * obtaining a good classical ground state from which to initiate our quantum
 * monte carlo from.
 */
class ClassicalMonteCarlo{
    public:
        ClassicalMonteCarlo(PotentialBase *, PotentialBase *, MTRand &, const
                Container *, DynamicArray <dVec,1> &);
        ~ClassicalMonteCarlo();

        void run(uint32, bool);                  ///< Perform the Monte Carlo equilibration

    private:

        PotentialBase *externalPtr;     // The external potential
        PotentialBase *interactionPtr;  // The interaction potential
        MTRand &random;                 // A reference to the RNG
        const Container *boxPtr;        // A constant reference to the container class
        DynamicArray <dVec,1> &config;  // The particle configurations

        double z;                       // The fugacity
        double energy;                  // The total potential energy
        double deltaV;                  // Change in potential energy
        dVec sep;                       // A particle separation
        

        int numParticles;               // The current number of particles

        uint32 numMoveTotal;            // Number of moves attempted
        uint32 numMoveAccept;           // Number of moves accepted

        uint32 numInsertTotal;          // Number of inserts attempted
        uint32 numInsertAccept;         // Number of inserts accepted

        uint32 numDeleteTotal;          // Number of deletes attempted
        uint32 numDeleteAccept;         // Number of deletes accepted

        /* Measured quantities */
        double aveEnergy;
        double aveNumParticles;
        double aveEoN;

        /* Move, create and destroy particles */
        void moveParticle();
        void insertParticle();
        void deleteParticle();

        /* Perform measurements */
        void measure(int &);

        /* Compute the total potential energy */
        double getTotalEnergy();
};

#endif

