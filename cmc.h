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
 * Pre-equilibration via classical monte carlo.
 */
class ClassicalMonteCarlo{
	public:
        ClassicalMonteCarlo(PotentialBase *,PotentialBase *,MTRand &,const Container *, 
                Array <dVec,1> &);
        ~ClassicalMonteCarlo();

        void run();                        ///< Perform the Monte Carlo simulation

	private:

		PotentialBase *externalPtr;		// The external potential
		PotentialBase *interactionPtr;	// The interaction potential
        MTRand &random;			        // A reference to the RNG
		const Container *boxPtr;		// A constant reference to the container class
        Array <dVec,1> config;          // The particle configurations

        double z;                       // The fugacity
        double energy;                  // The total potential energy

        int numParticles;               // The current number of particles

        uint32 numUpdateTotal;          
        uint32 numUpdateAccept;         

        uint32 numInsertTotal;          
        uint32 numInsertAccept;         

        uint32 numDeleteTotal;          
        uint32 numDeleteAccept;         


        void insertMove();
        void deleteMove();
        void updateMove();
        void measure(int);

        double getTotalEnergy();

        double aveEnergy;
        double aveNumParticles;
        double aveEoN;
};

#endif

