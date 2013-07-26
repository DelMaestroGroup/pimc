/**
 * @file pimc.h
 * @author Adrian Del Maestro
 * @date 10.27.2008
 *
 * @brief PathIntegralMonteCarlo class definition.
 */

#ifndef PIMC_H 
#define PIMC_H

#include "common.h"
#include "communicator.h"
#include "estimator.h"

class ActionBase;
class Path;
class LookupTable;
class MoveBase;
class EstimatorBase;

typedef boost::ptr_vector<EstimatorBase> estimator_vector;
typedef boost::ptr_vector<MoveBase> move_vector;

// ========================================================================  
// PathIntegralMonteCarlo Class
// ========================================================================  
/** 
 * The main 
 * driver class for the entire path integral monte carlo program.
 *
 * Holds the path, action, move and estimator objects along with methods
 * which actually perform the monte carlo sampling procedure.  
 */
class PathIntegralMonteCarlo {
	public:
		PathIntegralMonteCarlo (Path &, ActionBase *, MTRand &, 
                move_vector &, estimator_vector &, const bool);
		~PathIntegralMonteCarlo ();

		/* The equilibration relaxation */
		void equilStep(const uint32, const bool);

		/* The equilibration relaxation with off diagonal updates */
		void equilStepOffDiagonal(const uint32, const bool);

		/* The actual monte carlo step */
		void step();

		/* Output results of the simulation */
		void finalOutput();

		/* Save the PIMC state to disk */
		void saveState(); 

		/* Output the world-line configurations in the protein databank format */
		void outputPDB();

		int numStoredBins;			///< Number of stored estimators
		int numDiagonal;			///< Number of consecutive diagonal configs
		int numConfig;				///< Number of configurations;
		int numCoMAttempted;					///< Number of Center of Mass moves
		int numCoMAccepted;			///< Number of equil CoM moves accepted
		int numDisplaceAttempted;	///< Number of equil Displace moves
		int numDisplaceAccepted;    ///< Number of equil Displace moves accepted

	private:
		int configNumber;			// The output configuration number
        int numImagTimeSweeps;      // Partitioning used for updates
        int numSteps;               // The number of steps performed during equilibration
        int numUpdates;             // The maximum number of updates per step.
        int numParticles;           // The number of particles

        int comIndex;               // The index of the CoM move
        int diagIndex;              // The index of the diagonal path update (bisection or staging)
        int displaceIndex;          // The index of a possible displace move
        int obdmIndex;              // The index of the One Body Density Matrix estimator

		bool savedState;			// We have saved at least one state
		bool startWithState;		// Are we starting from a saved state
        bool success;               // Track move success/failure

		MTRand &random;				// The global random number generator
		Path &path;					// The actual wordline paths

        move_vector move;		    // The move array
        estimator_vector estimator;	// The estimator array

		vector <double> attemptDiagProb;		// The cumulative diagonal attempt Probabilities
		vector <double> attemptOffDiagProb;		// The cumulative off-diagonal attempt Probabilities

		/* Output estimators to disk */
		void output();

		/* Load the PIMC state from disk */
		void loadState();

		/* Load classical or quantum initial states */
        void loadClassicalState(Array <dVec,2>&, Array <unsigned int, 2>&, int);
        void loadQuantumState(Array <dVec,2>&, Array <beadLocator,2>&, Array<beadLocator,2>&, int, int);

		/* Shuffle the offDiag move list, returning it in moves*/
		void shuffleMoves();

		/* Perform a Metropolis update */
		string update(const double,const int);
};

#endif

