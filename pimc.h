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

class Path;
class LookupTable;
class MoveBase;
class EstimatorBase;

/** A vector containing measurable estimators */
typedef boost::ptr_vector<EstimatorBase> estimator_vector;

/** A vector containing Monte Carlo updates */
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
        PathIntegralMonteCarlo (boost::ptr_vector<Path> &,MTRand &, boost::ptr_vector<move_vector> &,
                                boost::ptr_vector<estimator_vector> &, const bool, uint32 binSize=100);
		~PathIntegralMonteCarlo ();

		/* The equilibration relaxation */
		void equilStep(const uint32, const bool, const bool);

		/* The equilibration relaxation with off diagonal updates */
		void equilStepOffDiagonal(const uint32, const bool);

		/* The actual monte carlo step */
		void step();

		/* Output results of the simulation */
		void finalOutput();

		/* Save the PIMC state to disk */
		void saveState(const int finalSave = 0);

		/* Output the world-line configurations in the protein databank format */
		void outputPDB();

		int numStoredBins;			///< Number of stored estimators
		int numDiagonal;			///< Number of consecutive diagonal configs
		int numConfig;				///< Number of configurations;
		int numCoMAttempted;		///< Number of Center of Mass moves
		int prevNumCoMAttempted;	///< Number of Center of Mass moves attempted
		int numCoMAccepted;			///< Number of equil CoM moves accepted
		int numDisplaceAttempted;	///< Number of equil Displace moves
		int numDisplaceAccepted;    ///< Number of equil Displace moves accepted
        int numMuAttempted;         ///< Number of moves between mu adjustments
        int numNAttempted;          ///< The number of particle measurements
        int cN;                     ///< The current total number of particles

        void printWormState();

	private:
		MTRand &random;				// The global random number generator

		int configNumber;			// The output configuration number
        int numImagTimeSweeps;      // Partitioning used for updates
        int numSteps;               // The number of steps performed during equilibration
        int numUpdates;             // The maximum number of updates per step.
        uint32 binSize;                // The maximum number measurements per bin.
        int numParticles;           // The number of particles

        uint32 Npaths;                 // Number of paths
    
        vector<string> stateStrings; // A vector of state strings from the last bin

		bool startWithState;		// Are we starting from a saved state
        bool success;               // Track move success/failure
    
        boost::ptr_vector<Path> &pathPtrVec; // The vector of all paths
        Path &path;                          // A reference to the first path in the path vector

        boost::ptr_vector<move_vector> &movePtrVec;  // Vector of all move_vectors
        move_vector &move;                           // A reference to the first move_vector
    
        boost::ptr_vector<estimator_vector> &estimatorPtrVec;  // Vector of estimators for 
                                                               // each path followed by multipath
        estimator_vector &estimator;                           // A reference to the first estimator_vector

		vector <double> attemptDiagProb;		// The cumulative diagonal attempt Probabilities
		vector <double> attemptOffDiagProb;		// The cumulative off-diagonal attempt Probabilities

		map <string,int> moveIndex;	            // A map to keep track of move names and indices
		map <string,int> estimatorIndex;	    // A map to keep track of estimator names and indices

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
		string update(const double,const int,const int);

};

#endif

