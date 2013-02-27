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
#include "move.h"
#include "estimator.h"

class ActionBase;
class Path;
class LookupTable;

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
		PathIntegralMonteCarlo (Path &, ActionBase *, MTRand &, const double, const bool);
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
		int numCoM;					///< Number of Center of Mass moves
		int numCoMAccepted;			///< Number of equil CoM moves accepted

	private:
		int configNumber;			// The output configuration number
        int numImagTimeSweeps;      // Partitioning used for updates
        int numSteps;               // The number of steps performed during equilibration
        int numUpdates;             // The maximum number of updates per step.
        int numParticles;           // The number of particles
		bool savedState;			// We have saved at least one state
		bool fixedNumLevels;		// Have we fixed the number of levels
		bool startWithState;		// Are we starting from a saved state
        bool success;               // Track move success/failure

		MTRand &random;				// The global random number generator
		Path &path;					// The actual wordline paths

		/* The actual moves to be performed */
		CenterOfMassMove centerOfMass;
		StagingMove    	 staging;
		OpenMove         open;
		CloseMove        close;
		InsertMove       insert;
		RemoveMove       remove;
		AdvanceHeadMove	 advanceHead;
		RecedeHeadMove 	 recedeHead;
		AdvanceTailMove	 advanceTail;
		RecedeTailMove 	 recedeTail;
		SwapHeadMove 	 swapHead;
		SwapTailMove 	 swapTail;

		/* The estimators that we will measure */
		DiagonalFractionEstimator       diagonalFraction;
		EnergyEstimator                 energy;
		NumberParticlesEstimator        numberParticles;
        ParticlePositionEstimator       particlePositions;
        PlaneParticlePositionEstimator  planeParticlePositions;
		SuperfluidFractionEstimator     superfluidFraction;
        LocalSuperfluidDensityEstimator localSuperfluidDensity;
        PlaneWindingSuperfluidDensityEstimator   planeWindingSuperfluidDensity;
        PlaneAreaSuperfluidDensityEstimator   planeAreaSuperfluidDensity;
        RadialWindingSuperfluidDensityEstimator  radialWindingSuperfluidDensity;
        RadialAreaSuperfluidDensityEstimator     radialAreaSuperfluidDensity;
		PermutationCycleEstimator       permutationCycle;
		OneBodyDensityMatrixEstimator   oneBodyDensityMatrix;
		PairCorrelationEstimator        pairCorrelation;
		RadialDensityEstimator          radialDensity;
		WormPropertiesEstimator         wormProperties;
		NumberDistributionEstimator     numberDistribution;

		/* The cylindrical estimators */
		CylinderEnergyEstimator               cylEnergy;
		CylinderNumberParticlesEstimator      cylNumberParticles;
		CylinderSuperfluidFractionEstimator   cylSuperFluidFraction;
		CylinderPairCorrelationEstimator      cylPairCorrelation;
		CylinderOneBodyDensityMatrixEstimator cylOneBodyDensityMatrix;
		CylinderNumberDistributionEstimator   cylNumberDistribution;
		CylinderRadialPotentialEstimator      cylRadialPotential;

		vector <EstimatorBase*> estimator;		// The estimator base array
		vector <MoveBase*> move;				// The move base array

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

		/* Perform all moves */
		string runMoves(const double,const int);
};

#endif

