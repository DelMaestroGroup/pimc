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

/** A std::vector containing measurable estimators */
typedef boost::ptr_vector<EstimatorBase> estimator_vector;

/** A std::vector containing Monte Carlo updates */
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
                                boost::ptr_vector<estimator_vector> &, const bool);
        ~PathIntegralMonteCarlo ();


        /* The equilibration relaxation */
        void equilStepDiagonal();

        bool  equilStepRelaxmu();
        bool  equilStepRelaxC0();

        void equilStep(const uint32, const bool, const bool);

        /* The actual monte carlo step */
        void step();

        /* Output results of the simulation */
        void finalOutput();

        /* Save the PIMC state to disk */
        void saveState(const int finalSave = 0);
        void saveStateHDF5(const int finalSave = 0);

        /* Output the world-line configurations in the protein databank format */
        void outputPDB();

        int numStoredBins;          ///< Number of stored estimators
        int numDiagonal;            ///< Number of consecutive diagonal configs
        int numConfig;              ///< Number of configurations;
        int numCoMAttempted;        ///< Number of Center of Mass moves
        int prevNumCoMAttempted;    ///< Number of Center of Mass moves attempted
        int numCoMAccepted;         ///< Number of equil CoM moves accepted
        int numDisplaceAttempted;   ///< Number of equil Displace moves
        int numDisplaceAccepted;    ///< Number of equil Displace moves accepted
        int numMuAttempted;         ///< Number of moves between mu adjustments
        int numNAttempted;          ///< The number of particle measurements
        int numStepsAttempted;      ///< Number of steps for relaxing C0

        void printWormState();
        std::string printHistogram();

    private:
        MTRand &random;             // The global random number generator

        int configNumber;           // The output configuration number
        int numImagTimeSweeps;      // Partitioning used for updates
        int numSteps;               // The number of steps performed during equilibration
        int numUpdates;             // The maximum number of updates per step.
        uint32 binSize;                // The maximum number measurements per bin.
        int numParticles;           // The number of particles

        int N0;                     ///< The initial/target number of particles
        bool foundmu;               ///< Have we found the optimal μ?
        bool foundC0;               ///< Have we found the optimal C0?
        double muFactor;            ///< Used for scaling muShift stepsize
        double bestmu;              ///< The previously best value of mu.
        int  bestPN;                ///< The previously best value of PN.
        double bestDiffAveN;        ///< The currently best value of aveN;
        bool inWindow;              ///< Is the number of particles in the window?

        int sgnAveN;                ///< What is the relative sign betwen aveN and N0?

	blitz::Array<int,1> PN;            ///< Number probability distribution (used in relaxmu)
        bool relaxmuMessage;        ///< For printing a message when μ relaxation begins
        bool relaxC0Message;        ///< For printing a message when μ relaxation begins
        bool equilMessage;          ///< For printing a message when equilibration begins
        bool equilODMessage;        ///< For printing a message when OD equilibration begins
        std::vector <double> C0Vals;     ///< Record of previously attempted values of C0
        std::vector <double> diagFracVals;   ///< Record fo diagonal fraction (for linear regression).

        double targetDiagFrac;      ///< The target diagonal fraction we are trying to obtain.
        int sgnDiagFrac;            ///< What is the relative sign for C0 relaxation
        double shiftC0;             ///< How much are we shifting C0?
        uint32 barStepSize;         ///< Used for printing a status bar.
        uint32 numRemainingSteps;   ///< The number of remaining steps in the equilibration stage.
        uint32 numBars;             ///< Used for printing out progress bar

        uint32 Npaths;                 // Number of paths
    
        std::vector<std::string> stateStrings; // A std::vector of state std::strings from the last bin
		std::vector<StateObject> stateVector; // A std::vector of state objects from the last bin

        bool startWithState;        // Are we starting from a saved state
        bool success;               // Track move success/failure
    
        boost::ptr_vector<Path> &pathPtrVec; // The std::vector of all paths
        Path &path;                          // A reference to the first path in the path std::vector

        boost::ptr_vector<move_vector> &movePtrVec;  // Vector of all move_vectors
        move_vector &move;                           // A reference to the first move_vector
    
        boost::ptr_vector<estimator_vector> &estimatorPtrVec;  // Vector of estimators for 
                                                               // each path followed by multipath
        estimator_vector &estimator;                           // A reference to the first estimator_vector

        std::vector <double> attemptDiagProb;        // The cumulative diagonal attempt Probabilities
        std::vector <double> attemptOffDiagProb;     // The cumulative off-diagonal attempt Probabilities

        std::map <std::string,int> moveIndex;             // A std::map to keep track of move names and indices
        std::map <std::string,int> estimatorIndex;        // A std::map to keep track of estimator names and indices

        /* Output estimators to disk */
        void output();

        /* perform an iterative linear regrssion for C0 calculation */
        double linearRegressionC0();

        /* Load the PIMC state from disk */
        void loadState();

        /* Load classical or quantum initial states */
        void loadClassicalState(blitz::Array <dVec,2>&, blitz::Array <unsigned int, 2>&, int);
        void loadQuantumState(blitz::Array <dVec,2>&, blitz::Array <beadLocator,2>&, blitz::Array<beadLocator,2>&, int, int);

        /* Shuffle the offDiag move list, returning it in moves*/
        void shuffleMoves();

        /* Perform a Metropolis update */
        std::string update(const double,const int,const int);

};

#endif

