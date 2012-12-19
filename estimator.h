/**
 * @file estimator.h
 * @author Adrian Del Maestro
 * @date 10.23.2008
 *
 * @brief Estimator class definitions.
 */

#include "common.h"

#ifndef ESTIMATOR_H 
#define ESTIMATOR_H

class Path;
class ActionBase;
class Potential;

// ========================================================================  
// EstimatorBase Class
// ========================================================================  
/**
 * The base class that all estimator classes will be derived from. 
 *
 * Contains the common elements of all estimators, including the ability
 * to initialize, accumulate, output and reset.  We assume a default measurment
 * frequency of every time.
 */
class EstimatorBase {

	public:
		EstimatorBase (const Path &, int _frequency=1, string _label="");
		virtual ~EstimatorBase() = 0;

		/* Sample the estimator */
		virtual void sample();		
		
		/* Reset the estimator */
		void reset();								

		/* Restart the estimator */
		void restart(const uint32, const uint32); 

		/* Output the estimator */
		virtual void output();					

		/** Ouptut the fooder to disk */
		virtual void outputFooter() {};

        /* Evaluate the basic conditions for sampling */
        bool baseSample();

		/** Get the total number of accumulated measurments */
		uint32 getTotNumAccumulated() { return totNumAccumulated; }

		/** Get the number of accumulated measurements since the last reset */
		uint32 getNumAccumulated() { return numAccumulated; }

		/** Get the number of samples since the last reset */
		uint32 getNumSampled() { return numSampled; }

		/** Get the name of the estimator */
		string getName() { return name; }

		/* Prepare the estimator for i/o */
		void prepare();

	protected:
		const Path &path;				///< A constant reference to the paths

		fstream *outFilePtr;			///< The output fie

		Array<double,1> estimator;		///< The estimator array
		Array<double,1> norm;			///< The normalization factor for each estimator

		int numEst;						///< The number of individual quantities measured
		string name;	                ///< The name of the estimator
        string label;                   ///< The label used for the output file
		uint32 numSampled;				///< The number of times we have sampled
		uint32 numAccumulated;			///< The number of accumulated values
		uint32 totNumAccumulated;		///< The total number of accumulated values
		int frequency;					///< The frequency at which we accumulate
		int numBeads0;					///< The target number of beads for the canonical ensemble

		bool diagonal;					///< Is this a diagonal estimator?
		bool endLine;					///< Should we output a carriage return?
		bool canonical;					///< Are we in the canonical ensemble?

		string header;					///< The data file header

		/** Accumulate the estimator */
		virtual void accumulate() = 0;	

		/* Initialize the estimator */
		void initialize(int);
};


// ========================================================================  
// Energy Estimator Class 
// ========================================================================  
/** 
 * Computes the total energy via the thermodynamic estimator.
 *
 * Measures the total potential and kinetic energy, as well as the
 * per-particle values using the thermodynamic estimator.
 *
 * @see S. Jang, S. Jang and G.A. Voth, J. Chem. Phys. 115, 7832 (2001).
 * @see W. Janke and T. Sauer, J. Chem. Phys. 107, 15 (1997).
 */
class EnergyEstimator: public EstimatorBase {

	public:
		EnergyEstimator(const Path &, ActionBase *, 
                int _frequency=1, string _label="estimator");
		~EnergyEstimator();

	private:
		ActionBase *actionPtr;
		void accumulate();		// Accumulate values

};

// ========================================================================  
// Number Particles Estimator Class 
// ========================================================================  
/**
 * Computes the average number of particles, as well as density.
 */
class NumberParticlesEstimator: public EstimatorBase {

	public:
		NumberParticlesEstimator(const Path &,
                int _frequency=1, string _label="estimator");
		~NumberParticlesEstimator();

	private:
		void accumulate();			// Accumulate values
};

// ========================================================================  
// Particle Density Estimator Class 
// ========================================================================  
/**
 * Create histogram of particle positions.
 * 
 */
class ParticlePositionEstimator: public EstimatorBase {

	public:
		ParticlePositionEstimator(const Path &, 
                int _frequency=1, string _label="position");
		~ParticlePositionEstimator();
        
        void output();              // overload the output

    private:
        void accumulate();			// Accumulate values
};


// ========================================================================  
// Number Distribution Estimator Class 
// ========================================================================  
/**
 * Computes the probability distribution function for the number of particles.
 */
class NumberDistributionEstimator: public EstimatorBase {

	public:
		NumberDistributionEstimator(const Path &,
                int _frequency=1, string _label="number");
		~NumberDistributionEstimator();

	private:
		int startParticleNumber;	// The smallest number of particles
		int endParticleNumber;		// The largest  number of particles
		int maxNumParticles;		// The particle Range
		int particleShift;			// The shift in the particle array

		void accumulate();			// Accumulate values
};

// ========================================================================  
// Superfluid Fraction Estimator Class 
// ========================================================================  
/**
 * Compute the superfluid fraction, as well as the winding number 
 * probability distribution.
 */
class SuperfluidFractionEstimator: public EstimatorBase {

	public:
		SuperfluidFractionEstimator(const Path &, 
                int _frequency=1, string _label="super");
		~SuperfluidFractionEstimator();

	private:
		int windMax;			// The maximum winding number considered

		void accumulate();		// Accumulate values

};

// ========================================================================  
// Local Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the local superfluid density.
 */
class LocalSuperfluidDensityEstimator: public EstimatorBase {

	public:
		LocalSuperfluidDensityEstimator(const Path &, 
                int _frequency=1, string _label="locsuper");
		~LocalSuperfluidDensityEstimator();

        void output();           ///< overload the output

	private:
		void accumulate();		// Accumulate values
        Array <double,1> locAz; // The local path area
        Array <dVec,1> locW;    // The local winding number

};

// ========================================================================  
// Diagonal Fraction Estimator Class 
// ========================================================================  
/** 
 * Compute the fraction of time we spend in the diagonal ensemble.
 */
class DiagonalFractionEstimator: public EstimatorBase {

	public:
		DiagonalFractionEstimator(const Path &,
                int _frequency=1, string _label="estimator");
		~DiagonalFractionEstimator();

        void sample();          // Overload to always-on sampling

	private:
		void accumulate();		// Accumulate values
};

// ========================================================================  
// Permutation Cycle Estimator Class 
// ========================================================================  
/** 
 * Computes the particle permutation cycle probability distribution.
 */
class PermutationCycleEstimator: public EstimatorBase {

	public:
		PermutationCycleEstimator(const Path &,
                int _frequency=1, string _label="pcycle");
		~PermutationCycleEstimator();

	private:
		Array <bool,1> doBead;		// Used for ensuring we don't double count beads
		int maxNumCycles;			// The maximum number of cycles to consider
		void accumulate();			// Accumulate values
};

// ========================================================================  
// One Body Density Matrix Estimator Class 
// ========================================================================  
/** 
 * Compute the one body density matrix n(r) which can be used to find the
 * momentum distribution function and structure factor. 
 *
 * We use a reference to the global random number generator, and thus 
 * this estimator can effect final results.  We also need a non-constant
 * local reference to the Path which is used for temporary updates.
 */
class OneBodyDensityMatrixEstimator: public EstimatorBase {

	public:
		OneBodyDensityMatrixEstimator(Path &, ActionBase *, MTRand &,
                int _frequency=20, string _label="obdm");
		~OneBodyDensityMatrixEstimator();

		void sample();				// Sample the estimator
		void outputFooter();		// Output the acceptance footer to disk

	private:
		Path &lpath;					// A non-constant local reference to the path
		ActionBase *actionPtr;			// The action pointer
		MTRand &random;					// A local reference to the random number generator

		double dR;						// The discretization
		int numReps;					// The number of measurments reps					
		uint32 numAccepted;				// The number of moves accepted
		uint32 numAttempted;			// The number of moves attempted


	 	dVec newTailPos,oldTailPos;		// The new and old tail position
	 	dVec newHeadPos;				// The new head position
		dVec newRanPos,neighborPos;		// The random shift

		double sqrt2LambdaTau;			// sqrt(2 * lambda * tau)
		double rho0Norm;				// Free density matrix
		double oldAction,newAction;		// The old and new action

		/* Get a random vector */
		dVec getRandomVector(const double);
		
		/* Returns a position used by the staging algorithm */
		dVec newStagingPosition(const beadLocator &, const beadLocator &, const int, const int);

		/* Accumulate values */
		void accumulate();	

};

// ========================================================================  
// Pair Correlation Estimator Class
// ========================================================================  
/** 
 * Compute the two-body pair-correlation function, g(r) ~ <rho(r)rho(0)>.
 */
class PairCorrelationEstimator: public EstimatorBase {

	public:
		PairCorrelationEstimator(const Path &, ActionBase *,
                int _frequency=1, string _label="pair");
		~PairCorrelationEstimator();

	private:
		ActionBase *actionPtr;			// A pointer that will be used to get the potential
		void accumulate();				// Accumulate values
		double dR;						// The discretization
};

// ========================================================================  
// Radial Density Estimator Class
// ========================================================================  
/** 
 * Compute the density as a function of position in the radial direction.
 */
class RadialDensityEstimator: public EstimatorBase {

	public:
		RadialDensityEstimator(const Path &, 
                int _frequency=1, string _label="radial");
		~RadialDensityEstimator();

	private:
		void accumulate();				// Accumulate values
		double dR;						// The discretization
};

// ========================================================================  
// Worm Properties Estimator Class
// ========================================================================  
/** 
 * Compute various properties related to the worm in the simulation.
 */
class WormPropertiesEstimator: public EstimatorBase {

	public:
		WormPropertiesEstimator(const Path &, 
                int _frequency=100, string _label="worm");
		~WormPropertiesEstimator();

	private:
		dVec sep;						// head-tail separation
		void accumulate();				// Accumulate values
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// CYLINDER ESTIMATORS ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// ========================================================================  
// Cylinder Energy Estimator Class 
// ========================================================================  
/** 
 * Computes the total energy via the thermodynamic estimator.
 *
 * Measures the total potential and kinetic energy, as well as the
 * per-particle values using the thermodynamic estimator.
 *
 * @see S. Jang, S. Jang and G.A. Voth, J. Chem. Phys. 115, 7832 (2001).
 * @see W. Janke and T. Sauer, J. Chem. Phys. 107, 15 (1997).
 */
class CylinderEnergyEstimator: public EstimatorBase {

	public:
		CylinderEnergyEstimator(const Path &, ActionBase *, double, 
                int _frequency=1, string _label="cyl_estimator");
		~CylinderEnergyEstimator();

	private:
		ActionBase *actionPtr;
		double maxR;			// The estimator cutoff radius

		void accumulate();		// Accumulate values

};

// ========================================================================  
// Cylinder Number Particles Estimator Class 
// ========================================================================  
/**
 * Computes the average number of particles, as well as density.
 */
class CylinderNumberParticlesEstimator: public EstimatorBase {

	public:
		CylinderNumberParticlesEstimator(const Path &, double, 
                int _frequency=1, string _label="cyl_estimator");
		~CylinderNumberParticlesEstimator();

	private:
		double maxR;			// The estimator cutoff radius

		void accumulate();		// Accumulate values
};

// ========================================================================  
// Cylinder Number Distribution Estimator Class 
// ========================================================================  
/**
 * Computes the probability distribution function for the number of particles.
 */
class CylinderNumberDistributionEstimator: public EstimatorBase {

	public:
		CylinderNumberDistributionEstimator(const Path &, double, 
                int _frequency=1, string _label="cyl_number");
		~CylinderNumberDistributionEstimator();

	private:
		double maxR;				// The estimator cutoff radius
		int maxNumParticles;		// The maximum number considered

		void accumulate();			// Accumulate values
};

// ========================================================================  
// Cylinder Superfluid Fraction Estimator Class 
// ========================================================================  
/**
 * Compute the superfluid fraction, as well as the winding number 
 * probability distribution.
 */
class CylinderSuperfluidFractionEstimator: public EstimatorBase {

	public:
		CylinderSuperfluidFractionEstimator(const Path &, double, 
                int _frequency=1, string _label="cyl_super");
		~CylinderSuperfluidFractionEstimator();

	private:
		Array <bool,1> doBead;	// Used for ensuring we don't double count beads
		int windMax;			// The maximum winding number considered
		double maxR;			// The estimator cutoff radius

		void accumulate();		// Accumulate values

};

// ========================================================================  
// Cylinder One Body Density Matrix Estimator Class 
// ========================================================================  
/** 
 * Compute the one body density matrix n(r) which can be used to find the
 * momentum distribution function and structure factor. 
 *
 * We use a reference to the global random number generator, and thus 
 * this estimator can effect final results.  We also need a non-constant
 * local reference to the Path which is used for temporary updates.
 */
class CylinderOneBodyDensityMatrixEstimator: public EstimatorBase {

	public:
		CylinderOneBodyDensityMatrixEstimator(Path &, ActionBase *, MTRand &, double, 
                int _frequency=20, string _label="cyl_obdm");
		~CylinderOneBodyDensityMatrixEstimator();

		void sample();				// Sample the estimator

	private:
		Path &lpath;					// A non-constant local reference to the path
		ActionBase *actionPtr;			// The action pointer
		MTRand &random;					// A local reference to the random number generator

		double dR;						// The discretization
		int numReps;					// The number of measurments reps					
		uint32 numAccepted;				// The number of moves accepted
		uint32 numAttempted;			// The number of moves attempted

		double maxR;					// The estimator cutoff radius


	 	dVec newTailPos,oldTailPos;		// The new and old tail position
	 	dVec newHeadPos;				// The new head position
		dVec newRanPos,neighborPos;		// The random shift

		double sqrt2LambdaTau;			// sqrt(2 * lambda * tau)
		double rho0Norm;				// Free density matrix
		double oldAction,newAction;		// The old and new action

		/* Get a random vector */
		dVec getRandomVector(const double);
		
		/* Returns a position used by the staging algorithm */
		dVec newStagingPosition(const beadLocator &, const beadLocator &, const int, const int);

		/* Accumulate values */
		void accumulate();	
};

// ========================================================================  
// Pair Correlation Estimator Class
// ========================================================================  
/** 
 * Compute the two-body pair-correlation function, g(r) ~ <rho(r)rho(0)>.
 */
class CylinderPairCorrelationEstimator: public EstimatorBase {

	public:
		CylinderPairCorrelationEstimator(const Path &, ActionBase *, double, 
                int _frequency=1, string _label="cyl_pair");
		~CylinderPairCorrelationEstimator();

		void sample();				// Sample the estimator

	private:
		ActionBase *actionPtr;			// A pointer that will be used to get the potential
		double maxR;					// The estimator cutoff radius

		void accumulate();				// Accumulate values
		double dR;						// The discretization
};

// ========================================================================  
// Cylinder Radial Potential Estimator Class
// ========================================================================  
/** 
 * Compute the effective radial potential in a cylinder.
 */
class CylinderRadialPotentialEstimator: public EstimatorBase {

	public:
		CylinderRadialPotentialEstimator(const Path &, ActionBase *, MTRand &, double, 
                int _frequency=1, string _label="cyl_potential");
		~CylinderRadialPotentialEstimator();

	private:
		ActionBase *actionPtr;			// A pointer that will be used to get the potential
		MTRand &random;					// A local reference to the random number generator
		
		double maxR;					// The estimator cutoff radius
		double dR;						// The discretization
		Array <double,1> radPot;		// Used for normalization

		void accumulate();				// Accumulate values
		void accumulate1();				// Accumulate values
};
#endif
