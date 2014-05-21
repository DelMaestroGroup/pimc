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

        virtual EstimatorBase* clone() const = 0;
    
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
    
        /* Append to default label */
        void appendLabel(string append);

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

inline EstimatorBase* new_clone(EstimatorBase const& other){
    return other.clone();
}

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
    
    EnergyEstimator* clone() const{ return new EnergyEstimator(*this); }

	private:
		ActionBase *actionPtr;
		void accumulate();		// Accumulate values

};


// ========================================================================  
// Thermodynamic and Virial Energy Estimator Class 
// ========================================================================  
/** 
 * Computes the total energy via the thermodynamic estimator.
 * Also computes total energy via the centroid virial estimator.
 *
 * Measures the total potential and kinetic energy, as well as the
 * per-particle values using the thermodynamic estimator.
 *
 * Measures the virial kinetic energy and subtracts this off of the 
 * total virial energy to get the potential --> useful for non-local
 * actions.
 *
 * Also measures the total virial energy and subtracts off the potential
 * which was computed via the operator method --> useful for local actions.
 *
 *
 * @see S. Jang, S. Jang and G.A. Voth, J. Chem. Phys. 115, 7832 (2001).
 * @see W. Janke and T. Sauer, J. Chem. Phys. 107, 15 (1997).
 */
class VirialEnergyEstimator: public EstimatorBase {

	public:
		VirialEnergyEstimator(const Path &, ActionBase *, 
                int _frequency=1, string _label="estimator");
		~VirialEnergyEstimator();

        VirialEnergyEstimator* clone() const{ return new VirialEnergyEstimator(*this); }

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
    
        NumberParticlesEstimator* clone() const{ return new NumberParticlesEstimator(*this); }

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
    
        ParticlePositionEstimator* clone() const{ return new ParticlePositionEstimator(*this); }
    
        void output();              // overload the output

    private:
        void accumulate();			// Accumulate values
};

// ========================================================================  
// BIPARTITION DENSITY ESTIMATOR CLASS 
// ========================================================================  
/** 
 * Compute density inside film and in bulk separately for Excluded
 * volume potentials.
 */
class BipartitionDensityEstimator: public EstimatorBase {

	public:
		BipartitionDensityEstimator(const Path &, ActionBase *,
                int _frequency=1, string _label="bipart_dens");
		~BipartitionDensityEstimator();

        BipartitionDensityEstimator* clone() const{ return new BipartitionDensityEstimator(*this); }

	private:
		ActionBase *actionPtr;
		void accumulate();		// Accumulate values
};

// ========================================================================  
// Plane Particle Density Estimator Class 
// ========================================================================  
/**
 * Create a 2d histogram of particle positions.
 * 
 */
class PlaneParticlePositionEstimator: public EstimatorBase {

	public:
		PlaneParticlePositionEstimator(const Path &, 
                int _frequency=1, string _label="planedensity");
		~PlaneParticlePositionEstimator();

        PlaneParticlePositionEstimator* clone() const{ return new PlaneParticlePositionEstimator(*this); }
    
    private:
        int numGrid;                // The number of grid points
        dVec dl;                    // The linear size of each spatial bin
        void accumulate();			// Accumulate values
        dVec side;                  // Local copy of container geometry
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
    
        NumberDistributionEstimator* clone() const{ return new NumberDistributionEstimator(*this); }

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

        SuperfluidFractionEstimator* clone() const{ return new SuperfluidFractionEstimator(*this); }
    
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

        LocalSuperfluidDensityEstimator* clone() const{ return new LocalSuperfluidDensityEstimator(*this); }
    
        void output();           ///< overload the output

	private:
        int numGrid;            // The number of grid points
        double dR;              // The size of the radial bin
        Array <double,1> locAz; // The local path area estimator
        Array <double,1> locA2; // The local area squared
        Array <double,1> locWz; // The local winding number estimator

		void accumulate();		// Accumulate values

};

// ========================================================================  
// Plane Winding Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class PlaneWindingSuperfluidDensityEstimator: public EstimatorBase {

	public:
		PlaneWindingSuperfluidDensityEstimator(const Path &, 
                int _frequency=1, string _label="planewind");
		~PlaneWindingSuperfluidDensityEstimator();
    
        PlaneWindingSuperfluidDensityEstimator* clone() const{ return new PlaneWindingSuperfluidDensityEstimator(*this); }

	private:
        dVec side;              // A local copy of the dimensions
        double dx;              // The linear x-size of the spatial bin
        double dy;              // The linear y-size of the spatial bin
        int numGrid;            // The number of grid points
        Array <double,1> locWz; // The local winding number estimator
		void accumulate();		// Accumulate values

};

// ========================================================================  
// Plane Area Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class PlaneAreaSuperfluidDensityEstimator: public EstimatorBase {

	public:
		PlaneAreaSuperfluidDensityEstimator(const Path &, 
                int _frequency=1, string _label="planearea");
		~PlaneAreaSuperfluidDensityEstimator();

        PlaneAreaSuperfluidDensityEstimator* clone() const{ return new PlaneAreaSuperfluidDensityEstimator(*this); }
    
	private:
        dVec side;              // A local copy of the dimensions
        double dx;              // The linear x-size of the spatial bin
        double dy;              // The linear y-size of the spatial bin
        int numGrid;            // The number of grid points
        Array <double,1> locAz; // The local area estimator
		void accumulate();		// Accumulate values

};

// ========================================================================  
// Radial Winding Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class RadialWindingSuperfluidDensityEstimator: public EstimatorBase {

	public:
		RadialWindingSuperfluidDensityEstimator(const Path &, 
                int _frequency=1, string _label="radwind");
		~RadialWindingSuperfluidDensityEstimator();

        RadialWindingSuperfluidDensityEstimator* clone() const{ return new RadialWindingSuperfluidDensityEstimator(*this); }
    
	private:
        double dR;              // The size of the radial bin
        int numGrid;            // The number of grid points
        Array <double,1> locWz; // The local winding number estimator
		void accumulate();		// Accumulate values

};

// ========================================================================  
// Radial Area Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class RadialAreaSuperfluidDensityEstimator: public EstimatorBase {

	public:
		RadialAreaSuperfluidDensityEstimator(const Path &, 
                int _frequency=1, string _label="radarea");
		~RadialAreaSuperfluidDensityEstimator();
    
        RadialAreaSuperfluidDensityEstimator* clone() const{ return new RadialAreaSuperfluidDensityEstimator(*this); }

	private:
        double dR;              // The size of the radial bin
        int numGrid;            // The number of grid points
        Array <double,1> locAz; // The local winding number estimator
		void accumulate();		// Accumulate values

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
    
        DiagonalFractionEstimator* clone() const{ return new DiagonalFractionEstimator(*this); }

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
    
        PermutationCycleEstimator* clone() const{ return new PermutationCycleEstimator(*this); }

	private:
		Array <bool,1> doBead;		// Used for ensuring we don't double count beads
		int maxNumCycles;			// The maximum number of cycles to consider
		void accumulate();			// Accumulate values
};

// ========================================================================  
// Local Permutation Cycle Estimator Class 
// ========================================================================  
/** 
 * Particle permutation number density histogram.
 * NORMALIZATION NEEDS TO BE FIXED!! 
 */
class LocalPermutationEstimator: public EstimatorBase {

	public:
		LocalPermutationEstimator(const Path &,
                int _frequency=1, string _label="locperm");
		~LocalPermutationEstimator();

        LocalPermutationEstimator* clone() const{ return new LocalPermutationEstimator(*this); }

        void output();

	private:
        Array <int, 1> numBeadInGrid;
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
		OneBodyDensityMatrixEstimator(Path &, ActionBase *, const MTRand &,
                int _frequency=20, string _label="obdm");
		~OneBodyDensityMatrixEstimator();
    
        OneBodyDensityMatrixEstimator* clone() const{ return new OneBodyDensityMatrixEstimator(*this); }

		void sample();				// Sample the estimator
		void outputFooter();		// Output the acceptance footer to disk

	private:
		Path &lpath;					// A non-constant local reference to the path
		ActionBase *actionPtr;			// The action pointer
		MTRand random;					// A local copy of the random number generator

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
    
        PairCorrelationEstimator* clone() const{ return new PairCorrelationEstimator(*this); }

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
    
        RadialDensityEstimator* clone() const{ return new RadialDensityEstimator(*this); }

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
    
        WormPropertiesEstimator* clone() const{ return new WormPropertiesEstimator(*this); }

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
    
        CylinderEnergyEstimator* clone() const{ return new CylinderEnergyEstimator(*this); }

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
    
        CylinderNumberParticlesEstimator* clone() const{ return new CylinderNumberParticlesEstimator(*this); }

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
    
        CylinderNumberDistributionEstimator* clone() const{ return new CylinderNumberDistributionEstimator(*this); }

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
    
        CylinderSuperfluidFractionEstimator* clone() const{ return new CylinderSuperfluidFractionEstimator(*this); }

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
		CylinderOneBodyDensityMatrixEstimator(Path &, ActionBase *, const MTRand &, double, 
                int _frequency=20, string _label="cyl_obdm");
		~CylinderOneBodyDensityMatrixEstimator();

        CylinderOneBodyDensityMatrixEstimator* clone() const{ return new CylinderOneBodyDensityMatrixEstimator(*this); }
    
		void sample();				// Sample the estimator

	private:
		Path &lpath;					// A non-constant local reference to the path
		ActionBase *actionPtr;			// The action pointer
		MTRand random;					// A local copy of the random number generator

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
    
        CylinderPairCorrelationEstimator* clone() const{ return new CylinderPairCorrelationEstimator(*this); }

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
    
        CylinderRadialPotentialEstimator* clone() const{ return new CylinderRadialPotentialEstimator(*this); }

	private:
		ActionBase *actionPtr;			// A pointer that will be used to get the potential
		MTRand &random;					// A local reference to the random number generator
		
		double maxR;					// The estimator cutoff radius
		double dR;						// The discretization
		Array <double,1> radPot;		// Used for normalization

		void accumulate();				// Accumulate values
		void accumulate1();				// Accumulate values
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// BEGIN PIGS ESTIMATORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// ========================================================================  
// Potential Energy Estimator Class 
// ========================================================================  
/** 
 * Computes the potential energy along the worldline.
 *
 */
class PotentialEnergyEstimator: public EstimatorBase {

	public:
		PotentialEnergyEstimator(const Path &, ActionBase *, 
                int _frequency=1, string _label="potential");
		~PotentialEnergyEstimator();
        PotentialEnergyEstimator* clone() const{ return new PotentialEnergyEstimator(*this); }
    
	private:
		ActionBase *actionPtr;
		void accumulate();		// Accumulate values

};

// ========================================================================
// Kinetic Energy Estimator Class
// ========================================================================
/**
 * Computes the total energy using a mixed estimator
 */
class KineticEnergyEstimator: public EstimatorBase {
    
public:
    KineticEnergyEstimator(const Path &, ActionBase *,
                           int _frequency=1, string _label="kinetic");
    ~KineticEnergyEstimator();
    
    KineticEnergyEstimator* clone() const{ return new KineticEnergyEstimator(*this); }
    
private:
    ActionBase *actionPtr;
    void accumulate();		// Accumulate values
    
};

// ========================================================================
// PIGS Energy Estimator Class
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
class PigsEnergyEstimator: public EstimatorBase {
    
public:
    PigsEnergyEstimator(const Path &, ActionBase *,
                    int _frequency=1, string _label="estimator");
    ~PigsEnergyEstimator();
    
    PigsEnergyEstimator* clone() const{ return new PigsEnergyEstimator(*this); }
    
private:
    ActionBase *actionPtr;
    void accumulate();		// Accumulate values
    
};


// ========================================================================
// Total Energy Estimator Class
// ========================================================================
/**
 * Computes the total energy using a mixed estimator
 */
class TotalEnergyEstimator: public EstimatorBase {
    
public:
    TotalEnergyEstimator(const Path &, ActionBase *,
                           int _frequency=1, string _label="energy");
    ~TotalEnergyEstimator();

    TotalEnergyEstimator* clone() const{ return new TotalEnergyEstimator(*this); }
    
private:
    ActionBase *actionPtr;
    void accumulate();		// Accumulate values
    
};

// ========================================================================
// Themodynamic Potential Energy Estimator Class
// ========================================================================
/**
 * Computes the total energy using a mixed estimator
 */
class ThermoPotentialEnergyEstimator: public EstimatorBase {
    
public:
    ThermoPotentialEnergyEstimator(const Path &, ActionBase *,
                         int _frequency=1, string _label="thpotenial");
    ~ThermoPotentialEnergyEstimator();
    
    ThermoPotentialEnergyEstimator* clone() const{ return new ThermoPotentialEnergyEstimator(*this); }
    
private:
    ActionBase *actionPtr;
    void accumulate();		// Accumulate values
    
};

//// ========================================================================
//// Local Energy Estimator Class
//// ========================================================================
///**
// * Computes the local energy along the worldline.
// *
// */
//class LocalEnergyEstimator: public EstimatorBase {
//    
//public:
//    LocalEnergyEstimator(const Path &, ActionBase *,
//                             int _frequency=1, string _label="local");
//    ~LocalEnergyEstimator();
//    
//private:
//    ActionBase *actionPtr;
//    void accumulate();		// Accumulate values
//    
//};
//
//// ========================================================================  
//// Total Energy Estimator Class 
//// ========================================================================  
///** 
// * Computes the total energy using a mixed estimator
// */
//class TotalEnergyEstimator: public EstimatorBase {
//
//	public:
//		TotalEnergyEstimator(const Path &, ActionBase *, const MTRand &,
//                int _frequency=1, string _label="energy");
//		~TotalEnergyEstimator();
//
//	private:
//		ActionBase *actionPtr;
// 		MTRand random;					// A local copy of the random number generator
//        dVec newPosition(const beadLocator &);
//		void accumulate();		// Accumulate values
//
//};

// ========================================================================  
// Position Estimator Class 
// ========================================================================  
/** 
 * Computes the average value of the position in 1D.
 *
 */
class PositionEstimator: public EstimatorBase {

	public:
		PositionEstimator(const Path &, int _frequency=1, 
                string _label="position");
		~PositionEstimator();

    PositionEstimator* clone() const{ return new PositionEstimator(*this); }

	private:
		void accumulate();		// Accumulate values

};

// ========================================================================
// Particle Resolved Position Estimator Class
// ========================================================================
/**
 * Computes the average position of each particle in 1D at the center time
 * slice
 */
class ParticleResolvedPositionEstimator: public EstimatorBase {
    
public:
    ParticleResolvedPositionEstimator(const Path &, int _frequency=1,
                      string _label="prposition");
    ~ParticleResolvedPositionEstimator();

    ParticleResolvedPositionEstimator* clone() const{ return new ParticleResolvedPositionEstimator(*this); }
    
private:
    void accumulate();		// Accumulate values
    
};


// ========================================================================
// Particle Correlation Estimator Class
// ========================================================================
/**
 * Computes the average position of each particle in 1D at the center time
 * slice
 */
class ParticleCorrelationEstimator: public EstimatorBase {
    
public:
    ParticleCorrelationEstimator(const Path &, int _frequency=1,
                                      string _label="prcorrelation");
    ~ParticleCorrelationEstimator();
    
    ParticleCorrelationEstimator* clone() const{ return new ParticleCorrelationEstimator(*this); }

private:
    void accumulate();		// Accumulate values
    
};


// ========================================================================
// Velocity Estimator Class
// ========================================================================
/**
 * Computes the imaginary time resolved "velocity" for the first particle .
 *
 */
class VelocityEstimator: public EstimatorBase {
    
public:
    VelocityEstimator(const Path &, int _frequency=1,
                      string _label="velocity");
    ~VelocityEstimator();

    VelocityEstimator* clone() const{ return new VelocityEstimator(*this); }
    
private:
    void accumulate();		// Accumulate values
    
};

// ========================================================================
// SubregionOccupation Estimator Class
// ========================================================================
/**
 * Computes the imaginary time resolved "velocity" for the first particle .
 *
 */
class SubregionOccupationEstimator: public EstimatorBase {
    
public:
    SubregionOccupationEstimator(const Path &,ActionBase *_actionPtr,
                                 int _frequency=1,
                      string _label="subregionocc");
    ~SubregionOccupationEstimator();
    
    SubregionOccupationEstimator* clone() const{ return new SubregionOccupationEstimator(*this); }
    
private:
    ActionBase *actionPtr;
    void accumulate();		// Accumulate values
    
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
class PIGSOneBodyDensityMatrixEstimator: public EstimatorBase {
    
public:
    PIGSOneBodyDensityMatrixEstimator(Path &, ActionBase *, const MTRand &,
                                  int _frequency=1, string _label="obdm");
    ~PIGSOneBodyDensityMatrixEstimator();
    
    PIGSOneBodyDensityMatrixEstimator* clone() const{ return new PIGSOneBodyDensityMatrixEstimator(*this); }
    
    void sample();				// Sample the estimator
    void outputFooter();		// Output the acceptance footer to disk
    
private:
    Path &lpath;					// A non-constant local reference to the path
    ActionBase *actionPtr;			// The action pointer
    MTRand random;					// A local copy of the random number generator
    
    double dR;						// The discretization
    int numReps;					// The number of measurments reps
    uint32 numAccepted;				// The number of moves accepted
    uint32 numAttempted;			// The number of moves attempted
    
    
    dVec newTailPos,oldTailPos;		// The new and old tail position
    dVec newRanPos,neighborPos;		// The random shift
    
    double sqrt2LambdaTau;			// sqrt(2 * lambda * tau)
    double rho0Norm;				// Free density matrix
    double oldAction,newAction;		// The old and new action
    
    /* Get a random vector */
    dVec getRandomVector(const double);
    
    /* Accumulate values */
    void accumulate();	
};


// ========================================================================
// Doubled Estimator Base Class
// ========================================================================
/**
 * Base class for estimators that use two paths
 */
class DoubledEstimator: public EstimatorBase {
    
public:
    DoubledEstimator(const Path &, const Path &,
                    int _frequency=1, string _label="");
    ~DoubledEstimator();
   
protected:
    const Path &path2;				///< A constant reference to the paths
    
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// END DOUBLED ESTIMATOR BASE CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


// ========================================================================
// Swap Estimator Class
// ========================================================================
/**
 * Computes the Swap Estimator between two paths.
 */
class SwapEstimator: public DoubledEstimator {
    
public:
    SwapEstimator(Path &, Path &,ActionBase *,ActionBase *,
                             int _frequency=1, string _label="swap");
    ~SwapEstimator();
    
    SwapEstimator* clone() const{ return new SwapEstimator(*this); }
    
private:
    Path &lpath;					// A non-constant local reference to path 1
    Path &lpath2;					// A non-constant local reference to path 2
    ActionBase *actionPtr;			// The action pointer for path 1
    ActionBase *actionPtr2;			// The action pointer for path 2
    
    void accumulate();			// Accumulate values
    void accumulateOpen();			// Accumulate values for open paths
    void accumulateClosed();			// Accumulate values for open paths
};


// ========================================================================
// Swap Estimator Class
// ========================================================================
/**
 * Computes the Swap Estimator between two paths.
 */
class EntPartEstimator: public DoubledEstimator {
    
public:
    EntPartEstimator(Path &, Path &,ActionBase *,ActionBase *,
                  int _frequency=1, string _label="entpart");
    ~EntPartEstimator();
    
    EntPartEstimator* clone() const{ return new EntPartEstimator(*this); }
    
private:
    Path &lpath;					// A non-constant local reference to path 1
    Path &lpath2;					// A non-constant local reference to path 2
    ActionBase *actionPtr;			// The action pointer for path 1
    ActionBase *actionPtr2;			// The action pointer for path 2
    
    void accumulate();			// Accumulate values
};



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// END PIGS ESTIMATORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
#endif
