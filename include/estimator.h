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

#ifdef USE_GPU
    #include "common_gpu.h"
#endif

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
        EstimatorBase(const Path & _path, ActionBase *_actionPtr, const MTRand &_random, 
                double _maxR, int _frequency=1, std::string _label="");
        virtual ~EstimatorBase();

        /* Sample the estimator */
        virtual void sample();      
        
        /* Reset the estimator */
        void reset();                               

        /* Restart the estimator */
        void restart(const uint32, const uint32); 

        /* Output the estimator */
        virtual void output();                  

        /* Output a flat average estimator */
        virtual void outputFlat();                  

        /* Output a histogram average estimator */
        virtual void outputHist();                  

        /** Ouptut the fooder to disk */
        virtual void outputFooter() {};

        /* Evaluate the basic conditions for sampling */
        bool baseSample();

        /** Get the total number of accumulated measurments */
        uint32 getTotNumAccumulated() const { return totNumAccumulated; }

        /** Get the number of accumulated measurements since the last reset */
        uint32 getNumAccumulated() const { return numAccumulated; }

        /** Get the number of samples since the last reset */
        uint32 getNumSampled() const { return numSampled; }

        /** Get the name of the estimator */
        virtual std::string getName() const { return "base"; }

        /** Prepare the estimator for i/o */
        void prepare();

        /** Add a carriage return to estimator files */
        void addEndLine() {endLine = true;}; 
    
        /** Append to default label */
        void appendLabel(std::string append);

        /** convert a dVec to std::string */
        std::string dVecToString(const dVec&);

        /** Get the estimator label */
        std::string getLabel() const {return label;};

    protected:
        const Path &path;               ///< A constant reference to the paths
        ActionBase *actionPtr;          ///< A pointer to the action
        MTRand random;                  // A local copy of the random number generator
        double maxR;                    // An estimator cutoff radius

        std::fstream *outFilePtr;            ///< The output fie

        std::map<std::string,int> estIndex;       ///< Map estimator labels to indices.

    DynamicArray<double,1> estimator;      ///< The estimator array
    DynamicArray<double,1> norm;           ///< The normalization factor for each estimator

        int numEst;                     ///< The number of individual quantities measured
        int frequency;                  ///< The frequency at which we accumulate
        int startSlice;                 ///< Where imaginary time averages begin
        int endSlice;                   ///< Where imaginary time averages end
        int endDiagSlice;               ///< Where imaginary time averages end for diagonal estimiators
        std::vector<double> sliceFactor;        ///< Used to properly incorporate end affects 
        std::string label;                   ///< The label used for the output file

        uint32 numSampled;              ///< The number of times we have sampled
        uint32 numAccumulated;          ///< The number of accumulated values
        uint32 totNumAccumulated;       ///< The total number of accumulated values
        int numBeads0;                  ///< The target number of beads for the canonical ensemble

        bool diagonal;                  ///< Is this a diagonal estimator?
        bool endLine;                   ///< Should we output a carriage return?
        bool canonical;                 ///< Are we in the canonical ensemble?

        std::string header;                  ///< The data file header

        /** Accumulate the estimator */
        virtual void accumulate() {}

        /* Initialize the estimator */
        void initialize(int);
        void initialize(std::vector<std::string>);

        /* generate q-vectors needed for momentum space estimators */
        void getQVectors(std::vector<dVec>&);
        void getQVectorsNN(std::vector<dVec>&);
        std::vector <std::vector<dVec> > getQVectors2(double, double, int&, std::string);
};

// ========================================================================  
// Time Estimator Class
// ========================================================================  
/** 
 * An estimator which tracks the ammount of time between bins, summing
 * them into a total at the end. 
 *
 */
class TimeEstimator : public EstimatorBase {

    public:
        TimeEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="estimator");
        ~TimeEstimator() {};
        
        void sample();          // Overload to always-on sampling

        static const std::string name;
        std::string getName() const {return name;}

        void output();              // overload the output

    private:
        void accumulate();      // Accumulate values
                std::chrono::high_resolution_clock::time_point time_begin;
                std::chrono::high_resolution_clock::time_point time_end;

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
        EnergyEstimator (const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="estimator");
        ~EnergyEstimator();

        static const std::string name;
        std::string getName() const {return name;}
    
    private:
        void accumulate();      // Accumulate values
                                
        uint32 numPPAccumulated; ///< The number of per particle (PP) accumulated values

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
        VirialEnergyEstimator (const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="virial");
        ~VirialEnergyEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
        uint32 numPPAccumulated; ///< The number of per particle (PP) accumulated values

};

// ========================================================================  
// Number Particles Estimator Class 
// ========================================================================  
/**
 * Computes the average number of particles, as well as density.
 */
class NumberParticlesEstimator: public EstimatorBase {

    public:
        NumberParticlesEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="estimator");
        ~NumberParticlesEstimator();

        static const std::string name;
        std::string getName() const {return name;}
    
    private:
        void accumulate();          // Accumulate values
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
        ParticlePositionEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="position");
        ~ParticlePositionEstimator();

        static const std::string name;
        std::string getName() const {return name;}
    
        void output() {outputFlat();}  // overload the output

    private:
        void accumulate();          // Accumulate values
        std::vector<std::string> diffLabels;  ///< The axis differential labels
};


#if NDIM > 2
// ========================================================================  
// Commensurate Order Parameter Class
// ========================================================================  
/**
 * A Commensurate/Incommensurate Order Parameter.  Eq. (12) of
 * @see https://journals.aps.org/prb/abstract/10.1103/PhysRevB.73.085422
 */
class CommensurateOrderParameterEstimator: public EstimatorBase {

    public:
        CommensurateOrderParameterEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="estimator");
        ~CommensurateOrderParameterEstimator();

        static const std::string name;
        std::string getName() const {return name;}
    
    private:
        void accumulate();          ///< Accumulate values
        std::vector<dVec> g;             ///< the g-vector set
};
#endif

#if NDIM > 2
// ========================================================================  
// BIPARTITION DENSITY ESTIMATOR CLASS 
// ========================================================================  
/** 
 * Compute density inside film and in bulk separately for Excluded
 * volume potentials.
 */
class BipartitionDensityEstimator: public EstimatorBase {

    public:
        BipartitionDensityEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="bipart_dens");
        ~BipartitionDensityEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
};
#endif

// ========================================================================  
// Linear Particle Density Estimator Class 
// ========================================================================  
/**
 * Create a 1d histogram of particle positions in the z-direction.
 * 
 */
class LinearParticlePositionEstimator: public EstimatorBase {

    public:
        LinearParticlePositionEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="lineardensity");
        ~LinearParticlePositionEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        int numGrid;                // The number of grid points
        double dz;                  // The size of each spatial bin in the z-direction
        void accumulate();          // Accumulate values
        dVec side;                  // Local copy of container geometry
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
        PlaneParticlePositionEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="planedensity");
        ~PlaneParticlePositionEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        int numGrid;                // The total number of grid boxes
        int numLinearGrid;          // The linear number of grid boxes
        dVec dl;                    // The linear size of each spatial bin
        void accumulate();          // Accumulate values
        dVec side;                  // Local copy of container geometry
};

// ========================================================================  
// Averaged Plane Particle Density Estimator Class 
// ========================================================================  
/**
 * Create a 2d histogram of particle positions but only store the average.
 * 
 */
class PlaneParticleAveragePositionEstimator: public EstimatorBase {

    public:
        PlaneParticleAveragePositionEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="planeavedensity");
        ~PlaneParticleAveragePositionEstimator();

        static const std::string name;
        std::string getName() const {return name;}

        void output() {outputFlat();}  // overload the output

    private:
        int numGrid;                // The total number of grid boxes
        int numLinearGrid;          // The linear number of grid boxes
        dVec dl;                    // The linear size of each spatial bin
        void accumulate();          // Accumulate values
        dVec side;                  // Local copy of container geometry
};

// ========================================================================  
// Averaged Plane External Potential Estimator Class 
// ========================================================================  
/**
 * Create a 2d histogram of particle positions but only store the average.
 * 
 */
class PlaneAverageExternalPotentialEstimator: public EstimatorBase {

    public:
        PlaneAverageExternalPotentialEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="planeaveVext");
        ~PlaneAverageExternalPotentialEstimator();

        static const std::string name;
        std::string getName() const {return name;}

        void output();  // overload the output


    private:
        int numGrid;                // The total number of grid boxes
        int numLinearGrid;          // The linear number of grid boxes
        dVec dl;                    // The linear size of each spatial bin
        void accumulate();          // Accumulate values
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
        NumberDistributionEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="number");
        ~NumberDistributionEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        int startParticleNumber;    // The smallest number of particles
        int endParticleNumber;      // The largest  number of particles
        int maxNumParticles;        // The particle Range
        int particleShift;          // The shift in the particle array

        void accumulate();          // Accumulate values
};

#if NDIM > 1
// ========================================================================  
// Superfluid Fraction Estimator Class 
// ========================================================================  
/**
 * Compute the superfluid fraction, as well as the winding number 
 * probability distribution.
 */
class SuperfluidFractionEstimator: public EstimatorBase {

    public:
        SuperfluidFractionEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="super");
        ~SuperfluidFractionEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        int windMax;             // The maximum winding number considered
        double W2Norm;           // A local normalizer for the winding superfluid fraction
        double ANorm;            // A local normalizer for the area superfluid fraction
        uint32 numPPAccumulated; // The number of per particle (PP) accumulated values

        void accumulate();       // Accumulate values

};
#endif

#if NDIM > 1
// ========================================================================  
// Local Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the local superfluid density.
 */
class LocalSuperfluidDensityEstimator: public EstimatorBase {

    public:
        LocalSuperfluidDensityEstimator(const Path &, ActionBase *, const MTRand &, double,
                int _frequency=1, std::string _label="locsuper");
        ~LocalSuperfluidDensityEstimator();

        void output();           ///< overload the output

        static const std::string name;
        std::string getName() const {return name;}

    private:
        int numGrid;            // The number of grid points
        double dR;              // The size of the radial bin
    DynamicArray <double,1> locAz; // The local path area estimator
    DynamicArray <double,1> locA2; // The local area squared
    DynamicArray <double,1> locWz; // The local winding number estimator

        void accumulate();      // Accumulate values

};
#endif

#if NDIM > 1
// ========================================================================  
// Plane Winding Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class PlaneWindingSuperfluidDensityEstimator: public EstimatorBase {

    public:
        PlaneWindingSuperfluidDensityEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="planewind");
        ~PlaneWindingSuperfluidDensityEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        dVec side;              // A local copy of the dimensions
        double dx;              // The linear x-size of the spatial bin
        double dy;              // The linear y-size of the spatial bin
        int numGrid;            // The number of grid points
    DynamicArray <double,1> locWz; // The local winding number estimator
        void accumulate();      // Accumulate values

};
#endif

#if NDIM > 1
// ========================================================================  
// Plane Area Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class PlaneAreaSuperfluidDensityEstimator: public EstimatorBase {

    public:
        PlaneAreaSuperfluidDensityEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="planearea");
        ~PlaneAreaSuperfluidDensityEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        dVec side;              // A local copy of the dimensions
        double dx;              // The linear x-size of the spatial bin
        double dy;              // The linear y-size of the spatial bin
        int numGrid;            // The number of grid points
    DynamicArray <double,1> locAz; // The local area estimator
        void accumulate();      // Accumulate values

};
#endif

#if NDIM > 1
// ========================================================================  
// Radial Winding Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class RadialWindingSuperfluidDensityEstimator: public EstimatorBase {

    public:
        RadialWindingSuperfluidDensityEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="radwind");
        ~RadialWindingSuperfluidDensityEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        double dR;              // The size of the radial bin
        int numGrid;            // The number of grid points
    DynamicArray <double,1> locWz; // The local winding number estimator
        void accumulate();      // Accumulate values

};
#endif

#if NDIM > 1
// ========================================================================  
// Radial Area Superfluid Density Estimator Class 
// ========================================================================  
/**
 * Compute the radially averaged local superfluid density.
 */
class RadialAreaSuperfluidDensityEstimator: public EstimatorBase {

    public:
        RadialAreaSuperfluidDensityEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="radarea");
        ~RadialAreaSuperfluidDensityEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        double dR;              // The size of the radial bin
        int numGrid;            // The number of grid points
    DynamicArray <double,1> locAz; // The local winding number estimator
        void accumulate();      // Accumulate values

};
#endif

// ========================================================================  
// Diagonal Fraction Estimator Class 
// ========================================================================  
/** 
 * Compute the fraction of time we spend in the diagonal ensemble.
 */
class DiagonalFractionEstimator: public EstimatorBase {

    public:
        DiagonalFractionEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="estimator");
        ~DiagonalFractionEstimator();
    
        void sample();          // Overload to always-on sampling

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
};

// ========================================================================  
// Permutation Cycle Estimator Class 
// ========================================================================  
/** 
 * Computes the particle permutation cycle probability distribution.
 */
class PermutationCycleEstimator: public EstimatorBase {

    public:
        PermutationCycleEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="pcycle");
        ~PermutationCycleEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
    DynamicArray <bool,1> doBead;      // Used for ensuring we don't double count beads
        int maxNumCycles;           // The maximum number of cycles to consider
        void accumulate();          // Accumulate values
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
        LocalPermutationEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="locperm");
        ~LocalPermutationEstimator();

        void output() {outputFlat();}

        static const std::string name;
        std::string getName() const {return name;}

    private:
    DynamicArray <int, 1> numBeadInGrid;
    DynamicArray <bool,1> doBead;      // Used for ensuring we don't double count beads
        int maxNumCycles;           // The maximum number of cycles to consider
        void accumulate();          // Accumulate values
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
                double, int _frequency=20, std::string _label="obdm");
        ~OneBodyDensityMatrixEstimator();
    
        void sample();              // Sample the estimator
        void outputFooter();        // Output the acceptance footer to disk

        static const std::string name;
        std::string getName() const {return name;}

    private:
        Path &lpath;                    // A non-constant local reference to the path

        double dR;                      // The discretization
        int numReps;                    // The number of measurments reps                   
        uint32 numAccepted;             // The number of moves accepted
        uint32 numAttempted;            // The number of moves attempted


        dVec newTailPos,oldTailPos;     // The new and old tail position
        //dVec newHeadPos;                // The new head position
        dVec newRanPos,neighborPos;     // The random shift

        double sqrt2LambdaTau;          // sqrt(2 * lambda * tau)
        double rho0Norm;                // Free density matrix
        double oldAction,newAction;     // The old and new action

        /* Get a random std::vector */
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
        PairCorrelationEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="pair");
        ~PairCorrelationEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
        double dR;                      // The discretization
};

// ========================================================================  
// Static Structure Factor Estimator Class
// ========================================================================  
/** 
 * Compute the static structure factor S(q)
 */
class StaticStructureFactorEstimator: public EstimatorBase {

    public:
        StaticStructureFactorEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="ssf");
        ~StaticStructureFactorEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
    DynamicArray <double,1> sf;            // structure factor
        std::vector <std::vector<dVec> > q;       // the q-vectors
};

#ifdef USE_GPU
// ========================================================================  
// GPU Accellerated Static Structure Factor
// ========================================================================  
/** 
 * Compute the intermediate scattering function F(q,\tau)
 */
class StaticStructureFactorGPUEstimator: public EstimatorBase {

    public: StaticStructureFactorGPUEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="ssfq");
        ~StaticStructureFactorGPUEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
        
        std::vector<dVec> qValues;              // Vector of q values
        DynamicArray<dVec,1> qValues_dVec;      // Vector of q values
        DynamicArray<double,1> ssf;             // local intermediate scattering function

        int numq;                        // the number of q std::vectors
        size_t bytes_beads;
        size_t bytes_ssf;
        size_t bytes_qvecs;
        double *d_beads;                // pointer to beads on gpu (device_beads)
        double *d_ssf;                  // pointer to ssf on gpu (device_ssf)
        double *d_qvecs;                // pointer to qvecs on gpu (device_qvecs)

        //FIXME stream handling needs to be moved out of estimators
        gpu_stream_t stream_array[MAX_GPU_STREAMS]; // Store Multiple GPU streams
};
#endif

#ifdef USE_GPU
// ========================================================================  
// GPU Accellerated Cylinder Static Structure Factor Estimator Class
// ========================================================================  
/** 
 * Compute the static structure factor S(q)
 */
class CylinderStaticStructureFactorGPUEstimator: public EstimatorBase {

    public: CylinderStaticStructureFactorGPUEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="cyl_ssfq");
        ~CylinderStaticStructureFactorGPUEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}
        void sample();              // Sample the estimator


    private:
        void accumulate();              // Accumulate values
        
        std::vector<dVec> qValues;              // Vector of q values
        DynamicArray<dVec,1> qValues_dVec;      // Vector of q values
        DynamicArray<double,1> ssf;             // local intermediate scattering function

        int numq;                        // the number of q std::vectors
        size_t bytes_beads;
        size_t bytes_ssf;
        size_t bytes_qvecs;
        double *d_beads;                // pointer to beads on gpu (device_beads)
        double *d_ssf;                  // pointer to ssf on gpu (device_ssf)
        double *d_qvecs;                // pointer to qvecs on gpu (device_qvecs)

        //FIXME stream handling needs to be moved out of estimators
        gpu_stream_t stream_array[MAX_GPU_STREAMS]; // Store Multiple GPU streams
};
#endif

// ========================================================================  
// Intermediate Scattering Function Estimator Class
// ========================================================================  
/** 
 * Compute the intermediate scattering function F(q,\tau)
 */
class IntermediateScatteringFunctionEstimator: public EstimatorBase {

    public:
        IntermediateScatteringFunctionEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="isf");
        ~IntermediateScatteringFunctionEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();                  // Accumulate values
        DynamicArray <double,1> isf;        // local intermediate scattering function
        std::vector<dVec> qValues;          // Vector of q values
        DynamicArray<dVec, 1> qValues_dVec; // Vector of q values

        int numq;                           // the number of q-magnitudes
        DynamicArray <int,1> numqVecs;      // the number of q-vectors with a given magnitude
        std::vector <std::vector<dVec> > q; // the q-vectors
};

#ifdef USE_GPU
// ========================================================================  
// Intermediate Scattering Function GPU Estimator Class
// ========================================================================  
/** 
 * Compute the intermediate scattering function F(q,\tau)
 */
class IntermediateScatteringFunctionEstimatorGpu: public EstimatorBase {

    public:
        IntermediateScatteringFunctionEstimatorGpu(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="isf");
        ~IntermediateScatteringFunctionEstimatorGpu();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
        std::vector<dVec> qValues;      // Vector of q values
    DynamicArray<dVec,1> qValues_dVec;     // Vector of q values
    DynamicArray<double,1> isf;           // local intermediate scattering function

        int numq;                        // the number of q std::vectors
        size_t bytes_beads;
        size_t bytes_isf;
        size_t bytes_qvecs;
        double *d_beads; // pointer to beads on gpu (device_beads)
        double *d_isf; // pointer to isf on gpu (device_isf)
        double *d_qvecs; // pointer to qvecs on gpu (device_qvecs)

        //FIXME stream handling needs to be moved out of estimators
        gpu_stream_t stream_array[MAX_GPU_STREAMS]; // Store Multiple GPU streams
        
};
#endif

#ifdef USE_GPU
// ========================================================================  
// Elastic Scattering GPU Estimator Class
// ========================================================================  
/** 
 * Compute the elastic scattering S(q, \omega = 0) //FIXME is this true?
 */
class ElasticScatteringEstimatorGpu: public EstimatorBase {

    public:
        ElasticScatteringEstimatorGpu(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="es");
        ~ElasticScatteringEstimatorGpu();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
        std::vector<dVec> qValues;      // Vector of q values
    DynamicArray<dVec,1> qValues_dVec;     // Vector of q values
    DynamicArray<double,1> es;           // local intermediate scattering function

        int numq;                        // the number of q std::vectors
        size_t bytes_beads;
        size_t bytes_es;
        size_t bytes_qvecs;
        double *d_beads; // pointer to beads on gpu (device_beads)
        double *d_es; // pointer to es on gpu (device_es)
        double *d_qvecs; // pointer to qvecs on gpu (device_qvecs)

        //FIXME stream handling needs to be moved out of estimators
        gpu_stream_t stream_array[MAX_GPU_STREAMS]; // Store Multiple GPU streams
};
#endif


// ========================================================================  
// Radial Density Estimator Class
// ========================================================================  
/** 
 * Compute the density as a function of position in the radial direction.
 */
class RadialDensityEstimator: public EstimatorBase {

    public:
        RadialDensityEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="radial");
        ~RadialDensityEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
        double dR;                      // The discretization
};

// ========================================================================  
// Worm Properties Estimator Class
// ========================================================================  
/** 
 * Compute various properties related to the worm in the simulation.
 */
class WormPropertiesEstimator: public EstimatorBase {

    public:
        WormPropertiesEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="worm");
        ~WormPropertiesEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        //dVec sep;                       // head-tail separation
        void accumulate();              // Accumulate values
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
        CylinderEnergyEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_estimator");
        ~CylinderEnergyEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:

        void accumulate();      // Accumulate values

};

// ========================================================================  
// Cylinder Number Particles Estimator Class 
// ========================================================================  
/**
 * Computes the average number of particles, as well as density.
 */
class CylinderNumberParticlesEstimator: public EstimatorBase {

    public:
        CylinderNumberParticlesEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_estimator");
        ~CylinderNumberParticlesEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
};

// ========================================================================  
// Cylinder Number Distribution Estimator Class 
// ========================================================================  
/**
 * Computes the probability distribution function for the number of particles.
 */
class CylinderNumberDistributionEstimator: public EstimatorBase {

    public:
        CylinderNumberDistributionEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_number");
        ~CylinderNumberDistributionEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        int maxNumParticles;        // The maximum number considered

        void accumulate();          // Accumulate values
};

#if NDIM > 1
// ========================================================================  
// Cylinder Linear Density Estimator Class 
// ========================================================================  
/**
 * Computes the density as a function of distance along the cylinder axis.
 */
class CylinderLinearDensityEstimator: public EstimatorBase {

    public:
        CylinderLinearDensityEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_linedensity");
        ~CylinderLinearDensityEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        double dz;                  // The bin-size in the z-direction
        double Lz;                  // The length of the cylinder

        void accumulate();          // Accumulate values
};
#endif

// ========================================================================  
// Cylinder Superfluid Fraction Estimator Class 
// ========================================================================  
/**
 * Compute the superfluid fraction, as well as the winding number 
 * probability distribution.
 */
class CylinderSuperfluidFractionEstimator: public EstimatorBase {

    public:
        CylinderSuperfluidFractionEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_super");
        ~CylinderSuperfluidFractionEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
    DynamicArray <bool,1> doBead;  // Used for ensuring we don't double count beads
        int windMax;            // The maximum winding number considered

        void accumulate();      // Accumulate values

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
        CylinderOneBodyDensityMatrixEstimator(Path &, ActionBase *, const MTRand &, 
                double, int _frequency=20, std::string _label="cyl_obdm");
        ~CylinderOneBodyDensityMatrixEstimator();

        void sample();              // Sample the estimator

        static const std::string name;
        std::string getName() const {return name;}

    private:
        Path &lpath;                    // A non-constant local reference to the path

        double dR;                      // The discretization
        int numReps;                    // The number of measurments reps                   
        uint32 numAccepted;             // The number of moves accepted
        uint32 numAttempted;            // The number of moves attempted


        dVec newTailPos,oldTailPos;     // The new and old tail position
        //dVec newHeadPos;                // The new head position
        dVec newRanPos,neighborPos;     // The random shift

        double sqrt2LambdaTau;          // sqrt(2 * lambda * tau)
        double rho0Norm;                // Free density matrix
        double oldAction,newAction;     // The old and new action

        /* Get a random std::vector */
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
        CylinderPairCorrelationEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_pair");
        ~CylinderPairCorrelationEstimator();
    
        void sample();              // Sample the estimator

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();              // Accumulate values
        double dR;                      // The discretization
};

#if NDIM > 1
// ========================================================================  
// Cylinder Linear Potential Estimator Class
// ========================================================================  
/** 
 * Compute the effective linear potential along the axis of the cylinder.
 */
class CylinderLinearPotentialEstimator: public EstimatorBase {

    public:
        CylinderLinearPotentialEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_linepotential");
        ~CylinderLinearPotentialEstimator();

        static const std::string name;
        std::string getName() const {return name;}

        void output() {outputHist();}  // overload the output

    private:
        
        double dz;                      // The discretization
        double Lz;                      // The length of the cylinder

        void accumulate();              // Accumulate values
        void accumulate1();              // Accumulate values
};
#endif

#if NDIM > 2
// ========================================================================  
// Cylinder Radial Potential Estimator Class
// ========================================================================  
/** 
 * Compute the effective radial potential in a cylinder.
 */
class CylinderRadialPotentialEstimator: public EstimatorBase {

    public:
        CylinderRadialPotentialEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="cyl_potential");
        ~CylinderRadialPotentialEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        
        double dR;                      // The discretization
    DynamicArray <double,1> radPot;        // Used for normalization

        void accumulate();              // Accumulate values
        void accumulate1();             // Accumulate values
};
#endif

// ========================================================================  
// Cylinder Static Structure Factor Estimator Class
// ========================================================================  
/** 
 * Compute the static structure factor S(q)
 */
class CylinderStaticStructureFactorEstimator: public EstimatorBase {

    public:
        CylinderStaticStructureFactorEstimator(const Path &, ActionBase *, 
                const MTRand &, double, int _frequency=1, std::string _label="cyl_ssf");
        ~CylinderStaticStructureFactorEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}
        void sample();              // Sample the estimator

    private:
        void accumulate();              // Accumulate values
    DynamicArray <double,1> sf;            // structure factor
        std::vector <std::vector<dVec> > q;       // the q-vectors
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
        PotentialEnergyEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="potential");
        ~PotentialEnergyEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values

};

// ========================================================================
// Kinetic Energy Estimator Class
// ========================================================================
/**
 * Computes the total energy using a mixed estimator
 */
class KineticEnergyEstimator: public EstimatorBase {
    
    public:
        KineticEnergyEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="kinetic");
        ~KineticEnergyEstimator();
    
        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
    
};

// ========================================================================
// Total Energy Estimator Class
// ========================================================================
/**
 * Computes the total energy using a mixed estimator
 */
class TotalEnergyEstimator: public EstimatorBase {
    
    public:
        TotalEnergyEstimator(const Path &, ActionBase *, const MTRand &, 
                double, int _frequency=1, std::string _label="energy");
        ~TotalEnergyEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values

};

// ========================================================================
// Themodynamic Potential Energy Estimator Class
// ========================================================================
/**
 * Computes the total energy using a mixed estimator
 */
class ThermoPotentialEnergyEstimator: public EstimatorBase {
    
    public:
        ThermoPotentialEnergyEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="thpotential");
        ~ThermoPotentialEnergyEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
    
};

// ========================================================================  
// Position Estimator Class 
// ========================================================================  
/** 
 * Computes the average value of the position in 1D.
 *
 */
class PositionEstimator: public EstimatorBase {

    public:
        PositionEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="position");
        ~PositionEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values

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
        ParticleResolvedPositionEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="prposition");
        ~ParticleResolvedPositionEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values
    
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
        ParticleCorrelationEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="prcorrelation");
        ~ParticleCorrelationEstimator();

        static const std::string name;
        std::string getName() const {return name;}


    private:
        void accumulate();      // Accumulate values

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
        VelocityEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="velocity");
        ~VelocityEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values

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
        SubregionOccupationEstimator(const Path &, ActionBase *, const MTRand &,
                double, int _frequency=1, std::string _label="subregionocc");
        ~SubregionOccupationEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        void accumulate();      // Accumulate values

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
                double, int _frequency=1, std::string _label="obdm");
        ~PIGSOneBodyDensityMatrixEstimator();

        void sample();              // Sample the estimator
        void outputFooter();        // Output the acceptance footer to disk

        static const std::string name;
        std::string getName() const {return name;}

    private:
        Path &lpath;                    // A non-constant local reference to the path

        double dR;                      // The discretization
        int numReps;                    // The number of measurments reps
        uint32 numAccepted;             // The number of moves accepted
        uint32 numAttempted;            // The number of moves attempted


        dVec newTailPos,oldTailPos;     // The new and old tail position
        //dVec newRanPos,neighborPos;     // The random shift

        double sqrt2LambdaTau;          // sqrt(2 * lambda * tau)
        double oldAction,newAction;     // The old and new action

        /* Get a random std::vector */
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
        DoubledEstimator(const Path &, const Path &, ActionBase*, ActionBase*, 
                const MTRand &, double, int _frequency=1, std::string _label="");
        ~DoubledEstimator();

        std::string getName() const {return "doubled base";}

    protected:
        const Path &path2;              ///< A constant reference to the paths
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
        SwapEstimator(Path &, Path &, ActionBase *, ActionBase *,
                const MTRand &, double, int _frequency=1, std::string _label="swap");
        ~SwapEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        Path &lpath;                    // A non-constant local reference to path 1
        Path &lpath2;                   // A non-constant local reference to path 2
        ActionBase *actionPtr;          // The action pointer for path 1
        ActionBase *actionPtr2;         // The action pointer for path 2

        void accumulate();          // Accumulate values
        void accumulateOpen();          // Accumulate values for open paths
        void accumulateClosed();            // Accumulate values for open paths
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
                const MTRand &, double, int _frequency=1, std::string _label="entpart");
        ~EntPartEstimator();

        static const std::string name;
        std::string getName() const {return name;}

    private:
        Path &lpath;                    // A non-constant local reference to path 1
        Path &lpath2;                   // A non-constant local reference to path 2
        ActionBase *actionPtr;          // The action pointer for path 1
        ActionBase *actionPtr2;         // The action pointer for path 2

        void accumulate();          // Accumulate values
};



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// END PIGS ESTIMATORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
#endif
