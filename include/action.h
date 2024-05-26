/**
 * @file action.h
 * @author Adrian Del Maestro
 * @date 10.14.2008
 *
 * @brief Action class definitions.
 */

#include "constants.h"

#ifndef ACTION_H 
#define ACTION_H

class Path;
class PotentialBase;
class LookupTable;
class WaveFunctionBase;

// ========================================================================  
// ActionBase Class
// ========================================================================  
/** 
 * Holds a base class that all action classes will be derived from. 
 *
 * Implements the details of the action, including the potential and kinetic
 * pieces.  Two types of actions, local and non-local derive from this base 
 * class and the actually used actions then derive from these.
 */
class ActionBase {

    public:
        ActionBase (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, WaveFunctionBase *, bool _local=true,
                std::string _name="Base", double _endFactor=1.0, int _period=1);
        virtual ~ActionBase();

        /** Returns the action name */
        std::string getActionName () { return name; }

        /** The full kinetic Action  */
        double kineticAction ();
        /** The kinetic Action at a single slice  */
        double kineticAction (const beadLocator &);
        /** The kinetic Action for wlLength slices */
        double kineticAction (const beadLocator &, int wlLength);

        /** The effective potential inter-ACTION for various pass conditions. */
        virtual double potentialAction() {return 0.0;}
        virtual double potentialAction (const beadLocator &, const beadLocator &);
        virtual double potentialAction (const beadLocator &) { return 0.0; }

        /* The bare potential action and its correction */
        virtual double barePotentialAction (const beadLocator &) { return 0.0; }
        virtual double potentialActionCorrection (const beadLocator &) {return 0.0; }
        virtual double potentialActionCorrection (const beadLocator &, const beadLocator &) { return 0.0; }

        /* Various derivatives of the potential action */
        virtual double derivPotentialActionTau (int) { return 0.0; }
        virtual double derivPotentialActionLambda (int) { return 0.0; }
        virtual double secondderivPotentialActionTau (int) { return 0.0; }

        /* Various derivatives of the potential action with a cutoff */
        virtual double derivPotentialActionTau (int,double) { return 0.0; }
        virtual double derivPotentialActionLambda (int,double) { return 0.0; }
        
        /* gradient of the potential action */
        virtual dVec gradPotentialAction(int) { return 0.0; }

        /* r \dot gradU -> for virial estimator.
         * It is necessary to split the terms in order to compute
         * the specific heat efficiently.*/
        virtual double rDOTgradUterm1(int) { return 0.0; }
        virtual double rDOTgradUterm2(int) { return 0.0; }

        /* (r-R) \dot gradU -> for centroid virial estimator.
         * It is necessary to split the terms in order to compute
         * the specific heat efficiently.*/
        virtual double deltaDOTgradUterm1(int) { return 0.0; }
        virtual double deltaDOTgradUterm2(int) { return 0.0; }

        /* return virial kinetic energy estimator term */
        virtual double virKinCorr(int) { return 0.0; }

        /* The bare local potential at a single time slice */
        virtual blitz::TinyVector<double,2> potential(int) { return blitz::TinyVector<double,2>(0.0); }
        virtual double potential(int,double) { return 0.0; }

        /** The public method that sets the tau scaling factor. */
        void setShift(int _shift) { shift = _shift; }
        
        /** Get the tau scaling factor. */
        int getShift() { return shift; }

        /** The free-particle density matrix */
        double rho0(const dVec&, const dVec&, int); 
        double rho0(const beadLocator&, const beadLocator&, int); 
        double rho0(const dVec&, const int); 
        
        /** The ensemble particle number weighting factor */
        double ensembleWeight(const int);

        /** Is the action local in imaginary time? */
        const bool local;

        /* The period of the action */
        const int period;

        PotentialBase *externalPtr;     ///< The external potential
        PotentialBase *interactionPtr;  ///< The interaction potential

	blitz::Array <int,1> sepHist;          ///< A histogram of separations
	blitz::Array <int,1> cylSepHist;       ///< A histogram of separations for a cylinder

    protected:
        std::string name;                    ///< The name of the action

        LookupTable &lookup;            ///< We need a non-constant reference for updates
        const Path &path;               ///< A reference to the paths
        WaveFunctionBase *waveFunctionPtr;   ///< A pointer to a trial wave function object
        double endFactor;               ///< Mutiplictive factor of the potential action on ends

        int shift;                      ///< The scaling factor for tau

        /* These definitions are needed for weighting in the canonical ensemble */
        bool canonical;                 ///< Are we in the canonical ensemble?
        int numBeads0;                  ///< The target number of beads
        bool window;                    // Whether or not to force particle numberin to window
        int windowWidth;                // Half width of particle number window
        bool gaussianEnsemble;          // Whether or not to use gaussian ensemble weight
        double gaussianEnsembleSD;        // Standard deviation of guassian ensemble weight

        /** The local shifted value of tau. */
        double tau() {return shift * constants()->tau();}

        beadLocator bead2,bead3;        // Bead indexers
        dVec sep,sep2;                  // The spatial separation between beads.
        double dSep;                    // The discretization for the separation histogram
        double dPerSep;                 // The PBC discretization for the separation histogram 

        /* Update the separation histogram */
        void updateSepHist(const dVec &);   
};

// ========================================================================  
// LocalAction Class
// ========================================================================  
/** 
 * A base class to be inherited by actions that are local in imaginary time.
 *
 */
class LocalAction : public ActionBase {

    public:
        LocalAction (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, WaveFunctionBase *, const blitz::TinyVector<double,2>&, 
                const blitz::TinyVector<double,2>&, bool _local=true, std::string _name="Local",
                double _endFactor=1.0, int _period=1);
        virtual ~LocalAction();

        /* The potential action */
        double potentialAction ();
        double potentialAction (const beadLocator &);

        /* The bare potential action and its correction */
        double barePotentialAction (const beadLocator &);
        double potentialActionCorrection (const beadLocator &);
        double potentialActionCorrection (const beadLocator &, const beadLocator &);

        /* Various derivatives of the potential action */
        double derivPotentialActionTau (int);
        double derivPotentialActionLambda (int);
        double secondderivPotentialActionTau (int);
        double derivPotentialActionTau (int,double);
        double derivPotentialActionLambda (int,double);
 
        /* gradient of the potential action */
        virtual dVec gradPotentialAction(int slice) { return gradU(slice); }

        /* R \dot gradU -> for virial estimators */
        double rDOTgradUterm1(int); 
        double rDOTgradUterm2(int); 

        /* delta \dot gradU -> for centroid virial estimators */
        virtual double deltaDOTgradUterm1(int slice) { return deltadotgradUterm1(slice); }
        virtual double deltaDOTgradUterm2(int slice) { return deltadotgradUterm2(slice); }

        /* Returns virial kinetic energy correction term. */
        virtual double virKinCorr(int slice) { return virialKinCorrection(slice); }

        /* The bare local potential at a single time slice */
        virtual blitz::TinyVector<double,2> potential(int slice) { return V(slice); }
        virtual double potential(int slice, double maxR) { return V(slice,maxR); }

    protected:
        int eo;                         ///< Is a slice even or odd?

	blitz::TinyVector <double,2> VFactor;      ///< The even/odd slice potential factor
	blitz::TinyVector <double,2> gradVFactor;  ///< The even/odd slice correction factor

        /* The full potential for a single bead and all beads at a single
         * time slice. */
        double V(const beadLocator&);   
        double V(const int, const double);

        /* For the potential at a given time slice we separate the interaction
         * and potential parts */
	blitz::TinyVector <double,2> V(const int);

        /* The gradient of the potential squared for a single bead and all beads
         * at a single time slice. */
        double gradVSquared(const beadLocator&);    
        double gradVSquared(const int);
        double gradVSquared(const int, const double);

        /* The full potential with the NN lookup table for a single bead and all
         * beads at a single time slice. */
        double Vnn(const beadLocator&);
        double Vnn(const int);

        /* The bare potential action for a trajectory */
        double bareUnn(const beadLocator &, const beadLocator &);

        /* The gradient of the potential squared for a single bead using the 
         * nearest neighbor lookup table */
        double gradVnnSquared(const beadLocator&);

        /* gradient and Laplacian of potential energy */
        dVec gradientV(const int);

        /* T-matrix needed for grad((gradV)^2) */
        dMat tMatrix(const int);

        /* r \dot gradU -> for virial estimator */
        double RdotgradUterm1(const int);
        double RdotgradUterm2(const int);

        /* delta \dot gradU -> for centroid virial estimator */
        double deltadotgradUterm1(const int);
        double deltadotgradUterm2(const int);

        /* virial kinetic energy term */
        double virialKinCorrection(const int);

        /* gradient of potential action */
        dVec gradU(const int);

};

// ========================================================================  
// NonLocalAction Class
// ========================================================================  
/** 
 * A base class to be inherited by actions that are non-local in imaginary 
 * time.
 *
 * @see Section IV.F in D. M. Ceperley, Rev. Mod. Phys. 67, 279–355 (1995).
 */
class NonLocalAction : public ActionBase {

    public:
        NonLocalAction (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, WaveFunctionBase *, bool _local=false,
                std::string _name="Non-local");
        virtual ~NonLocalAction();

        /* The potential Action */
        double potentialAction ();
        double potentialAction (const beadLocator &);

        /* Derivatives of the potential action */
        double derivPotentialActionTau (int);
        double derivPotentialActionLambda (int);
    
        /* The bare local external and interaction potential at a single time slice */
        virtual blitz::TinyVector<double,2> potential(int slice) { return U(slice); }

    protected:
	blitz::TinyVector<double,2> U(int);
    
    private:
        std::vector<bool> NNbead; //Records which beads where already visited for NN operation
};
#endif
