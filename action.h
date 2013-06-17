/**
 * @file action.h
 * @author Adrian Del Maestro
 * @date 10.14.2008
 *
 * @brief Action class definitions.
 */

#include "pimc.h"
#include "constants.h"

#ifndef ACTION_H 
#define ACTION_H

class Path;
class PotentialBase;
class LookupTable;

// ========================================================================  
// ActionBase Class
// ========================================================================  
/** 
 * Holds a base class that all action classes will be derived from. 
 *
 * Implements the details of the action, includign the potential and kinetic
 * pieces.  Two types of actions, local and non-local derive from this base 
 * class and the actually used actions then derive from these.
 */
class ActionBase {

	public:
		ActionBase (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, bool _local=true, string _name="Base");
		virtual ~ActionBase();

		/** Returns the action name */
		string getActionName () { return name; }

		/* The kinetic Action  */
		double kineticAction ();
		double kineticAction (const beadLocator &);
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

        /* Various derivatives of the potential action with a cutoff */
		virtual double derivPotentialActionTau (int,double) { return 0.0; }
		virtual double derivPotentialActionLambda (int,double) { return 0.0; }

        /* The bare local potential at a single time slice */
        virtual double potential(int) { return 0.0; }
        virtual double potential(int,double) { return 0.0; }

		/** The public method that sets the tau scaling factor. */
		void setShift(int _shift) { shift = _shift; }
		/** Get the tau scaling factor. */
		int getShift() { return shift; }

		/** The free-particle density matrix */
		double rho0(const dVec&, const dVec&, int); 
		double rho0(const beadLocator&, const beadLocator&, int); 
		
		/** The ensemble particle number weighting factor */
		double ensembleWeight(const int);

        /** Is the action local in imaginary time? */
        const bool local;

		PotentialBase *externalPtr;		///< The external potential
		PotentialBase *interactionPtr;	///< The interaction potential

		Array <int,1> sepHist;		    ///< A histogram of separations
		Array <int,1> cylSepHist;	    ///< A histogram of separations for a cylinder

	protected:
        string name;                    ///< The name of the action

		LookupTable &lookup;			///< We need a non-constant reference for updates
		const Path &path;				///< A reference to the paths

		int shift;						///< The scaling factor for tau

		TinyVector <double,2> eFactor;	///< The even/odd slice energy factor

		/* These definitions are needed for weighting in the canonical ensemble */
		bool canonical;					///< Are we in the canonical ensemble?
		int numBeads0;					///< The target number of beads
		double deltaNumBeads2;			///< The fluctuating weight
        bool window;                    // Whether or not to force particle numberin to window
        int windowWidth;                // Half width of particle number window
        bool gaussianEnsemble;          // Whether or not to use gaussian ensemble weight
        double gaussianEnsembleSD;        // Standard deviation of guassian ensemble weight

		/** The local shifted value of tau. */
		double tau() {return shift * constants()->tau();}

		beadLocator bead2,bead3;		// Bead indexers
		dVec sep,sep2;				    // The spatial separation between beads.
		double dSep;					// The discretization for the separation histogram

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
                PotentialBase *, string _name="Local");
		virtual ~LocalAction();

        /* The potential action */
		double potentialAction ();
		//double potentialAction (const beadLocator &, const beadLocator &);
		double potentialAction (const beadLocator &);

        /* The bare potential action and its correction */
		double barePotentialAction (const beadLocator &);
		double potentialActionCorrection (const beadLocator &);
		double potentialActionCorrection (const beadLocator &, const beadLocator &);

        /* Various derivatives of the potential action */
		double derivPotentialActionTau (int);
		double derivPotentialActionLambda (int);
		double derivPotentialActionTau (int,double);
		double derivPotentialActionLambda (int,double);

        /* The bare local potential at a single time slice */
        virtual double potential(int slice) { return V(slice); }
        virtual double potential(int slice, double maxR) { return V(slice,maxR); }

    protected:
		int eo;							///< Is a slice even or odd?

		TinyVector <double,2> VFactor;	    ///< The even/odd slice potential factor
		TinyVector <double,2> gradVFactor;	///< The even/odd slice correction factor

        /* Short-hands for potential action and correction terms */
        TinyVector <double,2> potentialFactor;  ///< The even/odd slice total potential factor
        TinyVector <double,2> fFactor;  ///< The even/odd slice force correction potential factor

        /* Setup the time-saving action factors */
        virtual void setFactor();

		/* The full potential for a single bead and all beads at a single
		 * time slice. */
		double V(const beadLocator&);	
		double V(const int);
		double V(const int, const double);

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
};

// ========================================================================  
// PrimitiveAction Class
// ========================================================================  
/** 
 * The basic primitive action, accurate to order O(\tau^2)
 */
class PrimitiveAction: public LocalAction {

	public:
		PrimitiveAction (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, string _name="Primitive");
		~PrimitiveAction();
};

// ========================================================================  
// LiBroughtonAction Class
// ========================================================================  
/** 
 * The Li-Broughton action which should be valid up to order tau^4.
 */
class LiBroughtonAction : public LocalAction {

	public:
		LiBroughtonAction (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, string _name="Li-Broughton");
		~LiBroughtonAction();
};

// ========================================================================  
// GSFAction Class
// ========================================================================  
/**
 * Implementation of the Generalized Suzuki Factorization action
 * accurate up to order tau^5. 
 * @see J. Chem. Phys. 115, 7832 (2001) for details of the implementation.
 */
class GSFAction : public LocalAction {

	public:
		GSFAction (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, string _name="Generalized Suzuki Factorization");
		~GSFAction();

	private:
		double alpha;			// The parameter of the GSF algorithm

};

// ========================================================================  
// NonLocalAction Class
// ========================================================================  
/** 
 * A base class to be inherited by actions that are non-local in imaginary 
 * time.
 *
 */
class NonLocalAction : public ActionBase {

    public:
		NonLocalAction (const Path &, LookupTable &, PotentialBase *, 
                PotentialBase *, string _name="Non-local");
        virtual ~NonLocalAction();

        /* The potential Action */
		double potentialAction ();
        double potentialAction (const beadLocator &);

        /* Derivatives of the potential action */
		double derivPotentialActionTau (int);
		double derivPotentialActionLambda (int);

    protected:
		double U(int);
};

// ========================================================================  
// PairProductAction Class
// ========================================================================  
/** 
 * The Pair Product Action without any external potential.
 *
 * @see Section IV.F in D. M. Ceperley, Rev. Mod. Phys. 67, 279–355 (1995).
 */
//class PairProductAction: public NonLocalAction {
//
//	public:
//		PairProductAction(const Path &, string _name="Pair Product Approximation"); 
//		~PairProductAction();
//};

#endif
