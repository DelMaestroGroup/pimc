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
#define ACtION_H

class Path;
class Potential;
class LookupTable;

// ========================================================================  
// ActionBase Class
// ========================================================================  
/** 
 * Holds a base class that all action classes will be derived from. 
 *
 * The free action, which includes no corrections invovling derivatives of 
 * the gradient of the potential.  Accurate up to tau^2.
 */
class ActionBase {

	public:
		ActionBase (const Path &, Potential*);
		virtual ~ActionBase();

		/* The total action */
		double totalAction();

		/** Returns the action name */
		virtual string getActionName () { return "Base Action"; }

		/* The kinetic Action  */
		double kineticAction ();
		double kineticAction (const beadLocator &beadIndex);
		double kineticAction (const beadLocator &beadIndex, int wlLength);

		/* The potential Action */
		virtual double potentialAction();
		virtual double potentialAction (const beadLocator &beadIndex);
		virtual double potentialAction (const beadLocator &beadIndex, int wlLength);

		/** The public method that sets the tau scaling factor. */
		void setShift(int _shift) { shift = _shift; }
		/** Get the tau scaling factor. */
		int getShift() { return shift; }

		/* The free-particle density matrix */
		double rho0(const dVec&, const dVec&, int); 
		double rho0(const beadLocator&, const beadLocator&, int); 
		
		/* The ensemble particle number weighting factor */
		double ensembleWeight(const int);

	protected:
		friend class EnergyEstimator;				// There are action corrections to the energy
		friend class OneBodyDensityMatrixEstimator;	// OBDM needs eFactor
		friend class PairCorrelationEstimator;		// Need to get at the separation histogram

		friend class CylinderEnergyEstimator;				// There are action corrections to the energy
		friend class CylinderOneBodyDensityMatrixEstimator;	// OBDM needs eFactor
		friend class CylinderPairCorrelationEstimator;		// Need to get at the separation histogram
		friend class CylinderRadialPotentialEstimator;		// Needed to compute the effective radial potential

		const Path &path;				///< A reference to the paths
		Potential *potentialPtr;		///< A pointer to a potential object

		int shift;						///< The scaling factor for tau
		int eo;							///< Is a slice even or odd?

		TinyVector <double,2> pFactor;	///< The even/odd slice potential factor
		TinyVector <double,2> cFactor;	///< The even/odd slice correction factor
		TinyVector <double,2> eFactor;	///< The even/odd slice energy factor

		/* These definitions are needed for weighting in the canonical ensemble */
		bool canonical;					///< Are we in the canonical ensemble?
		int numBeads0;					///< The target number of beads
		double deltaNumBeads2;			///< The fluctuating weight


		/** The local shifted value of tau. */
		double tau() {return shift * constants()->tau();}

};

// ========================================================================  
// LiBroughtonAction Class
// ========================================================================  
/** 
 * The Li-Broughton action which should be valid up to order tau^4.
 */
class LiBroughtonAction : public ActionBase {

	public:
		LiBroughtonAction(const Path &, Potential*); 
		~LiBroughtonAction();

		/** Return the action name */
		string getActionName () { return "Li-Broughton Action"; }
};

// ========================================================================  
// GSFAction Class
// ========================================================================  
/**
 * Implementation of the Generalized Suzuki Factorization action
 * accurate up to order tau^5. 
 * @see J. Chem. Phys. 115, 7832 (2001) for details of the implementation.
 */
class GSFAction : public ActionBase {

	public:
		GSFAction(const Path &, Potential*);
		~GSFAction();

		/** Return the action name */
		string getActionName () { return "Generalized Suzuki Factorization Action"; }

	private:
		double alpha;			// The parameter of the GSF algorithm

};

#endif
