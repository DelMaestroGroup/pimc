/**
 * @file potential.h
 * @author Adrian Del Maestro
 * @date 10.14.2008
 *
 * @brief All possible potential classes.
 */

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "common.h"
#include "constants.h"

class Path;
class LookupTable;
class Container;

// ========================================================================  
// PotentialBase Class
// ========================================================================  
/** 
 * The base class from which all specific potentials are derived from.
 *
 * This class contains methods which return the actual value of 
 * the potential, an effective potential related to the pair product
 * approximation,  one which returns the gradient of the potential and a final
 * one which generates a sensible initial particle configuration.
 * We require knowledge of both the interaction as well as external potential
 * to run the simulation.
 */
class PotentialBase {

	public:
		PotentialBase ();
		virtual ~PotentialBase();
	
		/** The potential */
		virtual double V(const dVec &r) { return 0.0; }

		/** The effective potential for the pair product approximation */
		virtual double V(const dVec &sep1, const dVec &sep2, double lamtau) { return 0.0; }

		/** The gradient of the potential*/
		virtual dVec gradV(const dVec &r) { return 0.0; }

        /** The derivative of the effective potential with respect to lambda
         *  and tau */
        virtual double dVdlambda(const dVec &sep1, const dVec &sep2, double lambda, double tau) {return 0.0;}
        virtual double dVdtau(const dVec &sep1, const dVec &sep2, double lambda, double tau) {return 0.0;}
		
		/** Default Initial configuration of particles*/
		virtual Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

		/** A debug method that output's the potential to a supplied separation */
		void output(const double);

		double tailV;		///< Tail correction factor.
};

// ========================================================================  
// TabulatedPotential Class
// ========================================================================  
/** 
 * Pre-tabulated potential for complicated functions.
 *
 * In order to speed up the evaluation of complicated potentials, we
 * use a 2-point Newton-Gregory spline fit to perform the actual 
 * interpolation.
 */
class TabulatedPotential {
	public:
		TabulatedPotential();
		virtual ~TabulatedPotential();

	protected:
		Array <double,1> lookupV;			///< A potential lookup table
		Array <double,1> lookupdVdr;		///< A lookup table for dVint/dr

		double dr;							///< The discretization for the lookup table
		int tableLength;					///< The number of elements in the lookup table

		TinyVector<double,2> extV;			///< Extremal value of V
		TinyVector<double,2> extdVdr;		///< Extremal value of dV/dr

		/* Initialize all data structures */
		void initLookupTable(const double, const double);

		/* Returns the 2-point spline fit to the lookup table */
		virtual double newtonGregory(const Array<double,1>&, const TinyVector<double,2>&, const double);

        /* Returns a bare lookup value */
		virtual double direct(const Array<double,1>&, const TinyVector<double,2>&, const double);

		/** The functional value of V */
		virtual double valueV (const double) = 0;				
		/** The functional value of dV/dr */
		virtual double valuedVdr (const double) = 0;					
};

// ========================================================================  
// FreePotential Class
// ========================================================================  
/** 
 * Free potential.
 */
class FreePotential: public PotentialBase {
	public:
		FreePotential();
		~FreePotential();

		/** The potential. */
		double V(const dVec &sep) { return 0.0*sep[0]; };

		/** The gradient of the potential. */
		dVec gradV(const dVec &pos) {
			return (0.0*pos);
		}
};

// ========================================================================  
// HarmonicPotential Class
// ========================================================================  
/** 
 * Computes the potential energy for an external harmonic potential.  
 *
 * We work in generalized units where hbar omega / k_B = 1.
 */
class HarmonicPotential : public PotentialBase {
	public:
		HarmonicPotential ();
		~HarmonicPotential ();

		/** The potential. */
		double V(const dVec &r) { 
			return (dot(r,r)/(4.0*constants()->lambda()));
		}

		/** The gradient of the potential. */
		dVec gradV(const dVec &r) {
			dVec tempr;
			tempr = r;
			return (tempr/(2.0*constants()->lambda()));
		}

		/** Initial configuration corresponding to Harmonic potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 
};

// ========================================================================  
// SingleWellPotential Class
// ========================================================================  
/**
 * Computes the potential energy for an external single well potential.
 */
class SingleWellPotential : public PotentialBase {
	public:
		SingleWellPotential ();
		~SingleWellPotential ();

		/** The potential */
		double V(const dVec &r) { 
			double r2 = dot(r,r);
			return ( 0.5*r2 + r2*r2 );
		}

		/** The gradient of the potential. */
		dVec gradV(const dVec &r) {
			double r2 = dot(r,r);
			dVec tempr;
			tempr = r;
			return ((1.0 + 4.0*r2)*tempr);
		}
};

// ========================================================================  
// HarmonicCylinderPotential Class
// ========================================================================  
/**
 * Computes the potential energy for an external harmonic potential with
 * axial symmetry.
 */
class HarmonicCylinderPotential : public PotentialBase {
	public:
		HarmonicCylinderPotential (const double);
		~HarmonicCylinderPotential ();

		/** The potential. */
		double V(const dVec &r) { 
			double r2 = 0.0;
			for (int i=0; i < NDIM-1; i++)
				r2 += r[i]*r[i];
			return ( 0.5 * c * constants()->m() * w * w * r2 ); 
		}

		/** The gradient of the potential. */
		dVec gradV(const dVec &r) {
			dVec tempr;
			tempr = 0.0;
			for (int i=0; i < NDIM-1; i++)
				tempr[i] = r[i];
			return ( c * constants()->m() * w * w * tempr );
		}

	private:
		double w;				// The confining frequency
		double c;				// A dimension-full constant
};

// ========================================================================  
// DeltaPotential Class
// ========================================================================  
/** 
 * Computes the potential energy for delta function interaction potential,
 * approximated here as the limit of a Cauchy distribution.
 */
class DeltaPotential : public PotentialBase  {
	public:
		DeltaPotential (double,double);
		~DeltaPotential ();

		/**
		 * Return the delta function potential with strength 2c approximated as
		 * the limit of a Guassian distribution.
		 * Tested and working in Mathematica.
		 */
		double V(const dVec &r) {
			return (norm*exp(-dot(r,r)/(4.0*a)));
		}

		/**
		 * Return the gradient of the delta function potential with strength 
		 * 2c approximated as the limit of a Guassian distribution.
		 * Tested and working in Mathematica.
		 */
		dVec gradV(const dVec &r) {
			return (-r*norm*exp(-dot(r,r)/(4.0*a))/(2.0*a));
		}

	private:
		double c;				// The strength of the delta function
		double norm;			// A normalization constant for fixed strength
		double a;				// The order of the limit
};

// ========================================================================  
// LorentzianPotential Class
// ========================================================================  
/**
 * Computes the potential energy for delta function interaction potential,
 * approximated here as the limit of a Cauchy distribution.
 */
class LorentzianPotential : public PotentialBase  {
	public:
		LorentzianPotential (double,double);
		~LorentzianPotential ();

		/** 
		 * Return the delta function potential with strength 2c approximated as
		 * the limit of a Lorentzian distribution.
		 * Tested and working in Mathematica.
		 */
		double V(const dVec &r) {
			return (norm / (a*a + dot(r,r)));
		}

		/** 
		 * Return the gradient of the delta function potential with strength 
		 * 2c approximated as the limit of a Lorentzian distribution.
		 * Tested and working in Mathematica.
		 */
		dVec gradV(const dVec &r) {
			double b = a*a + dot(r,r);
			return ((-(2.0*norm*a)/(b*b))*r);
		}

	private:
		double c;				// The strength of the delta function
		double norm;			// A normalization constant for fixed strength
		double a;				// The order of the limit
};


// ========================================================================  
// Hard Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a hard-walled 
 * cylindrical cavity.
 */
class HardCylinderPotential : public PotentialBase {
	public:
		HardCylinderPotential (const double);
		~HardCylinderPotential ();

		/** A step function at rho=R. */
		double V(const dVec &r) {
			if (sqrt(r[0]*r[0]+r[1]*r[1]) >= R)
				return BIG;
			else
				return 0.0;
		}

		/** A delta function at rho=R. */
		dVec gradV(const dVec &r) {
			dVec tempr;
			tempr = r;
			tempr[2] = 0.0;
			if (abs((r[0]*r[0]+r[1]*r[1])-R*R)<1.0E-3)
				return BIG*tempr;
			else
				return 0.0*tempr;
		}

	private:
		double R;		// Radius of the tube
};

// ========================================================================  
// LJ Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a cylindrical 
 * cavity. 
 */
class LJCylinderPotential : public PotentialBase, public TabulatedPotential {
	public:
		LJCylinderPotential (const double);
		~LJCylinderPotential ();

		/** The integrated LJ Wall potential. */
		double V(const dVec &r) {
			int k = static_cast<int>(sqrt(r[0]*r[0] + r[1]*r[1])/dR);
			if (k >= tableLength)
				return extV[1];
			else
				return lookupV(k);
		}

		/* The gradient of the LJ Wall potential */
		dVec gradV(const dVec &);

		/** Initial configuration corresponding to the LJ cylinder potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

	private:
		/* All the parameters needed for the LJ wall potential */
		double density;
		double sigma;
		double epsilon;

		double R;		// Radius of the tube
		double dR;		// Discretization for the lookup table

		double minV;	// The minimum value of the potential

		/* Used to construct the lookup tables */
		double valueV (const double);				
		double valuedVdr (const double);					
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of aziz potential for separation r using a 
 * lookup table. 
 */
inline dVec LJCylinderPotential::gradV(const dVec &r) {
	double rnorm = sqrt(r[0]*r[0] + r[1]*r[1]);
	dVec tempr;
	tempr = r;
	tempr[2] = 0.0;
	int k = static_cast<int>(rnorm/dR);
	dVec gV;
	if (k >= tableLength)
		gV = (extdVdr[1]/rnorm)*tempr;
	else
		gV = (lookupdVdr(k)/rnorm)*tempr;
	return gV;
}

// ========================================================================  
// Aziz Potential Class
// ========================================================================  
/** 
 * Computes the value of the semi-empircal Aziz potential that is known
 * to be accurate for He-4.
 */
class AzizPotential : public PotentialBase, public TabulatedPotential {
	public:
		AzizPotential (const dVec &);
		~AzizPotential ();

		/* The Aziz HFDHE2 Potential */
		double V(const dVec &);

		/* The gradient of the Aziz potential */
		dVec gradV(const dVec &);

	private:
		/* All the parameters of the Aziz potential */
		double rm, A, epsilon, alpha, D, C6, C8, C10;

		/* Used to construct the lookup tables */
		double valueV (const double);				
		double valuedVdr (const double);					

		/* The F-function needed for the Aziz potential */
		double F(const double x) {
			return (x < D ? exp(-(D/x - 1.0)*(D/x - 1.0)) : 1.0 );
		}

		/* The derivative of the F-function needed for the Aziz potential */
		double dF(const double x) {
			double ix = 1.0/x;
			double r = 2.0*D*ix*ix*(D*ix-1.0)*exp(-(D*ix - 1.0)*(D*ix - 1.0));
			return (x < D ? r : 0.0 );
		}
};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the aziz potential for separation r using a lookup table. 
 */
inline double AzizPotential::V(const dVec &r) {
	//double rnorm = sqrt(dot(r,r));
	//return newtonGregory(lookupV,extV,rnorm);
	return direct(lookupV,extV,sqrt(dot(r,r)));
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of aziz potential for separation r using a 
 * lookup table. 
 */
inline dVec AzizPotential::gradV(const dVec &r) {
	double rnorm = sqrt(dot(r,r));
	dVec gV;
	//gV = (newtonGregory(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
	gV = (direct(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
	return gV;
}

// ========================================================================  
// FixedAzizPotential Class
// ========================================================================  
/** 
 * Computes the potential energy resulting from a series of fixed helium
 * atoms that are not updated and provide a static 'external' potential.
 *
 * We require a Aziz potential interaction pointer to properly compute
 * the interaction with the static particles.
 */
class FixedAzizPotential : public PotentialBase  {

	public:
		FixedAzizPotential(const Container *);
		~FixedAzizPotential();

		/* Return the sum of the Aziz 'interaction energy' between the supplied
		 * particle and all fixed particles. */
		double V(const dVec &r);

		/* Return the gradient of the sum of the Aziz 'interaction energy' */
		dVec gradV(const dVec &r);

		/** Initial configuration corresponding to FixedAziz potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

	private:
		AzizPotential aziz;					// A copy of the aziz potential
		Array <dVec,1> fixedParticles;		// The location of the fixed particles
		Array <int,2> fixedBeadsInGrid;		// The local grid hash table
		Array <int,1> numFixedBeadsInGrid;	// The number of fixed particles in each grid box
		int numFixedParticles;				// The total number of fixed particles
		LookupTable *lookupPtr;				// A lookup table pointer
		double rc2;							// A local copy of the potential cutoff squared

};

// ========================================================================  
// Gasparini_1_Potential Class (Primitive Approximation)
// ========================================================================  
/**
 * Computes potential energy for Gasparini potential.
 * This is the case of Primitive Approx. (gradV = 0 for all positions)
 */
class Gasparini_1_Potential : public PotentialBase {
	public:
		Gasparini_1_Potential (double, double, const Container*);
		~Gasparini_1_Potential ();

        /** The potential */
        double V(const dVec &r) { 
            double r2 = 0.0;
            if ((r[2] >= -az) && (r[2] <= az) && (r[1] >= -ay) && (r[1] <= ay))
                r2 = 1.0*V0;
            return r2;
		}

		/** The gradient of the potential. */
		dVec gradV(const dVec &r) {
			return 0.0;
		}
    
        /** Initial configuration corresponding to FixedAziz potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /* parameters needed for Gasp Potential_1 */
        const double az;        //scales the length (z)
        const double ay;        //scales the width (y,x)
        const double V0;       //scales the potential step
};

// ========================================================================  
// Hard Sphere Potential Class
// ========================================================================  
/** 
 * Computes the effective potential from the exact two-body density matrix
 * for hard spheres.  
 *
 * @see: S. Pilati, K. Sakkos, J. Boronat, J. Casulleras, and 
 *       S. Giorgini, Phys Rev A 74, 043621 (2006).
 */
class HardSpherePotential : public PotentialBase  {
	public:
		HardSpherePotential (double);
		~HardSpherePotential ();

		/** The classical potential */
		virtual double V(const dVec &r) { 
            return ((sqrt(dot(r,r)) <= a) ? BIG : 0.0);
        }

		/** The effective potential */
		double V(const dVec &, const dVec &, double);
        double dVdlambda(const dVec &, const dVec &, double, double);
        double dVdtau(const dVec &, const dVec &, double, double);

	private:
		double a;				// The strength of the delta function
};

// ========================================================================  
// Carbon Nanotube Potential Class
// ========================================================================  
/** 
 * The smooth non-corregated version of the helium-carbon nanotube potential.
 * @see http://prb.aps.org/abstract/PRB/v62/i3/p2173_1
 */
//class CarbonNanotubePotential : public PotentialBase, public TabulatedPotential {
//	public:
//		CarbonNanotubePotential(const double);
//		~CarbonNanotubePotential();
//
//		/** The cylindrically symmetric potential. */
//		double V(const dVec &r) {
//			int k = int(sqrt(r[0]*r[0] + r[1]*r[1])/dR);
//			if (k >= tableLength)
//				return extV[1];
//			else
//				return lookupV(k);
//		}
//
//		/* The gradient of the CNT potential */
//		dVec gradV(const dVec &);
//
//		/** Initial configuration corresponding to the CNT potential */
//		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 
//
//	private:
//		/* All the parameters needed for the LJ wall potential */
//		double density;
//		double sigmaHe,sigmaC;
//		double epsilonHe,epsilonC;
//
//		double R;		// Radius of the tube
//		double dR;		// Discretization for the lookup table
//
//		double minV;	// The minimum value of the potential
//
//		/* Used to construct the lookup tables */
//		double valueV (const double);				
//		double valuedVdr (const double);					
//};
#endif
