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
		virtual double V(const dVec &) { return 0.0; }

		/** The effective potential for the pair product approximation */
		virtual double V(const dVec &, const dVec &) { return 0.0; }

		/** The gradient of the potential*/
		virtual dVec gradV(const dVec &) { return 0.0; }

        /** Grad^2 of the potential*/
		virtual double grad2V(const dVec &) { return 0.0; }

        /** The derivative of the effective potential with respect to lambda
         *  and tau */
        virtual double dVdlambda(const dVec &, const dVec &) {return 0.0;}
        virtual double dVdtau(const dVec &, const dVec &) {return 0.0;}
		
		/** Default Initial configuration of particles*/
		virtual Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

		/** A debug method that output's the potential to a supplied separation */
		void output(const double);

		double tailV;		///< Tail correction factor.

        /** Array to hold data elements*/
		virtual Array<double,1> getExcLen();

    protected:
        double deltaSeparation(double sep1,double sep2) const;
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
		Array <double,1> lookupd2Vdr2;		///< A lookup table for d2Vint/dr2

		double dr;							///< The discretization for the lookup table
		int tableLength;					///< The number of elements in the lookup table

		TinyVector<double,2> extV;			///< Extremal value of V
		TinyVector<double,2> extdVdr;		///< Extremal value of dV/dr
		TinyVector<double,2> extd2Vdr2;		///< Extremal value of d2V/dr2

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

        /** The functional value of d2V/dr2 */
		virtual double valued2Vdr2 (const double) = 0;					
	
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
        HarmonicPotential (double);
		~HarmonicPotential ();
    
        double omega2;       //The SHO frequency in units of hbar

		/** The potential. */
		double V(const dVec &r) { 
			return (omega2*dot(r,r)/(4.0*constants()->lambda()));
		}

		/** The gradient of the potential. */
		dVec gradV(const dVec &r) {
			dVec tempr;
			tempr = r;
			return (omega2*tempr/(2.0*constants()->lambda()));
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
 *
 * NB: Need to add an explanation to compare with the previous version of the
 *     Gaussian.
 */
class DeltaPotential : public PotentialBase  {
	public:
		DeltaPotential (double,double);
		~DeltaPotential ();

		/**
		 * Return the delta function potential with strength g approximated as
		 * the limit of a Gaussian distribution.
		 * Tested and working in Mathematica.
		 */
		double V(const dVec &r) {
			return (norm*exp(-dot(r,r)*inv2sigma2));
		}

		/**
		 * Return the gradient of the delta function potential with strength 
		 * g approximated as the limit of a Gaussian distribution.
		 * Tested and working in Mathematica.
		 */
		dVec gradV(const dVec &r) {
			return (-2.0*r*norm*inv2sigma2*exp(-dot(r,r)*inv2sigma2));
		}

	private:
		double norm;			// A normalization constant for fixed strength
        double inv2sigma2;      // 1/(2\sigma^2) where \sigma^2 is the variance
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
// SutherlandPotential Class
// ========================================================================  
/** 
 * Computes the potential energy for the periodic Sutherland model which
 * approximates long-range 1/r^2 interactions on a ring.
 * @see:  B. Sutherland, Phys Rev A 4, 2019 (1971) 
 */
class SutherlandPotential : public PotentialBase  {
	public:
		SutherlandPotential (double);
		~SutherlandPotential ();

		/**
		 * Return the Sutherland potential g/r^2.
		 */
		double V(const dVec &r) {
            double x = pioL*(sqrt(dot(r,r)) + EPS);
            /* if (x < EPS) */
            /*     return LBIG; */
            return g * pioL * pioL / (sin(x)*sin(x));
		}

		/**
		 * Return the gradient of the Sutherland potential.
		 */
		dVec gradV(const dVec &r) {
            double rnorm = sqrt(dot(r,r));
            double x = pioL*(rnorm + EPS);
            double s = sin(x);
            /* if (s < EPS) */
            /*     return 0.0; */
			return (-2.0* g * pioL * pioL * pioL * cos(x) / (s*s*s*rnorm)) * r;
		}

		/**
		 * Return the Laplacian of the Sutherland potential.
		 */
		double grad2V(const dVec &r) {
            double x = pioL*(sqrt(dot(r,r))+EPS);
            double s = sin(x);
			return 2.0* g * pioL * pioL * pioL * pioL * (2.0+cos(2*x)) / 
                (s*s*s*s);
		}

	private:
		double g;			    // The interaction constant g
		double pioL;		    // \pi/L
};

// ========================================================================  
// Dipolar Potential Class
// ========================================================================  
/** 
 * Computes the potential energy for polarized electrical dipoles with strength
 * D in reduced units where lengths are measured in units of a = m D / \hbar^2 and
 * energies in units of \hbar^2 / m a^2.
 * @see:  http://link.aps.org/doi/10.1103/PhysRevA.87.063604
 */
class DipolePotential : public PotentialBase  {
	public:
		DipolePotential ();
		~DipolePotential ();

		/**
		 * Return the dipole potential 1/r^3.
		 */
		double V(const dVec &r) {
            double x = sqrt(dot(r,r));
            if (x < EPS)
                return LBIG;
            return 1.0/(x*x*x);
		}

		/**
		 * Return the gradient of the dipole potential.
         * \partial (1/r^3) \partial r \hat{r} = -3r^{-5} \vec{r}
		 */
		dVec gradV(const dVec &r) {
            double x = sqrt(dot(r,r));
            if (x < EPS)
                return 0.0;
			return (-3.0/(x*x*x*x*x)) * r;
		}

		/**
		 * Return the Laplacian of the dipolar potential.
		 */
		double grad2V(const dVec &r) {
            double x = sqrt(dot(r,r));
            if (x < EPS)
                return 0.0;
			return 6.0/(x*x*x*x*x);
		}
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
				return LBIG;
			else
				return 0.0;
		}

		/** A delta function at rho=R. */
		dVec gradV(const dVec &r) {
			dVec tempr;
			tempr = r;
			tempr[2] = 0.0;
			if (abs((r[0]*r[0]+r[1]*r[1])-R*R)<1.0E-3)
				return LBIG*tempr;
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

        /* Laplacian of the LJ Wall potential */
        double grad2V(const dVec &);

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
        double valued2Vdr2 (const double);
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

/** 
 * Return the Laplacian of aziz potential for separation r using a 
 * lookup table. 
 */
inline double LJCylinderPotential::grad2V(const dVec &r) {
	double rnorm = sqrt(r[0]*r[0] + r[1]*r[1]);
	dVec tempr;
	tempr = r;
	tempr[2] = 0.0; // PBC in z-direction
	int k = static_cast<int>(rnorm/dR);
	double g2V;
	if (k >= tableLength)
		g2V = extd2Vdr2[1];
	else
		g2V = lookupd2Vdr2(k);
	return g2V;
}

// ========================================================================  
// LJ Hour Glass Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for an hour-glass shaped
 * cavity. 
 */
class LJHourGlassPotential : public PotentialBase { 

	public:
		LJHourGlassPotential (const double, const double, const double);
		~LJHourGlassPotential ();

		/** The integrated LJ hour glass potential. */
		double V(const dVec &); 

		/** The gradient of the potential. Use the prrimitive approx. */
		dVec gradV(const dVec &pos) {
			return (0.0*pos);
		}

        /* Laplacian of the LJ Wall potential. Use primitive approx. */
        double grad2V(const dVec & pos) {
            return 0.0;
        }

		/** Initial configuration corresponding to the LJ cylinder potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

	private:
		/* All the parameters needed for the LJ wall potential */
		double density;
		double sigma;
		double epsilon;

		double R;		// Radius of the tube at +/- L/2
		double dR;		// variation of radius

		/* double minV;	// The minimum value of the potential */

        double dz;      // The length discretization
        double L;       // The legnth of the pore
        double invd;    // The inverse of the variation length scale
        double R0;      // A constant used in the tanh radius function

        int M;          // How many slices we discretize in the z-direction

        /* Different functional forms for the hourglass radius */
        double Rlinear(double z) {
            return (2.0*dR*abs(z)/L + R - dR);
        }

        double Rtanh(double z) {
            double t = R0*tanh(z*invd);
            return (dR*t*t + R - dR);
        }
};

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

        /* The Laplacian of the Aziz potential */
		double grad2V(const dVec &);

	private:
		/* All the parameters of the Aziz potential */
		double rm, A, epsilon, alpha, D, C6, C8, C10;

		/* Used to construct the lookup tables */
		double valueV (const double);				
		double valuedVdr (const double);					
        double valued2Vdr2 (const double);

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
     
        /* The 2nd derivative of the F-function needed for the Aziz potential
         * Double checked with Mathematica --MTG */
		double d2F(const double x) {
			double ix = 1.0/x;
			double r = 2.0*D*ix*ix*ix*( 2.0*D*D*D*ix*ix*ix - 4.0*D*D*ix*ix 
                    - D*ix + 2.0) * exp(-(D*ix - 1.0)*(D*ix - 1.0));
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

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the Laplacian of aziz potential for separation r using a 
 * lookup table. 
 */

inline double AzizPotential::grad2V(const dVec &r) {
	double rnorm = sqrt(dot(r,r));
	double g2V;
	//g2V = (newtonGregory(lookupd2Vdr2,extd2Vdr2,rnorm)/rnorm)*r;
	//g2V = (direct(lookupd2Vdr2,extd2Vdr2,rnorm)/rnorm)*r;
	g2V = direct(lookupd2Vdr2,extd2Vdr2,rnorm);
	return g2V;
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
// Excluded Volume Class (volume excluded w/ large potential)
// ========================================================================  
/**
 * Computes potential energy for Gasparini potential.
 * This is the case of gradV = 0 for all positions in cell.
 */
class Gasparini_1_Potential : public PotentialBase {
	public:
		Gasparini_1_Potential (double, double, const Container*);
		~Gasparini_1_Potential ();

        /** The potential */
        double V(const dVec &r){ 
            double r2 = 0.0;
            if ((r[2] >= -excZ) && (r[2] <= excZ) && (r[1] >= -excY) && (r[1] <= excY))
                r2 = 1.0*V0;
            return r2;
		}

		/** The gradient of the potential. */
		dVec gradV(const dVec &) { return 0.0; }

        /** Laplacian of the potential. */
        double grad2V(const dVec &r) { return 0.0; }

    
        /** Initial configuration corresponding to FixedAziz potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /** get the exclusion lengths ay and az */
		Array<double,1> getExcLen(); 

        /* parameters needed for Gasp Potential_1 */
        const double excZ;      //half amt. of exclusion (z)
        const double excY;      //half amt. of exclusion (y)
        const double V0;        //scales the potential step
};

// ========================================================================  
// Hard Sphere Potential Class
// ========================================================================  
/** 
 * Computes the effective potential from the exact two-body density matrix
 * for hard spheres in 3D.  
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
            return ((sqrt(dot(r,r)) <= a) ? LBIG : 0.0);
        }

		/** The effective potential */
		double V(const dVec &, const dVec &);
        double dVdlambda(const dVec &, const dVec &);
        double dVdtau(const dVec &, const dVec &);

	private:
		double a;				// The strength of the delta function
};


// ========================================================================
// 1D Delta Potential Class
// ========================================================================
/**
 * Computes the effective potential from the exact two-body density matrix
 * for delta interactions in 1D
 *       
 */
class Delta1DPotential : public PotentialBase  {
public:
    Delta1DPotential (double);
    ~Delta1DPotential ();
    
    /** The classical potential */
    virtual double V(const dVec &r) {
        return (0.0);
    }
    
    /** The effective potential */
    double V(const dVec &, const dVec &);
    double dVdlambda(const dVec &, const dVec &);
    double dVdtau(const dVec &, const dVec &);
    
private:

    double erfCO;   //cutoff of erfc for asymptotic form
    double g;		// The strength of the delta function
    double l0;      // The COM kinetic length scale: 2\sqrt{\lambda\tau}
    double xi;      // The ratio l0/li
    double li;      // The COM interaction length scale: 2\lambda/g
    
    double xiSqOver2; // (xi^2)/2.0
    double xiSqrtPIOver2; // xi*\sqrt{\pi/2.0}
    
    /** Interaction weight and derivatives **/
    double Wint(double yt, double dxt);
    double dWdyt(double yt, double dxt);
    double dWdxi(double yt, double dxt);
    double dWddxt(double yt, double dxt);
};


// ========================================================================  
// Hard Rod Potential Class
// ========================================================================  
/** 
 * Computes the effective potential from the exact two-body density matrix
 * for hard rods in 1D.  
 *
 */
class HardRodPotential : public PotentialBase  {
	public:
		HardRodPotential (double);
		~HardRodPotential ();

		/** The classical potential */
		virtual double V(const dVec &r) { 
            return ((sqrt(dot(r,r)) <= a) ? LBIG : 0.0);
        }

		/** The effective potential */
		double V(const dVec &, const dVec &);
        double dVdlambda(const dVec &, const dVec &);
        double dVdtau(const dVec &, const dVec &);

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

// ========================================================================  
// GraphenePotential Class
// ========================================================================  
/** 
 * \brief Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GraphenePotential: public PotentialBase  {

	public:
		GraphenePotential(const double, const double, const double, const double, const double);
		~GraphenePotential();

		/* Return the sum of the van der Waals' interaction energy between the supplied
		 * particle and the fixed graphene lattice. */
		double V(const dVec &r);

		/** Initial configuration corresponding to graphene-helium vdW potential */
		Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 

	private:
		double sigma;
		double epsilon;

		double a1x;
		double a1y;
		double a2x;
		double a2y;

		double g1x;
		double g1y;
		double g2x;
		double g2y;

		double b1x;
		double b1y;
		double b2x;
		double b2y;

		double A;

};
#endif
