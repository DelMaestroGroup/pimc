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
        virtual dVec gradV(const dVec &) { return dVec{}; }

        /** Grad^2 of the potential*/
        virtual double grad2V(const dVec &) { return 0.0; }

        /** The derivative of the effective potential with respect to lambda
         *  and tau */
        virtual double dVdlambda(const dVec &, const dVec &) {return 0.0;}
        virtual double dVdtau(const dVec &, const dVec &) {return 0.0;}
        
        /** Default Initial configuration of particles*/
        virtual DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /** A debug method that output's the potential to a supplied separation */
        void output(const double);

        double tailV;       ///< Tail correction factor.

        /** Array to hold data elements*/
        virtual DynamicArray<double,1> getExcLen();

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
	DynamicArray <double,1> lookupV;           ///< A potential lookup table
	DynamicArray <double,1> lookupdVdr;        ///< A lookup table for dVint/dr
	DynamicArray <double,1> lookupd2Vdr2;      ///< A lookup table for d2Vint/dr2

        double dr;                          ///< The discretization for the lookup table
        int tableLength;                    ///< The number of elements in the lookup table

	std::array<double,2> extV;          ///< Extremal value of V
	std::array<double,2> extdVdr;       ///< Extremal value of dV/dr
	std::array<double,2> extd2Vdr2;     ///< Extremal value of d2V/dr2

        /* Initialize all data structures */
        void initLookupTable(const double, const double);

        /* Returns the 2-point spline fit to the lookup table */
        virtual double newtonGregory(const DynamicArray<double,1>&, const std::array<double,2>&, const double);

        /* Returns a bare lookup value */
        virtual double direct(const DynamicArray<double,1>&, const std::array<double,2>&, const double);

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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 
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
            dVec tempr{};
            for (int i=0; i < NDIM-1; i++)
                tempr[i] = r[i];
            return ( c * constants()->m() * w * w * tempr );
        }

    private:
        double w;               // The confining frequency
        double c;               // A dimension-full constant
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
        double norm;            // A normalization constant for fixed strength
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
        double c;               // The strength of the delta function
        double norm;            // A normalization constant for fixed strength
        double a;               // The order of the limit
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
            return g * pioL * pioL / (sin(x)*sin(x));
        }

        /**
         * Return the gradient of the Sutherland potential.
         */
        dVec gradV(const dVec &r) {
            double rnorm = sqrt(dot(r,r)) + EPS;
            double x = pioL*rnorm;
            double s = sin(x);
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
        double g;               // The interaction constant g
        double pioL;            // \pi/L
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


#if NDIM > 2
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
        double R;       // Radius of the tube
};
#endif

#if NDIM > 2
// ========================================================================  
// Plated LJ Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a plated cylindrical 
 * cavity. 
 */
class PlatedLJCylinderPotential : public PotentialBase, public TabulatedPotential {
    public:
        PlatedLJCylinderPotential (const double, const double, const double, const double, const double);
        ~PlatedLJCylinderPotential ();

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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        /** local methods for computing the potential of a LJ cylinder */
        double V_ (const double, const double, const double, const double, const double);
        double dVdr_ (const double, const double, const double, const double, const double);

        /* All the parameters needed for the LJ wall potential */
        double densityPlated;
        double sigmaPlated;
        double epsilonPlated;
        double density;
        double sigma;
        double epsilon;

        double Ri;      // Inner radius of the tube
        double Ro;      // Outer radius of the tube
        double dR;      // Discretization for the lookup table

        double minV;    // The minimum value of the potential

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
inline dVec PlatedLJCylinderPotential::gradV(const dVec &r) {
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
inline double PlatedLJCylinderPotential::grad2V(const dVec &r) {
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
#endif

#if NDIM > 2
// ========================================================================  
// LJ Cylinder Potential Class
// ========================================================================  
/** 
 * Computes the value of the external wall potential for a cylindrical 
 * cavity. 
 */
class LJCylinderPotential : public PotentialBase, public TabulatedPotential {
    public:
        LJCylinderPotential (const double, const double, const double, const double);
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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        /* All the parameters needed for the LJ wall potential */
        double density;
        double sigma;
        double epsilon;

        double R;       // Radius of the tube
        double dR;      // Discretization for the lookup table

        double minV;    // The minimum value of the potential

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
#endif

#if NDIM > 2
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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 
	DynamicArray<dVec,1> initialConfig1(const Container*, MTRand &, const int); 

    private:
        /* All the parameters needed for the LJ wall potential */
        double density;
        double sigma;
        double epsilon;

        double R;       // Radius of the tube at +/- L/2
        double dR;      // variation of radius

        /* double minV; // The minimum value of the potential */

        double L;       // The legnth of the pore
        double invd;    // The inverse of the variation length scale
        double R0;      // A constant used in the tanh radius function

        /* Different functional forms for the hourglass radius */
        double Rlinear(double z) {
            return (2.0*dR*abs(z)/L + R - dR);
        }

        double Rtanh(double z) {
            double t = R0*tanh(z*invd);
            return (dR*t*t + R - dR);
        }
};
#endif

// ========================================================================  
// Aziz Potential Class
// ========================================================================  
/** 
 * Computes the value of the semi-empircal Aziz potential that is known
 * to be accurate for He-4.
 */
class AzizPotential : public PotentialBase, public TabulatedPotential {
    public:
        AzizPotential (const Container *);
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
    g2V = direct(lookupd2Vdr2,extd2Vdr2,rnorm);
    return g2V;
}


// ========================================================================  
// Szalewicz Potential Class
// ========================================================================  
/** 
 * Computes the value of the semi-empircal Szalewicz potential that is known
 * to be accurate for He-4.
 */
class SzalewiczPotential : public PotentialBase, public TabulatedPotential {
    public:
        SzalewiczPotential (const Container *);
        ~SzalewiczPotential ();

        /* The Szalewicz HFDHE2 Potential */
        double V(const dVec &);

        /* The gradient of the Szalewicz potential */
        dVec gradV(const dVec &);

        /* The Laplacian of the Szalewicz potential */
        double grad2V(const dVec &);

    private:
        /* All the parameters of the Szalewicz potential */
        double rm = 2.9262186279335958;
        double Cn[17] = {   0.0,
                            0.0,
                            0.0,
                            0.000000577235,
                            -0.000035322,
                            0.000001377841,
                            1.461830,
                            0.0,
                            14.12350,
                            0.0,
                            183.7497,
                            -0.7674e2,
                            0.3372e4,
                            -0.3806e4,
                            0.8534e5,
                            -0.1707e6,
                            0.286e7
                            };
        double a = 3.64890303652830;
        double b = 2.36824871743591;
        double eta = 4.09423805117871;
        double Pn[3] = {-25.4701669416621, 269.244425630616, -56.3879970402079};
        double Qn[2] = {38.7957487310071, -2.76577136772754};
        long int factorials[17] = {     1,
                                        1,
                                        2,
                                        6,
                                        24,
                                        120,
                                        720,
                                        5040,
                                        40320,
                                        362880,
                                        3628800,
                                        39916800,
                                        479001600,
                                        6227020800,
                                        87178291200,
                                        1307674368000,
                                        20922789888000};
        /* Used to construct the lookup tables */
        double valueV (const double);               
        double valuedVdr (const double);                    
        double valued2Vdr2 (const double);

        /* The Tang-Toennies damping function needed for the Szalewicz potential */
        /* Can be described as upper incomplete gamma function.*/
        double fn(const double x, const int n) {
            double s1 = 0.0;
            for (int i = 0; i < n + 1; i++) {
            s1 += pow(x,i)/factorials[i];           
            }
            return 1.0 - (exp(-x)*s1);
        }

        /* The derivative of the Tang-Toennies damping function needed for the Szalewicz potential */
        double dfn(const double x, const int n) {
            return (exp(-x)*pow(x,n)/factorials[n]);
        }   
     
        /* The 2nd derivative of the Tang-Toennies needed for the Szalewicz potential*/
        double d2fn(const double x, const int n) {
            return (exp(-x)*(n-x)*pow(x,(n-1))/factorials[n]);
        }

        /* The Tang-Toennies damping function expanded for small r */
        /* needed for the Szalewicz potential */
        /* Can be described as upper incomplete gamma function. */
        double fn2(const double x, const int n) {
        
            double s1 = 0.0;
            for (int i = 0; i < 17; i++) {
                s1 += (pow(-1,i)*pow(x,(i+1)))/(factorials[i]*(n+i+1));         
            }
            s1 *= pow(x,n)/factorials[n];
            return s1;
        }

};

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// INLINE FUNCTION DEFINITIONS
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the Szalewicz potential for separation r using a lookup table. 
 */
inline double SzalewiczPotential::V(const dVec &r) {
    //double rnorm = sqrt(dot(r,r));
    //return newtonGregory(lookupV,extV,rnorm);
    return direct(lookupV,extV,sqrt(dot(r,r)));
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the gradient of Szalewicz potential for separation r using a 
 * lookup table. 
 */
inline dVec SzalewiczPotential::gradV(const dVec &r) {
    double rnorm = sqrt(dot(r,r));
    dVec gV;
    //gV = (newtonGregory(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
    gV = (direct(lookupdVdr,extdVdr,rnorm)/rnorm)*r;
    return gV;
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/** 
 * Return the Laplacian of Szalewicz potential for separation r using a 
 * lookup table. 
 */

inline double SzalewiczPotential::grad2V(const dVec &r) {
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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        AzizPotential aziz;                 // A copy of the aziz potential
	DynamicArray <dVec,1> fixedParticles;      // The location of the fixed particles
	DynamicArray <int,2> fixedBeadsInGrid;     // The local grid hash table
	DynamicArray <int,1> numFixedBeadsInGrid;  // The number of fixed particles in each grid box
        int numFixedParticles;              // The total number of fixed particles
        LookupTable *lookupPtr;             // A lookup table pointer
        double rc2;                         // A local copy of the potential cutoff squared

};

#if NDIM > 2
// ========================================================================  
// FixedPositionLJPotential Class
// ========================================================================  
/** 
 * @brief Returns Lennard-Jones potential between adatoms and fixed postions in FILENAME.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from interaction between adatom and fixed positions
 * given in FILENAME.
 */
class FixedPositionLJPotential: public PotentialBase  {

    public:
        FixedPositionLJPotential(const double, const double, const Container*);
        ~FixedPositionLJPotential();

        /* Return the sum of the Lennard-Jones potential between the supplied
         * particle and the fixed positions found in FILENAME. */
        double V(const dVec &r);

    private:
        const Container *boxPtr;
        double sigma;
        double epsilon;
        double Lz;
        
	DynamicArray <dVec,1> fixedParticles;      // The location of the fixed particles
        int numFixedParticles;              // The total number of fixed particles
};
#endif

#if NDIM > 2
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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

        /** get the exclusion lengths ay and az */
	DynamicArray<double,1> getExcLen(); 

        /* parameters needed for Gasp Potential_1 */
        const double excZ;      //half amt. of exclusion (z)
        const double excY;      //half amt. of exclusion (y)
        const double V0;        //scales the potential step
};
#endif

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
        double a;               // The strength of the delta function
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
    double g;       // The strength of the delta function
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
        double a;               // The strength of the delta function
};

// ========================================================================  
// Carbon Nanotube Potential Class
// ========================================================================  
/** 
 * The smooth non-corregated version of the helium-carbon nanotube potential.
 * @see http://prb.aps.org/abstract/PRB/v62/i3/p2173_1
 */
//class CarbonNanotubePotential : public PotentialBase, public TabulatedPotential {
//  public:
//      CarbonNanotubePotential(const double);
//      ~CarbonNanotubePotential();
//
//      /** The cylindrically symmetric potential. */
//      double V(const dVec &r) {
//          int k = int(sqrt(r[0]*r[0] + r[1]*r[1])/dR);
//          if (k >= tableLength)
//              return extV[1];
//          else
//              return lookupV(k);
//      }
//
//      /* The gradient of the CNT potential */
//      dVec gradV(const dVec &);
//
//      /** Initial configuration corresponding to the CNT potential */
//      Array<dVec,1> initialConfig(const Container*, MTRand &, const int); 
//
//  private:
//      /* All the parameters needed for the LJ wall potential */
//      double density;
//      double sigmaHe,sigmaC;
//      double epsilonHe,epsilonC;
//
//      double R;       // Radius of the tube
//      double dR;      // Discretization for the lookup table
//
//      double minV;    // The minimum value of the potential
//
//      /* Used to construct the lookup tables */
//      double valueV (const double);               
//      double valuedVdr (const double);                    
//};

#if NDIM > 2
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
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

    private:
        double sigma;
        double epsilon;
        double Lz;
        double Lzo2;

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

#if NDIM > 2
// ========================================================================  
// GrapheneLUTPotential Class
// ========================================================================  
/** 
 * \brief Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols & Adrian Del Maestro
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUTPotential: public PotentialBase  {

    public:
        GrapheneLUTPotential(const double, const double, const double, const double, const double, const Container*);
        ~GrapheneLUTPotential();

        /* Return the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        double V(const dVec &r);

        /* Return the gradient of the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        dVec gradV(const dVec &r);
        
        /** Initial configuration corresponding to graphene-helium vdW potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 

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

        double Lzo2;            ///< half the system size in the z-direction
        double Lz;              ///< The size of the system in the z-direction
        double zWall;           ///< The location of the onset of the "hard" wall 
        double invWallWidth;    ///< How fast the wall turns on.
        double V_zmin;          ///< A large potential value used for the cutoff
        
        double q = 2.0; // number of basis vectors
        
        static const int gnum = 3;
        static const int gtot = 16; //pow(gnum + 1,2);

        double dr = 1.0E-5;
        double zmin = 1.5;
        double zmax = 10.0;
        int tableLength;

        /* Array<int,2> karr; */
        /* double garr [gtot]; */

	DynamicArray<int,1> gMagID;    // g-magnitude lookup index
        
        /* The lookup tables */
	DynamicArray<double,2> vg;
	DynamicArray<double,2> gradvg;
};
#endif

#if NDIM > 2
// ========================================================================  
// GrapheneLUT3DPotential Class
// ========================================================================  
/** 
 * \brief Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUT3DPotential: public PotentialBase  {

    public:
        GrapheneLUT3DPotential(const std::string, const Container*);
        ~GrapheneLUT3DPotential();

        /* Return the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        double V(const dVec &);

        /* Return the gradient of the sum of the van der Waals' interaction energy between the supplied
         * particle and the fixed graphene lattice. */
        dVec gradV(const dVec &);
        double grad2V(const dVec &);
        
        /** Initial configuration corresponding to graphene-helium vdW potential */
	DynamicArray<dVec,1> initialConfig(const Container*, MTRand &, const int); 
	DynamicArray<dVec,1> initialConfig1(const Container*, MTRand &, const int); 

        void put_in_uc( dVec &, double, double);
        void cartesian_to_uc( dVec &, double, double, double, double);
        double trilinear_interpolation(DynamicArray<double,3>,dVec,double,double,double);
        double direct_lookup(DynamicArray<double,3>,dVec,double,double,double);

    private:
        double Lzo2;      ///< half the system size in the z-direction
        double zWall;     ///< The location of the onset of the "hard" wall 
        double invWallWidth; ///< How fast the wall turns on.
        
        /* spacing of the lookup tables */
        double dx;
        double dy;
        double dz;

        /* dimensions of the lookup tables */
        double cell_length_a;
        double cell_length_b;
        double zmin;
        double zmax;
        double V_zmin;
        
        /* transfer matrix */
        double A11;
        double A12;
        double A21;
        double A22;

	DynamicArray<double,3> V3d; // Potential lookup table
	DynamicArray<double,3> gradV3d_x; // gradient of potential x direction lookup table
	DynamicArray<double,3> gradV3d_y; // gradient of potential y direction lookup table
	DynamicArray<double,3> gradV3d_z; // gradient of potential z direction lookup table
	DynamicArray<double,3> grad2V3d; // Laplacian of potential
	DynamicArray<double,1> LUTinfo;

};

/****
 * Put inside unit cell
 ***/
inline void GrapheneLUT3DPotential::put_in_uc(dVec &r, double cell_length_a, double cell_length_b) {
    r[0] -= cell_length_a * floor(r[0]/cell_length_a);
    r[1] -= cell_length_b * floor(r[1]/cell_length_b);
}

/****
 * transfer from cartesian to unit cell coordinates
 ***/

inline void GrapheneLUT3DPotential::cartesian_to_uc( dVec &r, double A11, double A12, double A21, double A22) {
    double _x = A11 * r[0] + A12*r[1];
    double _y = A21 * r[0] + A22*r[1];
    r[0] = _x;
    r[1] = _y;
}

inline double GrapheneLUT3DPotential::trilinear_interpolation(DynamicArray<double,3> P,dVec r,double dx,double dy,double dz) {
    double x = r[0];
    double y = r[1];
    double z = r[2];

    double _xidx = floor(x/dx);
    double _yidx = floor(y/dy);
    double _zidx = floor(z/dz);
    int xidx = static_cast<int>(_xidx);
    int yidx = static_cast<int>(_yidx);
    int zidx = static_cast<int>(_zidx);

    double Delta_x = x/dx - _xidx;
    double Delta_y = y/dy - _yidx;
    double Delta_z = z/dz - _zidx;

    double c00 = P(xidx,yidx,zidx)*(1 - Delta_x) + P(xidx + 1,yidx,zidx)*Delta_x;
    double c01 = P(xidx,yidx,zidx + 1)*(1 - Delta_x) + P(xidx + 1,yidx,zidx + 1)*Delta_x;
    double c10 = P(xidx,yidx + 1,zidx)*(1 - Delta_x) + P(xidx + 1,yidx + 1,zidx)*Delta_x;
    double c11 = P(xidx,yidx + 1,zidx + 1)*(1 - Delta_x) + P(xidx + 1,yidx + 1,zidx + 1)*Delta_x;

    double c0 = c00*(1 - Delta_y) + c10*Delta_y;
    double c1 = c01*(1 - Delta_y) + c11*Delta_y;

    double c = c0*(1 - Delta_z) + c1*Delta_z;
    return c;
}

inline double GrapheneLUT3DPotential::direct_lookup(DynamicArray<double,3> P,dVec r,double dx,double dy,double dz) {
    double x = r[0];
    double y = r[1];
    double z = r[2];

    int xidx = static_cast<int>(x/dx);
    int yidx = static_cast<int>(y/dy);
    int zidx = static_cast<int>(z/dz);

    double c0 = P(xidx,yidx,zidx);
    return c0;
}
#endif

#if NDIM > 2
// ========================================================================  
// GrapheneLUT3DPotentialGenerate Class
// ========================================================================  
/** 
 * \brief FIXME Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal space.
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUT3DPotentialGenerate: public PotentialBase  {

    public:
        GrapheneLUT3DPotentialGenerate(
                const double, const double, const double, const double,
                const double, const int, const int, const int, const int, const Container*);
        ~GrapheneLUT3DPotentialGenerate();
        
    private:

        double Lzo2;    // half the system size in the z-direction
        
        /* resolution of the lookup tables */
        int xres;
        int yres;
        int zres;

        /* dimensions of the lookup tables */
        /* double zmin; */
        double zmax;
	int k_max;
        /* double V_zmin; */
        
        double Vz_64( double, double, double, int );
        double Vz_64( double, double, double );
        double Vz_64( double, double );
        
        double gradVz_x_64( double, double, double, int );
        double gradVz_x_64( double, double, double );
        double gradVz_x_64( double, double );
        
        double gradVz_y_64( double, double, double, int );
        double gradVz_y_64( double, double, double );
        double gradVz_y_64( double, double );
        
        
        double gradVz_z_64( double, double, double, int );
        double gradVz_z_64( double, double, double );
        double gradVz_z_64( double, double );

        double grad2Vz_64( double, double, double, int );
        double grad2Vz_64( double, double, double );
        double grad2Vz_64( double, double );
        
        double Vg_64( double, double, double, double, double, double, double,
                double, double, double );
        double Vg_64( double, double, double, double, double, double, double,
                double, double );
        
        double gradVg_x_64( double, double, double, double, double, double, double,
                double, double, double );
        double gradVg_x_64( double, double, double, double, double, double, double,
                double, double );

        double gradVg_y_64( double, double, double, double, double, double, double,
                double, double, double );
        double gradVg_y_64( double, double, double, double, double, double, double,
                double, double );

        double gradVg_z_64( double, double, double, double, double, double, double,
                double, double, double );
        double gradVg_z_64( double, double, double, double, double, double, double,
                double, double );

        double grad2Vg_64( double, double, double, double, double, double, double,
                double, double, double );
        double grad2Vg_64( double, double, double, double, double, double, double,
                double, double );

        double V_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double gradV_x_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double gradV_y_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double gradV_z_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );

        double grad2V_64( double, double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>, DynamicArray<int,1>,
                DynamicArray<int,1>, DynamicArray<double,1> );
        
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors();
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors( double );
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors( double, double, double );
        std::tuple<
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>,
            std::array<double,2>, std::array<double,2>
                > get_graphene_vectors_old( double, double, double );
        std::tuple< DynamicArray<int,1>, DynamicArray<int,1>, DynamicArray<double,1>
            > get_g_magnitudes( std::array<double,2>, std::array<double,2> );

        template <class T> double calculate_magnitude( T vec ) {
            return sqrt(dot(vec,vec));
        }
        
        template <class T1, class T2> double calculate_angle(
                T1 x, T2 y, double magnitude_x, double magnitude_y ) {
            return acos(dot(x,y)/magnitude_x/magnitude_y);
        }
        
        template <class T1, class T2> double calculate_angle( T1 x, T2 y ) {
            return calculate_angle(x,y,calculate_magnitude(x),calculate_magnitude(y));
        }

        void calculate_V3D_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_gradV3D_x_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_gradV3D_y_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_gradV3D_z_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        void calculate_grad2V3D_64(
                DynamicArray<double,3>, DynamicArray<double,2>, DynamicArray<double,2>,
                DynamicArray<double,1>, double, double,
                double, std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>,
                DynamicArray<double,1> );

        std::pair<double, double> get_z_min_V_min( 
                double, double, double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>, DynamicArray<double,1> );

        std::pair<double, double> get_z_V_to_find( 
                double, double, double,
                std::array<double,2>, std::array<double,2>,
                std::array<double,2>, std::array<double,2>,
                DynamicArray<int,1>, DynamicArray<int,1>, DynamicArray<double,1> );

	DynamicArray<double,3> get_V3D(
                double, double, double, int, int, int, double, double ); 
        std::pair<DynamicArray<double,3> , DynamicArray<double,1>> get_V3D(
                double, double, double, int, int, int, double );
        std::tuple<
            DynamicArray<double,3>, DynamicArray<double,3>, DynamicArray<double,3>,
            DynamicArray<double,3>, DynamicArray<double,3>, DynamicArray<double,2>,
            DynamicArray<double,2>, DynamicArray<double,1>
            > get_V3D_all( double, double, double, int, int, int, double );
};
#endif

#if NDIM > 2
// ========================================================================
// GrapheneLUT3DPotentialToBinary Class
// ========================================================================
/**
 * \brief FIXME Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal sp
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUT3DPotentialToBinary: public PotentialBase  {

    public:
        GrapheneLUT3DPotentialToBinary(const std::string, const Container*);
        ~GrapheneLUT3DPotentialToBinary();

    private:
	DynamicArray<double,3> V3d; // Potential lookup table
	DynamicArray<double,3> gradV3d_x; // gradient of potential x direction lookup table
	DynamicArray<double,3> gradV3d_y; // gradient of potential y direction lookup table
	DynamicArray<double,3> gradV3d_z; // gradient of potential z direction lookup table
	DynamicArray<double,3> grad2V3d; // Laplacian of potential
	DynamicArray<double,1> LUTinfo;

};
#endif

#if NDIM > 2
// ========================================================================
// GrapheneLUT3DPotentialToText Class
// ========================================================================
/**
 * \brief FIXME Returns van der Waals' potential between a helium adatom and a graphene sheet using summation in reciprocal sp
 *
 * Author: Nathan Nichols
 * Returns the potential energy resulting from a van der Waals' interaction
 * between a helium adatom and a fixed infinite graphene lattice
 */
class GrapheneLUT3DPotentialToText: public PotentialBase  {

    public:
        GrapheneLUT3DPotentialToText(const std::string, const Container*);
        ~GrapheneLUT3DPotentialToText();

    private:
	DynamicArray<double,3> V3d; // Potential lookup table
	DynamicArray<double,3> gradV3d_x; // gradient of potential x direction lookup table
	DynamicArray<double,3> gradV3d_y; // gradient of potential y direction lookup table
	DynamicArray<double,3> gradV3d_z; // gradient of potential z direction lookup table
	DynamicArray<double,3> grad2V3d; // Laplacian of potential
	DynamicArray<double,1> LUTinfo;

};
#endif

#endif
